rm(list = ls())
cat("\f")
##Sleep disorder 
setwd("C:/Users/rsh001/OneDrive - University of Bergen/Project files/Rajesh/Sleep_DATA")
library(RColorBrewer)
library(readr)
library(openxlsx)
library(tidyverse)
library(microbiome)
library(vegan) # adonis
library(RVAideMemoire) # Pairwise permanova
library(compositions)
library(magrittr)
library(qwraps2)
library(dplyr)
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(ANCOMBC)
library(DT)
library(microbiomeutilities)
library(readxl)
library(lubridate)
library(haven)
library(mia)
library(kableExtra)
library(DT)
library(openxlsx)
library(gtsummary)


# read the metadata 

sleep_metadata = read.xlsx("Rh_sleep.xlsx",
                           sep=",", colNames = TRUE)



#calculate age by using the date of birth and date of participating 

sleep_metadata <- sleep_metadata  %>% 
  mutate(date_parti= make_datetime(qyearo, qmontho,qdayo))

sleep_metadata <- sleep_metadata  %>% 
  mutate(date_birth= make_datetime(byearo, bmontho,bdayo))

sleep_metadata$N_age = as.numeric(difftime(sleep_metadata$date_parti, sleep_metadata$date_birth, units = "weeks"))/52.25


selected<- select(sleep_metadata, "#SampleID", N_age, 
                  adult,  sexo, heighto, weighto, gero,
                  hoursleepo,emao, dmso, diso, "osaso.y","osaso.x" , educationo, antibiotics, "snoringo.x","snoringo.y","smoke",minutessleepo,)
ls (selected) 
table(selected$snoringo.x)

#Heart burn and belching

selected$gero[selected$gero == 99] <- NA
selected$gero[selected$gero==1] <- "never or almost never"
selected$gero[selected$gero ==2]<- "less than once a week"
selected$gero[selected$gero == 3]<- "once or twice a week"
selected$gero[selected$gero == 4]<- "3-5 nights/days a week" # merge lower two group
selected$gero[selected$gero == 5]<- "3-5 nights/days a week"



selected$gero = factor(selected$gero, c( "never or almost never", "less than once a week", "once or twice a week", "3-5 nights/days a week"))


#Difficulty getting to sleep at night
selected$diso[selected$diso == 99] <- NA
selected$diso[selected$diso ==1] <- "never or almost never"
selected$diso[selected$diso ==2]<- "less than once a week"
selected$diso[selected$diso == 3]<- "once or twice a week"
selected$diso[selected$diso == 4]<- "3-5 nights/days a week"
selected$diso[selected$diso == 5]<- "3-5 nights/days a week"

selected$diso = factor(selected$diso, c( "never or almost never", "less than once a week", "once or twice a week", "3-5 nights/days a week"))

selected$snoringo.x[selected$snoringo.x == 99] <- NA
selected$snoringo.x[selected$snoringo.x ==1] <- "never or almost never"
selected$snoringo.x[selected$snoringo.x ==2]<- "less than once a week"
selected$snoringo.x[selected$snoringo.x == 3]<- "once or twice a week"
selected$snoringo.x[selected$snoringo.x == 4]<- "3-5 nights/days a week"
selected$snoringo.x[selected$snoringo.x == 5]<- "3-5 nights/days a week"

selected$snoringo.x = factor(selected$snoringo.x, c( "never or almost never", "less than once a week", "once or twice a week", "3-5 nights/days a week"))

selected = selected%>%mutate(bmi = as.numeric(weighto) / (as.numeric(heighto) / 100)^2, age = as.numeric(N_age), gender =sexo , 
                             height =heighto, weight=weighto, snoring=snoringo.x, sleep_difficutly = diso, heart_burn =gero )




table (selected$sleep_difficutly)##Difficulty getting to sleep at night

table (selected$snoring)#snoring 

table (selected $heart_burn)##Heart burn and belching


#OTU table 

otu_data <-  read_tsv("feature-table.txt", skip = 1)

ls (otu_data)

##remove the column otu since it is now used as a row name
otu_mat <- otu_data %>% select (-"#OTU ID")

##Taxonomy table
tax <- read.table(file = 'taxonomy.tsv', sep = '\t', header = TRUE)
taxtable<-tax %>% as.tibble() %>% separate(Taxon, sep="; ", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxtable


##  Idem for the two other matrices

row.names(taxtable) <- taxtable$Feature.ID
tax_mat <- taxtable %>% select (-Feature.ID)

### tree
##phylo-class (tree)   # i cant connect the tree 

library("ape")
tree <- read.tree("tree.nwk",text = NULL, tree.names = NULL, skip = 0,comment.char = "", keep.multi = FALSE)



##Transform into matrixes otu and tax tables (sample table can be left as data frame)

otu_mat <- as.matrix(otu_mat)

tax_mat <- as.matrix(tax_mat)

dim(tax_mat)
se <- SummarizedExperiment(assays = list(counts = otu_mat),
                           colData = selected,
                           rowData = tax_mat)

se
assays(se)$counts[1:3, 1:3]
names(rowData(se)) <- c("Kingdom", "Phylum", "Class", "Order", 
                        "Family", "Genus")
head(rowData(se))


pseq <- makePhyloseqFromTreeSummarizedExperiment(se)


# Exclude all participants who have used antibiotics in the last four weeks
pseq = subset_samples(pseq, antibiotics != "last4weeks")

# keep only adults 
otu_data_a = subset_samples(pseq, adult == "1")
otu_data_a

#Removing missing 
otu_data_a = subset_samples(otu_data_a, snoring != "na")

#extrating metadata
m <- meta(otu_data_a)
table(m$snoring)
ls(m)


#===============================================================================
# snoring 
#===============================================================================
# Running differential abundance analysis using ANCOM-BC, at the Genus level

# Aggregate to Genus level
Genus_data = aggregate_taxa(otu_data_a, "Genus")

tax_mat = as(tax_table(Genus_data), "matrix")
# Run ANCOM-BC
feature_table = abundances(Genus_data); meta_data = meta(Genus_data)
# ANCOM-BC requires an id column for metadata
out = ancombc(phyloseq = Genus_data, formula = "snoring + bmi+ age + gender+ smoke + educationo",
              p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
              group = "snoring", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
out


# Outputs
res = out$res
res_global = out$res_global
res_global

res_direct = out$res_direct
res_direct

res_pattern = out$res_pattern
res_pattern



tab_coef = res$beta
tab_coef
col_name = c("less than once a week", "once or twice a week", "3-5 nights/days a week",  "weight", "smoke",
             "smokeprevious", "educationo2")
colnames(tab_coef) = col_name
tab_coef %>% datatable(caption = "Coefficients from the Primary Result") %>%
  formatRound(col_name, digits = 2)




# SEs
tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name, digits = 2)

#Test statistics
tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name, digits = 2)


#P-values

tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name, digits = 2)

# Adjusted p-values
tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name, digits = 2)

#Differentially abundant taxa

tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>%
  datatable(caption = "Differentially Abundant Taxa
            from the Primary Result")



# Waterfall plot
tab_coef = res$beta
col_name = c("less than once a week", "once or twice a week", "3-5 nights/days a week",  "weight", "smoke",
             "smokeprevious", "educationo2")
colnames(tab_coef) = col_name

tab_se = res$se
colnames(tab_se) = col_name

tab_w = res$W
colnames(tab_w) = col_name

tab_p = res$p_val
colnames(tab_p) = col_name

tab_q = res$q
colnames(tab_q) = col_name

tab_diff = res$diff_abn
colnames(tab_diff) = col_name

tab_zero = out$zero_ind
tab_zero = tab_zero - tab_zero[, 1]
tab_zero = abs(tab_zero[, -1])
tab_zero = data.frame(tab_zero, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")


# Coerce the SE of structural zero to be zero
tab_se[, grepl("snoring", colnames(tab_se))] = tab_se[, grepl("snoring", colnames(tab_se))] *
  (1 - tab_zero[, grepl("snoring", colnames(tab_zero))])

# less than once a week
df1 = data.frame(genus = rownames(tab_coef),
                 lfc = tab_coef$`less than once a week`,
                 se = tab_se$`less than once a week`,
                 p = tab_p$`less than once a week`,
                 q = tab_q$`less than once a week`) %>%
  filter(q < 0.05) %>%
  arrange(desc(lfc)) %>%
  mutate(type = if_else(lfc > 0, "g1", "g2"),
         star = case_when(p < .001 ~ "***",
                          p < .01 ~ "**",
                          TRUE ~ "*"),
         pos = if_else(type == "g1", 
                       lfc + se + 0.2,
                       lfc - se - 0.2)
  )

df1$genus = factor(df1$genus, levels = df1$genus)
p1 = df1 %>%
  ggplot(aes(x = genus, y = lfc, 
             fill = type, color = type)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc - se, 
                    ymax = lfc + se), 
                width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  geom_text(aes(y = pos, label = star), 
            vjust = .7, color = "black", 
            position = position_dodge(width = 0.05)) +
  labs(x = NULL, y = "Log fold change") +
  guides(color = FALSE, drop = FALSE) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 10, face = "bold",
                                   angle = 60, hjust = 1)) + 
  labs(title = "less than once a week") +
  scale_fill_brewer(palette = "Set1", drop = FALSE,
                    name = NULL,
                    label = c("g1" = "Increased abundance",
                              "g2" = "Reduced abundance")) +
  scale_color_brewer(palette = "Set1", drop = FALSE) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 

ggsave(plot = p1, filename = "less than once a week.jpeg", height = 5, width = 6.25, dpi = 300)



# once or twice a week


df2 = data.frame(genus = rownames(tab_coef),
                 lfc = tab_coef$`once or twice a week`,
                 se = tab_se$`once or twice a week`,
                 p = tab_p$`once or twice a week`,
                 q = tab_q$`once or twice a week`) %>%
  filter(q < 0.05) %>%
  arrange(desc(lfc)) %>%
  mutate(type = if_else(lfc > 0, "g1", "g2"),
         star = case_when(p < .001 ~ "***",
                          p < .01 ~ "**",
                          TRUE ~ "*"),
         pos = if_else(type == "g1", 
                       lfc + se + 0.2,
                       lfc - se - 0.2)
  )

df2$genus = factor(df2$genus, levels = df2$genus)
p2 = df2 %>%
  ggplot(aes(x = genus, y = lfc, 
             fill = type, color = type)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc - se, 
                    ymax = lfc + se), 
                width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  geom_text(aes(y = pos, label = star), 
            vjust = .7, color = "black", 
            position = position_dodge(width = 0.05)) +
  labs(x = NULL, y = "Log fold change") +
  guides(color = FALSE, drop = FALSE) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 10, face = "bold",
                                   angle = 60, hjust = 1)) + 
  labs(title = "once or twice a week") +
  scale_fill_brewer(palette = "Set1", drop = FALSE,
                    name = NULL,
                    label = c("g1" = "Increased abundance",
                              "g2" = "Reduced abundance")) +
  scale_color_brewer(palette = "Set1", drop = FALSE) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 

ggsave(plot = p2, filename = "once or twice a week.jpeg", height = 5, width = 6.25, dpi = 300)





#3-5 nights/days a week



df3 = data.frame(genus = rownames(tab_coef),
                 lfc = tab_coef$`3-5 nights/days a week`,
                 se = tab_se$`3-5 nights/days a week`,
                 p = tab_p$`3-5 nights/days a week`,
                 q = tab_q$`3-5 nights/days a week`) %>%
  filter(q < 0.05) %>%
  arrange(desc(lfc)) %>%
  mutate(type = if_else(lfc > 0, "g1", "g2"),
         star = case_when(p < .001 ~ "***",
                          p < .01 ~ "**",
                          TRUE ~ "*"),
         pos = if_else(type == "g1", 
                       lfc + se + 0.2,
                       lfc - se - 0.2)
  )

df3$genus = factor(df3$genus, levels = df3$genus)
p3 = df3 %>%
  ggplot(aes(x = genus, y = lfc, 
             fill = type, color = type)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc - se, 
                    ymax = lfc + se), 
                width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  geom_text(aes(y = pos, label = star), 
            vjust = .7, color = "black", 
            position = position_dodge(width = 0.05)) +
  labs(x = NULL, y = "Log fold change") +
  guides(color = FALSE, drop = FALSE) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 10, face = "bold",
                                   angle = 60, hjust = 1)) + 
  labs(title = "3-5 nights/days a week") +
  scale_fill_brewer(palette = "Set1", drop = FALSE,
                    name = NULL,
                    label = c("g1" = "Increased abundance",
                              "g2" = "Reduced abundance")) +
  scale_color_brewer(palette = "Set1", drop = FALSE) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 
p3

ggsave(plot = p3, filename = "3-5 nights/days a week.jpeg", height = 3, width = 2, dpi = 300) # cant save this figure 

