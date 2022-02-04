


#Alpha and beta diversity 
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
head(rowData(se))


pseq <- makePhyloseqFromTreeSummarizedExperiment(se)


# Exclude all participants who have used antibiotics in the last four weeks
pseq = subset_samples(pseq, antibiotics != "last4weeks")

# Split by age
otu_data_a = subset_samples(pseq, adult == "1")
otu_data_a


#Removing missing 
otu_data_a = subset_samples(otu_data_a, snoring != "na")


# Aggregate to Genus level
Genus_data = aggregate_taxa(otu_data_a, "Genus")

meta_data_a = meta(otu_data_a)



ls(meta_data_a)




#Rerifraction 
ps.rarefied = rarefy_even_depth(Genus_data, rngseed=1, sample.size=1103, replace=F)


#

hmp.div <- alpha(ps.rarefied, index = "all")

hmp.div

datatable(hmp.div)

# get the metadata out as seprate object
hmp.meta <- meta(ps.rarefied)

# Add the rownames as a new colum for easy integration later.
hmp.meta$sam_name <- rownames(hmp.meta)

# Add the rownames to diversity table
hmp.div$sam_name <- rownames(hmp.div)

# merge these two data frames into one
div.df <- merge(hmp.div,hmp.meta, by = "sam_name")

# check the tables
colnames(div.df)


# convert phyloseq object into a long data format.

div.df2 <- div.df[, c("snoring",  "diversity_shannon", "observed" , "chao1" , "evenness_pielou" )]

# the names are not pretty. we can replace them

colnames(div.df2) <- c("snoring",  "Shannon", "Observed"  , "Chao1" ,"Pielous evenness")

# check
colnames(div.df2)


div_df_melt <- reshape2::melt(div.df2)
head(div_df_melt)



# Now use this data frame to plot
p <- ggboxplot(div_df_melt, x = "snoring", y = "value",
               fill = "snoring",
               palette = "lancet",
               legend= "right",
               facet.by = "variable",
               scales = "free")

p <- p + rotate_x_text()
# we will remove the x axis lables

p <- p + rremove("x.text")
p




lev <- levels(div_df_melt$snoring) # get the variables
lev
# make a pairwise list that we want to compare.
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])


pval <- list(
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
  symbols = c("****", "***", "**", "*", "0.05")
)

p2 <- p + stat_compare_means(
  comparisons = L.pairs,
  label = "p.signif",
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("****", "***", "**", "*", "0.05")
  )
)
p2



ggsave("alphadiversity_snoring.png", height = 10, width = 10, units = 'in', dpi = 300)


## snoring beta diveristy

epseq = subset_samples(ps.rarefied, !is.na(get("snoring")))
set.seed(123)

# PERMANOVA
permanova = adonis(t(abundances(epseq)) ~ snoring,
                   data = meta(epseq),
                   permutations = 999, method = "bray")$aov.tab

# PERMDISP
dis = vegdist(t(abundances(epseq)), method = "bray")
groups = meta(epseq)[, "snoring"]
mod = betadisper(d = dis, group = groups, type = "median")

# Draw the Plot
p1 = signif(permanova$`Pr(>F)`[1], 2)
p2 = signif(permutest(mod)$tab$`Pr(>F)`[1], 2)

labs = paste0("PCoA", 1:2, " (", signif(100 * mod$eig / sum(mod$eig), 3), "%)")

# Export
# Consider change the palette of "#FF99FF", "#FF0000FF", "#CCFF00"
# since it is hard to read
# but make sure the palette of beta diversity plot is the same alpha diversity plot
jpeg(filename = "betadiversity_snoring.jpeg",
     height = 5, width = 8, res = 300, units = "in")
plot(mod, pch = 15:18, cex.lab = 1.25, cex = 0.7, label.cex = 0.6,
     main = "PCoA plot of Bray-Curtis dissimilarity", xlab = labs[1], ylab = labs[2],
     ylim = c(-0.5, 0.5), xlim = c(-0.6, 0.6), xaxt = "n",
     col = c("blue", "#FF0000FF", "black", "green"), sub = NULL,
     hull = FALSE, ellipse = TRUE, conf = 0.68) # 68% data coverage for data ellipses
axis(1, at = round(seq(-0.6, 0.6, by = 0.2), 1), las = 1)
legend(0.5, 0.4, legend = c ("never", "< 1/week", "1-2/week", "3-5 nights/days a week"),
       
       col = c("blue", "#FF0000FF", "black", "green"),
       pch = 15:16, cex = 0.8)
legend(x = 0.2, y = -0.3, cex = 0.7,
       legend = c(paste0("p (PERMANOVA) = ", p1),
                  paste0("p (PERMDISP) = ", p2)))
dev.off()


