
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
library("phyloseq")
library(tidyverse)
library(ggpubr)
library(ANCOMBC)
library(DT)
library(microbiomeutilities)
library("readxl")
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

#Heart burn and belching #asked by using questionnaires 

selected$gero[selected$gero == 99] <- NA
selected$gero[selected$gero==1] <- "never or almost never"
selected$gero[selected$gero ==2]<- "less than once a week"
selected$gero[selected$gero == 3]<- "once or twice a week"
selected$gero[selected$gero == 4]<- "3-5 nights/days a week"  #merged lower two group 
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


table (selected$diso)##Difficulty getting to sleep at night


table (selected$snoringo.y)#snoring 

table (selected $gero)##Heart burn and belching


#OTU table 

otu_data <-  read_tsv("feature-table.txt", skip = 1)

ls (otu_data)

##remove the column otu since it is now used as a row name
otu_mat <- otu_data %>% select (-"#OTU ID")

##Taxonomy table
tax <- read.table(file = 'taxonomy.tsv', sep = '\t', header = TRUE)
taxtable<-tax %>% as.tibble() %>% separate(Taxon, sep="; ", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxtable


##  Idem for the two other matrixes

row.names(taxtable) <- taxtable$Feature.ID
tax_mat <- taxtable %>% select (-Feature.ID)

### tree
##phylo-class (tree)

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

# Exclude all participants who have used antibiotics in the last four weeks
otu_data = subset_samples(otu_data, antibiotics != "last4weeks")

# Split by age
otu_data_a = subset_samples(otu_data, adult == "1")
otu_data_a

meta_data_a = meta(otu_data_a)
ls(meta_data_a)


#renaming the variables 
meta_data_a = meta_data_a%>%mutate(bmi = as.numeric(weighto) / (as.numeric(heighto) / 100)^2, age = as.numeric(N_age), gender =sexo , 
                                   height =heighto, weight=weighto, snoring=snoringo.x, sleep_difficutly = diso, heart_burn =gero )



#Descriptive table for total sample 
meta_data_a %>%
  select( "age",  "gender",  "height", "weight", bmi, "smoke" , "educationo", snoring, heart_burn, sleep_difficutly) %>%
  tbl_summary(
    type = all_continuous() ~ "continuous2",
    statistic = all_continuous() ~ c("{N_nonmiss}",
                                     "{mean} ({sd})",
                                     "{min}, {max}"),
    missing = "no"
  )

#stratified based on snoring
meta_data_a %>%
  select("age",  "gender",  "height", "weight", bmi, "smoke" , "educationo",snoring, heart_burn, sleep_difficutly) %>%
  tbl_summary(
    by = "snoring",
    type = all_continuous() ~ "continuous2",
    statistic = all_continuous() ~ c("{N_nonmiss}",
                                     "{mean} ({sd})",
                                     "{min}, {max}"),
    missing = "no"
  ) %>%
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2))


table(meta_data_a$snoring, useNA = "always")









