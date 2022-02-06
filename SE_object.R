
rm(list = ls())

## clear the Console

cat("\f")

setwd("C:/Users/rsh001/OneDrive - University of Bergen/Project files/Rajesh/Sleep_DATA")

#Install the packages 

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


#Select the required variables 

selected<- select(sleep_metadata, "#SampleID", N_age, 
                  adult,  sexo, heighto, weighto, gero,
                  hoursleepo,emao, dmso, diso, "osaso.y","osaso.x" , educationo, antibiotics, "snoringo.x","snoringo.y","smoke",minutessleepo,)

#Heart burn and belching
selected$gero[selected$gero == 99] <- NA
selected$gero[selected$gero==1] <- "never or almost never"
selected$gero[selected$gero ==2]<- "less than once a week"
selected$gero[selected$gero == 3]<- "once or twice a week"
selected$gero[selected$gero == 4]<- "3-5 nights/days a week" # merge lower two group
selected$gero[selected$gero == 5]<- "3-5 nights/days a week"


#Difficulty getting to sleep at night
selected$diso[selected$diso == 99] <- NA
selected$diso[selected$diso ==1] <- "never or almost never"
selected$diso[selected$diso ==2]<- "less than once a week"
selected$diso[selected$diso == 3]<- "once or twice a week"
selected$diso[selected$diso == 4]<- "3-5 nights/days a week"
selected$diso[selected$diso == 5]<- "3-5 nights/days a week"

#snoring
selected$snoringo.x[selected$snoringo.x == 99] <- NA
selected$snoringo.x[selected$snoringo.x ==1] <- "never or almost never"
selected$snoringo.x[selected$snoringo.x ==2]<- "less than once a week"
selected$snoringo.x[selected$snoringo.x == 3]<- "once or twice a week"
selected$snoringo.x[selected$snoringo.x == 4]<- "3-5 nights/days a week"
selected$snoringo.x[selected$snoringo.x == 5]<- "3-5 nights/days a week"


#Reordering the categories of the variables

selected$snoringo.x = factor(selected$snoringo.x, c( "never or almost never", "less than once a week", "once or twice a week", "3-5 nights/days a week"))

selected$gero = factor(selected$gero, c( "never or almost never", "less than once a week", "once or twice a week", "3-5 nights/days a week"))



#Renaming the variables 

selected = selected%>%mutate(bmi = as.numeric(weighto) / (as.numeric(heighto) / 100)^2, age = as.numeric(N_age), gender =sexo , 
                                   height =heighto, weight=weighto, snoring=snoringo.x, sleep_difficutly = diso, heart_burn =gero )


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

##Transform into matrixes otu and tax tables (sample table can be left as data frame)

otu_mat <- as.matrix(otu_mat)

tax_mat <- as.matrix(tax_mat)

dim(tax_mat)
se <- SummarizedExperiment(assays = list(counts = otu_mat),
                           colData = selected,
                           rowData = tax_mat)

names(rowData(se)) <- c("Kingdom", "Phylum", "Class", "Order", 
                        "Family", "Genus")
head(rowData(se))

#Tree file (i cant merge it to make tse)

library("ape")
tree <- read.tree("tree.nwk",text = NULL, tree.names = NULL, skip = 0,comment.char = "", keep.multi = FALSE)

tse <- as(se, "TreeSummarizedExperiment")

# Add tree to rowTree
rowTree(tse) <- tree
# Check
tse































