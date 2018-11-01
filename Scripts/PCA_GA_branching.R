### This script is intended to visualize traits differences between SAM lines in a PCA
### but also compare within branching groups. SO SPECIFY BRANCHING AT INCDICATED (***) LOCATIONS!!!
### Also when switching between branching datasets, reread the "ind" file!


### Set working directory
setwd(dir = "/volumes/pbio-burke/users/ijameel/inam_SELECTION_239/")

library(ggplot2)
library(PerformanceAnalytics)
library(dplyr)
library(ade4)

### Read in your trait data
traits <- read.csv("Analysis_GA/Data/239_trait_means_102618.csv", header = TRUE)
traits <- traits[complete.cases(traits),] #remove individuals with ANY missing data

# make the subset branching datasets. specify the select fuction to include any rows you like 
branched_t <- subset(traits, Branching > 0, select=Lineage:X._Oil_.g.kg.)
unbranched_t <- subset(traits, Branching == 0, select=Lineage:X._Oil_.g.kg.)
All_t <- subset(traits, select=Lineage:X._Oil_.g.kg.)

### PCA stuff
ind <- read.csv(file = "239_PCA_102318.csv", header = TRUE)
# remove lines in the ind data that are not in the data set. 
# specify dataset! ***
ind <- subset(ind,Plant %in% 
                #All_t
                unbranched_t
                #branched_t
              $Plant)

# specify the category from the ind file that you want to color the dots according to  
ggplot(ind, aes(x=PC1, y=PC2, color=Oil)) +
  geom_point(shape=1)    # Use hollow circles

# write csv if needed
write.csv(ind, "PCA_Unbranched_GA2010.csv")
