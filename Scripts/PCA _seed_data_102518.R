### This script is intended to visualize traits differences between SAM lines in a PCA
### but also compare within branching groups. SO SPECIFY BRANCHING AT INCDICATED (***) LOCATIONS!!!
### Also when switching between branching datasets, reread the "ind" file!


### Set working directory
setwd(dir = "/volumes/pbio-burke/users/ijameel/inam_SELECTION_239/")

library(ggplot2)
library(PerformanceAnalytics)
library(dplyr)
library(ade4)

### read in individual data with information for heterotic groups/oil cont etc and the trait data
ind <- read.csv(file = "239_PCA_102318.csv", header = TRUE)
traits <- read.csv("Analysis_Seeds/Data/239_shape_mean_102318.csv", header = TRUE)

# remove lines with any missing data, we will do this for all data sets
traits <- traits[complete.cases(traits),] 
str(traits)

# using subset function, make two data sets for branched plants and unbranched plants
branched_t <- subset(traits, Branching > 0, select=Line:AverageChroma)
unbranched_t <- subset(traits, Branching == 0, select=Line:AverageChroma)
All_t <- subset(traits, select=Line:AverageChroma)


# if you are just interested in plotting the genetic PCAs, skip the tpca step and go to the 
# ind <- subset step, and skip the binding columns and straight to plotting

# run the pca and set 4 for axes, remove traits that are not comparable between lines and the branching
# Specify the data set you would like ***!

tpca <- dudi.pca(
  All_t
  #unbranched_t
  #branched_t
  [,-c(1:3)])
  4

# make data frames for the "li" (row coordinates) and "co" (principal components). 
# refer to the dudi.pca help for more info 
li <- tpca$li
co <- tpca$co

# remove lines in the ind data that are not in the data set. 
# specify dataset! ***
ind <- subset(ind,Plant %in% 
                #All_t
                unbranched_t
                #branched_t
              $Plant)

# add the li data frame to the individual data
ind <- bind_cols(ind,li)

# check to make sure the axes you specified are added to the end of the file
tail (ind)


#############################################################################################################################
#############################################################################################################################
########################################### PCA and box&whisker #############################################################
#############################################################################################################################
#############################################################################################################################

# specify the category from the ind file that you want to color the dots according to  
ggplot(ind, aes(x=PC1, y=PC2, color=Oil)) +
  geom_point(shape=1)    # Use hollow circles


### Now see the relationship between traits of interest identifyied in the PCA and Branched/unbranched
### we will do a t-test and follow it up with a boxplot

t.test(traits$Circular ~ traits$Branching)

traits$Branching <- as.character(traits$Branching)

ggplot(traits, aes(x=Branching, y=Circular)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4)








