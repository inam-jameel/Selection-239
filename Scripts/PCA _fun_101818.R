### Set working directory
setwd(dir = "/volumes/pbio-burke/users/ijameel/inam_SELECTION_239/")

library(ggplot2)
library(PerformanceAnalytics)
library(dplyr)
library(ade4)

### read in individual data with information for heterotic groups/oil cont etc
ind <- read.csv(file = "239_PCA.csv", header = TRUE)

### Read in your trait data and only keep original measuments and their differences by environment
### Note: this removes the first row as well
### also, the Salty..trimmed_traits.csv file is already trimmed from the original (SALTYHEL17_traits)
### so refer to that if you need to start will all the data

traits <- read.csv("Analysis_SaltyHel17/Data/SaltyHel17_trimmed_traits.csv", header = TRUE)
salt_traits <- select(traits, contains("salt"))
control_traits <- select(traits, contains("water"))
diff_traits <- select(traits, contains("_diff"))
Salt_control_diff_traits <- bind_cols(salt_traits,control_traits,diff_traits)
Salt_control_diff_traits <- select(Salt_control_diff_traits, -contains("ln"))

# add the SAM line info from the data set back in and rename to "traits", to keep it simple
KeptSam <- read.csv("Analysis_SaltyHel17/Data/SaltyHel17_trimmed_traits.csv", header = TRUE)
KeptSam <- select(KeptSam, (2))
traits <- bind_cols(KeptSam,Salt_control_diff_traits)

# remove lines with any missing data
traits <- traits[complete.cases(traits),] 

# run the pca and set 4 for axes, remove the line info 
tpca <- dudi.pca(traits[,-c(1)])
4
# make data frames for the "li" (row coordinates) and "co" (principal components). 
# refer to the dudi.pca help for more info 
li <- tpca$li
co <- tpca$co

# remove lines in the 239_PCA data that are not in traits 
ind <- subset(ind,SAM_LINE %in% traits$Plant)

# add the li data frame to the individual data
ind <- bind_cols(ind,li)

tail (ind)


################################random ggplot##############

ggplot(ind, aes(x=Axis2, y=Axis4, color=Lineage)) +
  geom_point(shape=1)    # Use hollow circles

