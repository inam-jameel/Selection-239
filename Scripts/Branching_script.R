### Looking at Flower number from DREC experiment and infer branching from it

setwd(dir = "/volumes/pbio-burke/users/ijameel/inam_SELECTION_239/")

library(ggplot2)
library(dplyr)

### read in flowering number data, average branching above 1.5 is considered branching
f_num <- read.csv(file = "/volumes/pbio-burke/users/ijameel/inam_SELECTION_239/DREC_flower_data/Flower_num_dry.csv", stringsAsFactors = FALSE)
f_num <- as.data.frame(f_num)

### now we need to filter the flowering number data with the lines used in the 239 lines exps
# read in the line discription data
description <- read.csv(file = "239_PCA.csv")
description <- as.data.frame(description)

# then filter out the lines in the flowering number data that arent in the 239 line exps
filtered_f_num <- f_num %>% filter(Line %in% discription$Line)

# now attach the branching information and PPN info for sanity check to the description file
branching <- select(filtered_f_num, "Branching","Line")
f_num_description <- bind_cols(discription,branching)
write.csv(f_num_description, file = "239_PCA_102318.csv")

### after i made the new description "239_PCA" file, I manually went in excel and checked the data
### I also moved the branching column to the begining of the file, deleted the imported PPN info,
### replaced "0" and "1" with unbranching and branching respectively

### Now we need to add this branching column to other data sets we are interested in
### we ultimately want to analyze traits within  lines that are only branched or only unbranched  
### NOTE: the values in the branching column is going to be 0 and 1 for unbranched
### and branched respectively. 

### First will be the seed dataset. The "og" is short for original
seeds_og <- read.csv(file = "Analysis_Seeds/shape_mean.csv", stringsAsFactors = FALSE)

# this file has not filtered out the lines in the 239 exp, we will first add the unfiltered 
# flower data, then we will filter out the extra lines
branching_all_lines <- select(f_num, "Branching","Line")
seeds_branching <- bind_cols(seeds_og,branching_all_lines)
write.csv(seeds_branching, file = "Analysis_Seeds/shape_mean_102318.csv")

# filter the lines.
filtered_seeds_branching<- seeds_branching %>% filter(Line %in% discription$Line)
write.csv(filtered_seeds_branching, file = "Analysis_Seeds/239_shape_mean_102318.csv")

# now you should be able to use the new trait data and descriptoin file in your analysis! 