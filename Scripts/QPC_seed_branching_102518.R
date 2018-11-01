### This script is intended to apply the QPC analysis from Joesphs et al 2018 to trait/genetic data
### but also compare within branching groups. SO SPECIFY BRANCHING AT INDICATED (***) LOCATIONS!!!

### "unbranched_" indicates datasets with unbranched plants 
### "branched_"   indicates datasets with branched plants
### "All_"        indicates datasets with all plants

setwd(dir = "/volumes/pbio-burke/users/ijameel/inam_SELECTION_239/")

library(ggplot2)
library(dplyr)
#library(LaCroixColoR)
library(viridis)

### load Emily Josephs' scripts
source("quaint-master/R/calcQpc.R")
source("quaint-master/R/make_k.R")
source("quaint-master/R/var0.R")

### I skipped the data wrangling part since this was done in the earlier script
#############################################################################################################################
#############################################################################################################################
########################################### THE ANALYSIS ####################################################################
#############################################################################################################################
#############################################################################################################################

### Read in your genetic data
gidata <- read.csv("QPC/239_qpc.csv", header = TRUE)

### Read in your trait data
traits <- read.csv("Analysis_Seeds/Data/239_shape_mean_102318.csv", header = TRUE)
traits <- traits[complete.cases(traits),] #remove individuals with ANY missing data

# make the subset branching datasets. specify the select fuction to include any rows you like 
branched_t <- subset(traits, Branching > 0, select=Line:AverageChroma)
unbranched_t <- subset(traits, Branching == 0, select=Line:AverageChroma)
All_t <- subset(traits, select=Line:AverageChroma)

### Remove individuals in SNP dataset that were removed from the branched 
### and unbranched trait dataset 
branched_g <- subset(gidata, LINE %in% branched_t$Plant)
unbranched_g <- subset(gidata, LINE %in% unbranched_t$Plant)
All_g <- subset(gidata, LINE %in% All_t$Plant)

### Replace 2's with 1's
branched_g <- as.matrix(branched_g[,-c(1,2)])
branched_g[branched_g=="2"] <- 1 # we recoded to presence/absence of minor allele (i.e., 2's were converted to 1)

unbranched_g <- as.matrix(unbranched_g[,-c(1,2)])
unbranched_g[unbranched_g=="2"] <- 1

All_g <- as.matrix(All_g[,-c(1,2)])
All_g[All_g=="2"] <- 1

### estimate kinship matrix, specify branching!!! ***
myK <- make_k(as.matrix(
  #unbranched_g
  #branched_g
  All_g
  ))

### draw a heatmap
heatmap(myK, col=inferno(100))

### do the PCA
eigF <- eigen(myK)
myU <- eigF$vectors
myLambdas <- eigF$values

### calculate the PC cutoffs
varexp <- myLambdas/sum(myLambdas)
sumexp <- sapply(1:length(varexp), function(x){sum(varexp[1:x])})

### estimate PCs to retain at sum 30 or 50% variance explained. I usually start with 50% to check things out
pcmax <- which(sumexp > 0.3)[1]
pcmax

### estimate # PCs to retain using more sophisticated methods
library(nFactors) # for the scree-plot test
myscree <- nScree(myLambdas)
myscree$Components[1] 

library(AssocTests) #for the tracy-widom test
mytw <- tw(myLambdas, eigenL = length(myLambdas))
mytw$SigntEigenL

### plot the variance explained by each PC and identify the cutoffs. Adjust this to match earlier value
pcmax <- which(sumexp > 0.3)[1]

plot(myLambdas, bty="n", xlab = "PCs", ylab = "Eigenvalues")
abline(v = pcmax, col = viridis(6)[3], lwd=2)
abline(v = myscree$Components[1], col = viridis(6)[4], lwd=2)
abline(v = mytw$SigntEigenL, col = viridis(6)[5], lwd=2)

# specify the PC percentage thing here 
legend('topright',c('30% var explained','Castells rule','Tracy Widom'), col = viridis(6)[3:5], lwd=2, bty="n")

#remove the last end of PCs since these are likely to be extra variable
tailCutoff <- round(.9*length(myLambdas))

#############################################################################################################################
#############################################################################################################################
########################################### Analysis and heatmap ############################################################
#############################################################################################################################
#############################################################################################################################
library(qvalue)

### extract trait names and turn into a vector. Have to remove any columns in the begining that arent measured traits
traitdata <- read.csv("Analysis_Seeds/Data/239_shape_mean_102318.csv", header = FALSE, stringsAsFactors = FALSE)
traitlabels <- as.character(traitdata[1,-c(1:4)])

# specify branching ***
mydfs <- apply(
  #unbranched_t
  #branched_t
  All_t
  [,-c(1:3)], 2, function(x){calcQpc(
  myZ = x, myU = eigF$vectors, myLambdas = eigF$values, myPCcutoff = 0.3
)})

getqvalues <- function(ptable){
  qobj = qvalue(p = c(ptable))
  myqvals = matrix(qobj$qvalues, nrow=dim(ptable)[1])
  return(myqvals)
}

# specify the dfs.., 1 : n 
allpvals <- sapply(1:22, function(x){mydfs[[x]]$pvals})
myqvals <- getqvalues(allpvals)

# now generate the heatmap
layout(matrix(1, nrow=1, ncol=1))
mysig2 <- cut((1:1000/1000), c(0,0.001,0.01,0.05,0.1,1)) #for legend
par(mar=c(8,14,2,2), xpd=T, mfrow=c(1,1))
mycol <- c(viridis(6, direction=1)[1:4], "white")
image(allpvals, col=mycol, xaxt="n", yaxt="n", bty="n", breaks=c(0,0.001,0.01,0.05,0.1,1))
axis(1, at=seq(0,1, length = pcmax), las=2, label = 1:pcmax)
# make sure you update the y axis trait range, mydfs value (0:n-1)/n-1 
axis(2, at=(0:21)/21, labels = traitlabels, las=2)
# adjust the legend coordinates here, (x,y)
legend(-0.65,-0.15, levels(mysig2), fill=mycol, bty="n", horiz=T)


#############################################################################################################################
#############################################################################################################################
########################################### PLOT BY TRAIT *IN PROGRESS ######################################################
#############################################################################################################################
#############################################################################################################################
pnames <- as.character(traits$Plant)

calcCIs <- function(myName, myBlups=bluptable, myU=eigF$vectors, myLambdas=eigF$values){
  myZ = myBlups[,myName][1:239]
  myZ = myZ - mean(myZ)
  myBm = myZ %*% myU
  myCm = myBm/sqrt(myLambdas)
  myVa = var0(myCm[(tailCutoff-50):tailCutoff])
  myCI = sqrt(myVa*myLambdas)
  return(myCI)}

### combine trait data with eigenvectors 
mydf <- data.frame(traits[-nrow(traits),], eigF$vectors, stringsAsFactors = F)

### estimate confidence intervals for traits
variance <- var0(mydfs$Plant_Height_.cm.$cm) #change this hardcoded bit for your character of interest
myCI <- 1.96*sqrt(variance*eigF$values)

### chose colors
palette(c('#F6511D', "#F6511D", "#00A6ED", "#7FB800")) # choose colors

par(mar=c(5,5,1,1), mfrow=c(1,1))

plot(mydf$X10, mydf$Plant_Height_.cm., bty="n", xlab = "PC 10", ylab = "Plant Height (cm)", yaxt="n", col=as.factor(mydf$Lineage), lwd=2)
axis(2, las=2)
abline(lm(mydf$Plant_Height_.cm. ~ mydf$X20), col=viridis(6)[1], lwd=2)
abline(a=mean(mydf$Plant_Height_.cm.), b = 1.96*myCI[10], lty=2, col=viridis(6)[3], lwd=2)
abline(a=mean(mydf$Plant_Height_.cm.), b = -1.96*myCI[10], lty=2, col=viridis(6)[3], lwd=2)





#############################################################################################################################
#############################################################################################################################
########################################### DONT RUN ANY OF THIS CODE #######################################################
#############################################################################################################################
#############################################################################################################################
### doing the eigen decomposition
myEig <- eigen(myK)
plot(myEig$vectors[,1], myEig$vectors[,2], bty="n", xlab = "PC1", ylab = "PC2", col = lacroix_palette('Mango')[1])
plot(myEig$values/sum(myEig$values)*100, col = lacroix_palette('Mango')[3], bty="n", ylab = "% variation explained by each PC", xlab = "PC")

### Running the QPC analysis
myQpc <- calcQpc(myZ = traits$Plant_Height_cm,
                myU = myEig$vectors, 
                myLambdas = myEig$values,
                myR = 40,
                myPCcutoff = 0.5)

### Plot log10 p values by PC
plot(-log10(myQpc$pvals), bty="n", xlab = "PCs", ylab = "-log10(p value)", col = lacroix_palette('Mango')[4], lwd=2, xaxt="n")
abline(h = -log10(0.05/length(myQpc$pvals)), col = lacroix_palette('Mango')[1], lwd=2)
axis(1, at = c(1:length(myQpc$pvals)))

### Estimate the confidence intervals
myVaest <- var0(myQpc$cm)
myCI <- 1.96*sqrt(myVaest*myEig$values)

#plot
palette(c('#999999', '#E69F00', '#56B4E9', "#009E73"))
par(mar = c(5,8,5,14), xpd=T)

plot(myEig$vectors[,20], traits$Plant_Height_cm[-nrow(traits)], bty="n", col = traits$Lineage, lwd=2, ylab = "", yaxt="n",xlab = "PC20", cex.lab=2, cex.axis=2, xaxt="n")
axis(1, cex.axis=1.5, lwd=2)
axis(2, las=2, cex.axis=1.5, lwd=2)
mtext("Plant Height (cm)" ,side=2, line=5, cex=2)
legend(0.06, levels(traits$Lineage), pch=1, pt.lwd = 2,col = palette(), bty="n", text.width = 0.04)
par(xpd=F)
abline(lm(traits$Plant_Height_cm[-nrow(traits)]~myEig$vectors[,20]), lwd=2, col = "#0072B2")
abline(a=mean(traits$Plant_Height_cm), b = 1.96*myCI[20], lty=2, col='#56B4E9', lwd=2)
abline(a=mean(traits$Plant_Height_cm), b = -1.96*myCI[20], lty=2, col='#56B4E9', lwd=2)