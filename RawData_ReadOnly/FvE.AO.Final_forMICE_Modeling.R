# === Description ==================================================================
# This script (and accompanying data) should be sufficient to recreate the analyses in Siddiqui et al. (under review, Proceedings B)
# Built 13 June 2023

# === 0) Preparation, Housekeeping, etc. ======
## empty environment
rm(list = ls())

## Set working directory. Hint, open this R script in the same folder as the data.
working.dir <- getwd()
setwd(working.dir)

## libraries
library(lme4)


##Global Variables
#variables
clrs <- c("#000000", "#FF0000") # c("black", "red")
options(max.print=1000000)
#Lake abbreviations in the study
LS_lake_IDs <- c("Com", "Fre", "Joe", "Nor", "Swa", "The", "Vib")
SP_lake_IDs <- c("Eno", "LQx", "Pax", "Pri", "Emi")
SP_lake_names <- c("Enos", "Little Quarry", "Paxton", "Priest", "Emily")
time.since.colonization <- read.csv("a_raw.data.readonly/KSampleMeanTimes.csv", header = TRUE, stringsAsFactors = FALSE)[, c(1,3)]
#function for later
f.sample.size.noNAs <- function(x) length(x) - sum(is.na(x))
#end 0)
# **********

# === 1) Data clean up and and Size Correction ======
## Read in data
#Using the scripts 2a.raw.data.manip.FvE.R and 2b.outliers.size.correction.FvE.R, and the data called in those scripts, we previously checked data for outliers greater than 3sd from the mean. Outliers were checked. If they were part of a set of outliers contained in a single indvidual, they were left alone. If they were obvious data collection or entry errors, they were fixed. Any outliers that remained were set to NA.
data.for.sc <- read.csv("b_cleaned.data/AO.sexualdimorphismPruned.for.sizecorrection.csv", header = TRUE, stringsAsFactors = FALSE)

#checking lake names
unique(substr(data.for.sc$univID, 1, 3))

# running into the problem of having L stand for both Limnetic and Lake. Change L for Lake to a G, to represent the unimodal, Generalist populations found in the non-species pair lakes.
substr(data.for.sc$univID, 5,5)[substr(data.for.sc$univID, 1, 3) %in% LS_lake_IDs] <- "G"
data.for.sc[,c(1:3)]

#remove lake-stream lake fish whose sex is uncertain (i.e., represented by lowercase "m" and "f", keeping only fish whose sex was certain from dissection)
data.for.sc <- data.for.sc[-which(substr(data.for.sc$univID, 7, 7) %in% c("m", "f")), ]

# === * 1a) Check for NA vs. 0 errors =======
#Double checking that that NAs for continuous traits that are biologically absent (rather than broken, etc) are actually zeroes. That is, tpg for a fish with ips = 0 should be 0 ditto pelvic spine lengths, etc. This is because evolved loss could be treated as zeros, the result of a developmental process that aborted before a trait became measurble. Thus, an observed absence is a data point of 0, which is different thing than missing data (NA). This problem is pervasive in the fossil data, where MAB did not make this distinction. Thus, focus only on the fossil data here. Many of these corrections were done during early cleaning steps but one last check here.

#For fossils, when pelvic score is less than 3 (ips < 3), pelvic spines are lost and should have lengths of 0.
#pelvic spine lengths
data.for.sc[which(substr(data.for.sc$univID, 1, 3) == "Fos" & data.for.sc$ips < 3) , c(1, 3, 7)] #for ips<3, set lps to 0. for ips = 0.0, set tpg = 0.0.
sum(data.for.sc$lps[substr(data.for.sc$univID, 1, 3) == "Fos" & data.for.sc$ips < 3], na.rm = T) #0. check that ips<3 corresponds to lps of 0.
sum(data.for.sc$tpg[substr(data.for.sc$univID, 1, 3) == "Fos" & data.for.sc$ips == 0], na.rm = T) #0. check. ips = 0 is complete pelvic loss, so tpg should be 0.

#mds and dorsal spine lengths
##when mds == 0, all three dorsal spine lengths should be set to 0.
data.for.sc[data.for.sc$mds == 0 , c(1, 4, 16, 17, 19) ] #when mds = 0, and ds1, ds2, and ds3 all equal NA. Change those to 0. No extant fish have mds = 0, so this is limited to fossils.
#ds1
ds1NA.mds0 <- which(data.for.sc$mds == 0)[which(data.for.sc$mds == 0) %in% which(is.na(data.for.sc$ds1)) ]
data.for.sc$ds1[ds1NA.mds0 ] <- 0
#ds2
ds2NA.mds0 <- which(data.for.sc$mds == 0)[which(data.for.sc$mds == 0) %in% which(is.na(data.for.sc$ds2)) ]
data.for.sc$ds2[ds2NA.mds0 ] <- 0
#ds3
ds3NA.mds0 <- which(data.for.sc$mds == 0)[which(data.for.sc$mds == 0) %in% which(is.na(data.for.sc$ds3)) ]
data.for.sc$ds3[ds3NA.mds0 ] <- 0

# when mds == 0 but at least one spine length is greater than 0, change mds to the number of values > 0
data.for.sc[data.for.sc$mds == 0 , c(1, 4, 16, 17, 19) ] #ds3 for K.03.147 > 0 but mds = 0
which(data.for.sc$mds == 0 & (data.for.sc$ds1>0 | data.for.sc$ds2>0 | data.for.sc$ds3>0)) #Fos.K.03.147
data.for.sc$mds[data.for.sc$univID == "Fos.K.03.147"] <- 1

#when mds = 2, do we ever have three ds lengths/values? 
data.for.sc[which(data.for.sc$mds ==2 & data.for.sc$ds1 > 0 & data.for.sc$ds2 > 0 & data.for.sc$ds3 > 0), c(1, 4, 16, 17, 19)]
#Yes. 5 extant fish.
data.for.sc$mds[which(data.for.sc$mds ==2 & data.for.sc$ds1 > 0 & data.for.sc$ds2 > 0 & data.for.sc$ds3 > 0)] <- rep(3,5)

#when mds = 1, do we ever have more than 1 ds length/value?
data.for.sc[data.for.sc$mds == 1 , c(1, 4, 16, 17, 19) ]
data.for.sc$mds[data.for.sc$univID == "Fos.K.08.176"] <- 2
data.for.sc$mds[data.for.sc$univID == "Fos.K.10.091"] <- 2 ; data.for.sc$ds1[data.for.sc$univID == "Fos.K.10.091"] <- 0

#when mds = 1 and neither ds2 or ds1 has a length, check that ds1=ds2=0.
data.for.sc[which(data.for.sc$mds == 1 & data.for.sc$ds3 >= 0) , c(1, 4, 16, 17, 19) ]

#fixing other errors when mds = 1
data.for.sc[which(data.for.sc$mds == 1), c(1, 4, 16, 17, 19)]
data.for.sc$ds1[data.for.sc$univID=="Fos.K.10.130"] = 0 ; data.for.sc$ds3[data.for.sc$univID=="Fos.K.10.130"] = 0
data.for.sc$ds1[data.for.sc$univID=="Fos.K.11.050"] = 0 ; data.for.sc$ds3[data.for.sc$univID=="Fos.K.11.050"] = 0
data.for.sc$ds1[data.for.sc$univID=="Fos.K.13.318"] = 0 ; data.for.sc$ds3[data.for.sc$univID=="Fos.K.13.318"] = 0
data.for.sc$ds2[data.for.sc$univID=="Fos.K.14.147"] = 0 ; data.for.sc$ds3[data.for.sc$univID=="Fos.K.14.147"] = 0

#when mds = 2, fixing errors
data.for.sc[which(data.for.sc$mds == 2), c(1,4,15, 16, 17)]
data.for.sc$ds2[data.for.sc$univID=="Fos.K.T0.063"] = 0
data.for.sc$ds2[data.for.sc$univID=="Fos.K.16.345"] = 0
data.for.sc$ds2[data.for.sc$univID=="Fos.K.16.372"] = 0
data.for.sc$ds2[data.for.sc$univID=="Fos.K.03.077"] = 0
data.for.sc$ds1[data.for.sc$univID=="Eno.L.X.118"] = 0

#when mds = 3, fixing errors
data.for.sc[which(data.for.sc$mds == 3), c(1,4,15, 16, 17)]
#none

#Other dorsal spine errors
data.for.sc$mds[data.for.sc$univID=="Fos.K.10.145"] = 1 ; data.for.sc$ds1[data.for.sc$univID=="Fos.K.10.145"] = 0 ; data.for.sc$ds2[data.for.sc$univID=="Fos.K.10.145"] = 0
data.for.sc$mds[data.for.sc$univID=="Fos.K.10.112"] = 1 ; data.for.sc$ds1[data.for.sc$univID=="Fos.K.10.112"] = 0 ; data.for.sc$ds2[data.for.sc$univID=="Fos.K.10.112"] = 0

str(data.for.sc)#check
unique(data.for.sc$mds) #check

head(data.for.sc)
tail(data.for.sc)
unique(substr(data.for.sc$univID, 1, 8))
names(data.for.sc)

#remove ips and ips.prop since they aren't comparable between fossil and extant populations. That is, they were collected in slightly different ways, so that a value of 1 for fossils is not the same as a value of 1 for extant fish, even when normalized to proportions.

data.for.sc <- data.for.sc[ , -which(names(data.for.sc) %in% c("ips", "ips.prop"))]
#end 1a)
# **********

# === * 1b) Size Correction ======================================
names(data.for.sc)
# size correct the following continuous traits
continuous.traits.to.size.correct <- c("lps", "ect", "tpg", "cle", "pmx", "ds1", "ds2", "lpt", "ds3")
# remove biological zeros from the following traits so as not to disrupt the size correction model. Then put back the 0s
continuous.traits.with.zero.class <- c("tpg", "lps", "ds1", "ds2", "lpt", "ds3")

#function
# === * f.sizecorrect ==========
# from: OKE ET AL COMMON GARDEN MANUSCRIPT METHOD (with small changes made)
# "All relative warps and univariate shape traits were allometrically 
# standardized to a common body size (Reist 1985; Lleonart et al. 2000) based 
# on Ms = M0 (Ls / L0)^b, where Ms is the standardized trait value (mm), M0 is 
# the non-standardized trait value (mm), Ls is the overall mean centroid size 
# (for RWs) or standard length (mm, for univariate shape traits), and L0 is the 
# centroid size or standard length (mm) of the individual. The common 
# within-group slope, b, was calculated from a linear mixed model of log10(M0)
# regressed on log10(L0), with group included as a random factor (Reist 1985;
# Lleonart et al. 2000)."
# - Oke et al.  The size correction formula described above: Ms = M0(Ls/L0)^b

# variables:
# standard.length = vector with lengths to be used for correcting other vars 
# random.effect.substr.char = the start and end of the substring that contains the grouping variables (species pair vs. fossils)
# dataframe.to.use = dataframe name with all variables
# vector.of.columns = vector with position of variables to be corrected. do not include standard length.
# zero.class.traits = the trait names for traits that have lots of evolutionarily lost data (zeros) that will affect the regression.

# The output will have all the original columns from dataframe.to.use with the size-corrected columns added on at the end.
# requires package lme4 

#example data frame: data.for.sc
#col.with.standard.length = 2
#random.effect.substr.char = c(1, 3)
#dataframe.to.use = data.for.sc
#vector.of.columns = which(names(data.for.sc) %in% continuous.traits.to.size.correct)
#zero.class.traits = continuous.traits.with.zero.class

f.sizecorrect <- function(col.with.standard.length, vector.of.columns, dataframe.to.use, random.effect.substr.char, zero.class.traits) {
  # Calculate common, within-group slope, b, for each continuous trait, defined by vector of columns
  b.vector.lmm <- vector()
  for(i in 1:length(vector.of.columns)){
    abcd <- dataframe.to.use[ , c(1, col.with.standard.length, vector.of.columns[i])]
    #plot data with 0s for check
    plot(abcd[,3] ~ abcd[,2], main = paste(colnames(abcd)[3], " with 0s"))
    #if a trait with zero-class that is biologically meaningful, remove zeros.
    if(colnames(abcd)[3] %in% zero.class.traits & length(which(abcd[ , 3] == 0)) > 0) {
      abcd <- abcd[-which(abcd[ , 3] == 0) , ]
    }
    #plot data after 0s removed for check.
    plot(abcd[, 3] ~ abcd[, 2], main = paste(colnames(abcd)[3], " without 0s"))
    
    # Here, I am treating each Fos, species pair, generalist lake as groups for random factor
    random.effect <- substr(x = abcd[ , 1], start = random.effect.substr.char[1], stop = random.effect.substr.char[2])
    log10.trait <- log10(abcd[ , 3])
    b.model <- lmer(log10.trait ~ log10(abcd[ , 2]) + (1|random.effect))
    plot(log10.trait ~ log10(abcd[ , 2]), main = colnames(abcd)[3]) ; abline(b = coef(summary(b.model))[2,1], a = coef(summary(b.model))[1,1]) #check that line goes through data.
    b <- coef(summary(b.model))[2,1]
    
    #append b for trait i
    b.vector.lmm <- c(b.vector.lmm, b)
  }
  
  # size correct
  xx <- dataframe.to.use  
  columnnames <- colnames(xx)
  j <- 1
  for (i in 1:length(vector.of.columns)) {
    #Calculate overall mean standard length (Ls), not including individuals with 0s from zero.class.traits
    if(colnames(xx)[vector.of.columns[i]] %in% zero.class.traits) Ls <- mean(xx[which(xx[vector.of.columns[i]] != 0) , 2], na.rm = TRUE) else Ls = mean(xx[ , 2], na.rm = TRUE)
    # Call individual standard length vector
    L0 <- xx[ , 2]
    #grab the appropriate column of data
    M0 <- xx[,vector.of.columns[i]] 
    #size correction formula
    Ms = M0 * ((Ls/L0)^b.vector.lmm[j]) 
    #rename, append
    j <- j + 1
    columnnames <- c(columnnames, paste(colnames(xx[vector.of.columns[i]]), "sc", sep = "."))
    xx <- cbind(xx, Ms)
  }
  colnames(xx) <- columnnames # Rename the columns in the temporary dataframe xx
  return(xx) # Output a new dataframe with the name provided in "outputfilename"
}
# end f.size.correct

## Size correction
data.sc <- f.sizecorrect(col.with.standard.length = 2, vector.of.columns = which(names(data.for.sc) %in% continuous.traits.to.size.correct), dataframe.to.use = data.for.sc, random.effect.substr.char = c(1,3), zero.class.traits = continuous.traits.with.zero.class)

#view
AO.analysis <- data.sc
head(AO.analysis); str(AO.analysis) ; unique(AO.analysis$mds)
unique(substr(AO.analysis$univID, 5, 5)) #benthic (B), limnetic (L),  fossil (K), generalist (G)
unique(substr(AO.analysis$univID, 7, 7)) #M, F, and some T,1,0 from fossil fish that don't have a character assigned to denote sex.
unique(substr(AO.analysis$univID, 1, 3)) #"Joe" "Nor" "Swa" "The" "Vib" "Emi" "LQx" "Pax" "Pri" "Fos" #Recall: no sex information for Enos Lake because those data weren't provided/collected by Schluter. Thus, only four benthic-limnetic lakes with sex information.

#end 1b)
# **********

#=== * 1c) Write data set for analysis =======
#20 June 2023
#now we have a combined file, Fossil + benthic-limnetic species pair plus + lake fish from lake-stream pairs. For extant fish, sex is noted in the 7th character of the ID. Ecotype (Benthic/limnetic/generalist) is noted in the 5th character. Lake ID is noted in the first three characters.
#write.csv(x = AO.analysis, file = "b_cleaned.data/AO.sexualdimorphismPruned.sc.forAnalysis_forMICE_Modeling_20June2023.csv", row.names = FALSE)
