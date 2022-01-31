#Association between Methylation and PSS: EWAS
#Analyst: Ritika Sharma
#Updated: 31.01.2022

#####################################
options(stringsAsFactors = FALSE)
rm(list=ls())

##Packages################
library(openxlsx)
library(dplyr)
library(tibble)
library(ggplot2)
library(parallel)
library(readxl)
library(affycoretools)
library(BiocStyle)
library(xtable)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
data(SNPs.147CommonSingle)
library(sva)
library(limma)
library(RColorBrewer)
library(DMRcate)
library(Homo.sapiens)
library(missMethyl)

##define paths###################

idatdir <-  "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/EPIC data/infos/ScanData/"
rdir <- "N:/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/Results_models_FELICITy/" ##results 
basedir <- "N:/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/EPIC data/"  #basedirectory

##loading the phenotype data######################
samps <-  read.csv(paste0(basedir,"Samps_FELICITy.csv"),fill=T , stringsAsFactors=F) 

##loading normalised methylation data
load(paste0(idatdir, "methdataall.Rdata"))
dat 

annoEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(annoEPIC)

#matching cpg site ronames of mvals and annoepic
annoEPICSub <- annoEPIC[match(rownames(getM(eset.nosex)),annoEPIC$Name),                    
                        c(1:4,12:19,22:ncol(annoEPIC))]                       


################DMPs:- PSS##########################################################################################################
##checking variables##
sum(is.na(samps$PSS))#0

##Creating a design model matrix
dim(model.matrix(~ PSS + Gender + Autoimmunediseases + Gestational_Diabetes  + Smoking +  Gestational_age_birth + Illumina.Plate, samps))

design <- model.matrix(~ PSS + Gender + Autoimmunediseases + Gestational_Diabetes + Smoking + Gestational_age_birth +Illumina.Plate, samps)

##Generating surrogate variables
sv <- sva(getM(eset.nosex)[,!is.na(samps$Cortisol)], design,vfilter = 5e4) #3

##combining the deisgn matrix and surrogate variables
design <- cbind(design, sv$sv)

##Finally fitting in the limma model
fit <- lmFit(getM(eset.nosex), design)
fit2 <- eBayes(fit)

colnames(fit2)

topTable(fit2, 2) 

DMPs <- topTable(fit2, coef=2,confint = TRUE, genelist=annoEPICSub, number= 808554) 

bonfc <- 0.05/nrow(fit2)
bonfc
#[1] 6.183879e-08
sum(DMPs$P.Value < bonfc) #0
sum(DMPs$adj.P.Val < 0.05) #0

##saving the results
write.table(DMPs,paste0(rdir, "/PSS/PSS_FELICITy_DMPs.csv"), sep=",", row.names=FALSE)

save(DMPs, samps,  file = paste0(rdir,"DMPS_PSS_adjustedmod.Rdata"))
#################################################################################################
