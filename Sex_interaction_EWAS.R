#Gender interaction analysis: EWAS: PSS; PDQ; Cortisol and FSI
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

##############################################################
###PSS:- Interaction analysis with Gender########################################################################
dim(model.matrix(~ Gender*PSS   + Smoking + Autoimmunediseases + Gestational_Diabetes + Gestational_age_birth + Illumina.Plate, samps))
#114 , 7
design <- model.matrix(~ Gender*PSS + Smoking + Autoimmunediseases + Gestational_Diabetes + Gestational_age_birth +Illumina.Plate, samps)

sv <- sva(getM(eset.nosex), design, vfilter = 5e4)

design <- cbind(design, sv$sv)
fit <- lmFit(getM(eset.nosex), design)
fit2 <- eBayes(fit)

colnames(fit2)
#[1] "(Intercept)"           "Gender2"               "PSS"
#[4] "Smoking1"              "Autoimmunediseases1"   "Gestational_Diabetes1"
#[7] "Gestational_age_birth" "Illumina.Plate2"       "Gender2:PSS"
#[10] ""                      ""                      ""


topTable(fit2, 9, confint = TRUE)
#                 logFC        CI.L        CI.R    AveExpr         t
#cg09723184   0.03522617  0.02054088  0.04991146  3.9491537  4.755296
#cg27293447 -0.04252191 -0.06038872 -0.02465510 -0.6118012 -4.718026
#                P.Value adj.P.Val        B
#cg09723184 6.225521e-06 0.8742161 3.048696
#cg27293447 7.246262e-06 0.8742161 2.902317

  
DMPs <- topTable(fit2, coef=9,confint = TRUE, genelist=annoEPICSub, number= 808554) 
DMPs10 <- topTable(fit2, 9,confint = TRUE, genelist = annoEPICSub)

#savingtop 10 dmps##
write.table(DMPs10, paste0(rdir, "/PSS/PSSGenderinteractiontop10_FELICITy_DMPs.csv"), sep=",", row.names=FALSE)

###############################################
###PDQ:- Interaction analysis with Gender####################################################################

sum(is.na(samps$PDQ))#0

dim(model.matrix(~ Gender*PDQ   + Smoking + Gestational_age_birth + Autoimmunediseases + Gestational_Diabetes +  Illumina.Plate, samps))
#114 , 7
design <- model.matrix(~ Gender*PDQ + Smoking + Gestational_age_birth + Autoimmunediseases + Gestational_Diabetes + Illumina.Plate, samps)

sv <- sva(getM(eset.nosex), design, vfilter = 5e4)

design <- cbind(design, sv$sv)
fit <- lmFit(getM(eset.nosex), design)
fit2 <- eBayes(fit)

colnames(fit2)
#[[1] "(Intercept)"           "Gender2"               "PDQ"
#[4] "Smoking1"              "Gestational_age_birth" "Autoimmunediseases1"
#[7] "Gestational_Diabetes1" "Illumina.Plate2"       "Gender2:PDQ"
#[10] ""                      ""                      ""

topTable(fit2, 9)
#                logFC   AveExpr         t      P.Value  adj.P.Val        B
#cg03756940 -0.05546401 -3.056117 -5.695319 1.095793e-07 0.08327872 7.022675
#cg00008621 -0.04407063  2.720152 -5.554016 2.059942e-07 0.08327872 6.408309
#cg02216363  0.03598264  1.751909  5.161649 1.139019e-06 0.19114349 4.748194

bonfc <- 0.05/nrow(fit2)
bonfc
#[1] 6.183879e-08
sum(DMPs$P.Value < bonfc) #0
sum(DMPs$adj.P.Val < 0.05) #0
DMPs10 <- topTable(fit2, 9,confint = TRUE, genelist = annoEPICSub)

#savingtop 10 dmps##
write.table(DMPs10, paste0(rdir, "/PDQ/PDQGenderinteractiontop10_FELICITy_DMPs.csv"), sep=",", row.names=FALSE)


################FSI:- Interaction analysis with Gender#############################

dim(model.matrix(~ Gender*HRV_FSI + Smoking + Gestational_age_birth + Autoimmunediseases + Gestational_Diabetes + Illumina.Plate, samps))
#[1] 96 9
design <- model.matrix(~ Gender*HRV_FSI + Smoking + Gestational_age_birth + Autoimmunediseases + Gestational_Diabetes + Illumina.Plate, samps)

sv <- sva(getM(eset.nosex)[,!is.na(samps$HRV_FSI)], design,vfilter = 5e4)

design <- cbind(design, sv$sv)
fit <- lmFit(getM(eset.nosex)[,!is.na(samps$HRV_FSI)], design)
fit2 <- eBayes(fit)

colnames(fit2)
#[1] "(Intercept)"           "Gender2"               "HRV_FSI"
#[4] "Smoking1"              "Gestational_age_birth" "Autoimmunediseases1"
#[7] "Gestational_Diabetes1" "Illumina.Plate2"       "Gender2:HRV_FSI"
#[10] ""                      ""                      ""

dim(fit2)
#[1] 808554     12

topTable(fit2, 9)
#               logFC    AveExpr         t      P.Value adj.P.Val        B
#cg24715106  0.5442109 -4.2179441  5.586518 2.511550e-07 0.2030724 4.099577
#cg23782719 -0.4247221  3.6535857 -5.030064 2.539982e-06 0.5635324 2.596753
#cg04409895  0.4667505  0.4614848  4.972631 3.203975e-06 0.5635324 2.445075


DMPs <- topTable(fit2, coef=9,confint = TRUE, genelist=annoEPICSub, number= 808554) 

bonfc <- 0.05/nrow(fit2)
bonfc
#[1] 6.183879e-08
sum(DMPs$P.Value < bonfc) #0
sum(DMPs$adj.P.Val < 0.05) #0
DMPs10 <- topTable(fit2, 9,confint = TRUE, genelist = annoEPICSub)

write.table(DMPs10,paste0(rdir, "/FSI/FSItop10interactiongender_FELICITy_DMPs.csv"), sep=",", row.names=FALSE)

################Cortisol:- Interaction analysis with Gender#############################

dim(model.matrix(~ Gender*Cortisol + Smoking + Gestational_age_birth + Autoimmunediseases + Gestational_Diabetes + Illumina.Plate, samps))
#[1] 107 7
design <- model.matrix(~ Gender*Cortisol + Smoking + Gestational_age_birth + Autoimmunediseases + Gestational_Diabetes + Illumina.Plate, samps)

sv <- sva(getM(eset.nosex)[,!is.na(samps$Cortisol)], design,vfilter = 5e4)

design <- cbind(design, sv$sv)
fit <- lmFit(getM(eset.nosex)[,!is.na(samps$Cortisol)], design)
fit2 <- eBayes(fit)

colnames(fit2)
#[1] "(Intercept)"           "Gender2"               "Cortisol"
#[4] "Smoking1"              "Gestational_age_birth" "Autoimmunediseases1"
#[7] "Gestational_Diabetes1" "Illumina.Plate2"       "Gender2:Cortisol"
#[10] ""                      ""                      ""

dim(fit2)
#[1] 808554     12

topTable(fit2, 9)
#                     logFC   AveExpr         t      P.Value adj.P.Val         B
#cg18197866     0.003298714  6.372074  5.532416 2.547095e-07 0.1685130 3.7826033
#cg20460797    -0.006006267 -7.026308 -5.418944 4.168257e-07 0.1685130 3.3000527
#ch.3.2741911R -0.004988868 -4.673911 -5.146369 1.332801e-06 0.3592138 2.1636729


DMPs <- topTable(fit2, coef=9,confint = TRUE, genelist=annoEPICSub, number= 808554) 

bonfc <- 0.05/nrow(fit2)
bonfc
#[1] 6.183879e-08
sum(DMPs$P.Value < bonfc) #0
sum(DMPs$adj.P.Val < 0.05) #0
DMPs10 <- topTable(fit2, 9,confint = TRUE, genelist = annoEPICSub)

write.table(DMPs10,paste0(rdir, "/Cortisol/Cortisoltop10interactiongendr_FELICITy_DMPs.csv"), sep=",", row.names=FALSE)

###########################################################################################################################


