# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### set-up working environment ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#clear environment
rm(list = ls())

#load libraries
library(tidyverse)
library(ggdendro)
library(gdata)
library(UpSetR)
library(VennDiagram)

# load multiplot function
source("multiplot.R")

#set color palette
myPalette <- c("#000000","#999999", "#009E73", "#E69F00", "#0072B2", "#D55E00","#56B4E9","#CC79A7")
# grey, black, green, dark green, orange, dark orange, red, blue, light blue, pink

#load datasets
PAO1_vitro_ess <- read.csv('Dataset_S1.csv')
PA14_vitro_ess <- read.csv('Dataset_S2.csv')
PAO1_model_ess <- read.csv('Dataset_S3.csv')
PA14_model_ess <- read.csv('Dataset_S4.csv')
PAO1_updated_model_ess <- read.csv('Dataset_S10.csv')
PA14_updated_model_ess <- read.csv('Dataset_S11.csv')

#load model genes
PAO1_model <- read.csv('modelGenes_PAO1.csv')
PA14_model <- read.csv('modelGenes_PA14.csv')

colnames(PAO1_model) <- c('Genes')
colnames(PA14_model) <- c('Genes')

# #create a new dataframe that combines in vitro and in silico essentiality information for PAO1
PAO1_vUnique <- data.frame(Genes = setdiff(PAO1_vitro_ess$Genes, PAO1_model_ess$Genes), LB_Model = NA, Sputum_Model = NA, Pyruvate_Model = NA, Succinate_Model = NA)
PAO1_mUnique <- data.frame(Genes = setdiff(PAO1_model_ess$Genes, PAO1_vitro_ess$Genes), PAO1.LB.913 = NA, PAO1.LB.201 = NA, PAO1.LB.335 = NA, PAO1.Sputum.224 = NA, PAO1.Sputum.405 = NA, PAO1.Pyruvate.179 = NA, PAO1.Succinate.640 = NA)

PAO1_vMerge <- rbind(PAO1_vitro_ess,PAO1_mUnique)
PAO1_mMerge <- rbind(PAO1_model_ess, PAO1_vUnique)

PAO1_ess_complete <- merge(PAO1_vMerge, PAO1_mMerge)

PAO1_ess_complete[is.na(PAO1_ess_complete)] <- 0

# #create a new dataframe that combines in vitro and in silico essentiality information for PA14
PA14_vUnique <- data.frame(Genes = setdiff(PA14_vitro_ess$Genes, PA14_model_ess$Genes), LB_Model = NA, Sputum_Model = NA)
PA14_mUnique <- data.frame(Genes = setdiff(PA14_model_ess$Genes, PA14_vitro_ess$Genes), PA14.LB.1544 = NA, PA14.LB.634 = NA, PA14.Sputum.510 = NA)

PA14_vMerge <- rbind(PA14_vitro_ess,PA14_mUnique)
PA14_mMerge <- rbind(PA14_model_ess, PA14_vUnique)

PA14_ess_complete <- merge(PA14_vMerge, PA14_mMerge)

PA14_ess_complete[is.na(PA14_ess_complete)] <- 0

# create a new dataframe that filters the combined dataframe based on presence in the model for PAO1
PAO1_ess_model <- subset(PAO1_ess_complete, Genes %in% PAO1_model$Genes)
PAO1_mUnique <- data.frame(Genes = setdiff(PAO1_model$Genes, PAO1_ess_model$Genes), PAO1.LB.913 = NA, PAO1.LB.201 = NA, PAO1.LB.335 = NA, PAO1.Sputum.224 = NA, PAO1.Sputum.405 = NA, PAO1.Pyruvate.179 = NA, PAO1.Succinate.640 = NA, LB_Model = NA, Sputum_Model = NA, Pyruvate_Model = NA, Succinate_Model = NA)

PAO1_ess_model <- rbind(PAO1_ess_model, PAO1_mUnique)
PAO1_ess_model[is.na(PAO1_ess_model)] <- 0

# create a new dataframe that filters the combined dataframe based on presence in the model for PA14
PA14_ess_model <- subset(PA14_ess_complete, Genes %in% PA14_model$Genes)
PA14_mUnique <- data.frame(Genes = setdiff(PA14_model$Genes, PA14_ess_model$Genes), PA14.LB.1544 = NA, PA14.LB.634 = NA, PA14.Sputum.510 = NA, LB_Model = NA, Sputum_Model = NA)

PA14_ess_model <- rbind(PA14_ess_model, PA14_mUnique)
PA14_ess_model[is.na(PA14_ess_model)] <- 0

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Calculate statistics ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate statistics for PAO1
TP_PAO1.LB.913 <- which(ifelse(PAO1_ess_model$PAO1.LB.913==PAO1_ess_model$LB_Model & PAO1_ess_model$PAO1.LB.913==1, 1,0)!=0)
TN_PAO1.LB.913 <- which(ifelse(PAO1_ess_model$PAO1.LB.913==PAO1_ess_model$LB_Model & PAO1_ess_model$PAO1.LB.913==0, 1,0)!=0)
FP_PAO1.LB.913 <- which(ifelse(PAO1_ess_model$PAO1.LB.913!=PAO1_ess_model$LB_Model & PAO1_ess_model$PAO1.LB.913==0, 1,0)!=0)
FN_PAO1.LB.913 <- which(ifelse(PAO1_ess_model$PAO1.LB.913!=PAO1_ess_model$LB_Model & PAO1_ess_model$PAO1.LB.913==1, 1,0)!=0)

TP_PAO1.LB.201 <- which(ifelse(PAO1_ess_model$PAO1.LB.201==PAO1_ess_model$LB_Model & PAO1_ess_model$PAO1.LB.201==1, 1,0)!=0)
TN_PAO1.LB.201 <- which(ifelse(PAO1_ess_model$PAO1.LB.201==PAO1_ess_model$LB_Model & PAO1_ess_model$PAO1.LB.201==0, 1,0)!=0)
FP_PAO1.LB.201 <- which(ifelse(PAO1_ess_model$PAO1.LB.201!=PAO1_ess_model$LB_Model & PAO1_ess_model$PAO1.LB.201==0, 1,0)!=0)
FN_PAO1.LB.201 <- which(ifelse(PAO1_ess_model$PAO1.LB.201!=PAO1_ess_model$LB_Model & PAO1_ess_model$PAO1.LB.201==1, 1,0)!=0)

TP_PAO1.LB.335 <- which(ifelse(PAO1_ess_model$PAO1.LB.335==PAO1_ess_model$LB_Model & PAO1_ess_model$PAO1.LB.335==1, 1,0)!=0)
TN_PAO1.LB.335 <- which(ifelse(PAO1_ess_model$PAO1.LB.335==PAO1_ess_model$LB_Model & PAO1_ess_model$PAO1.LB.335==0, 1,0)!=0)
FP_PAO1.LB.335 <- which(ifelse(PAO1_ess_model$PAO1.LB.335!=PAO1_ess_model$LB_Model & PAO1_ess_model$PAO1.LB.335==0, 1,0)!=0)
FN_PAO1.LB.335 <- which(ifelse(PAO1_ess_model$PAO1.LB.335!=PAO1_ess_model$LB_Model & PAO1_ess_model$PAO1.LB.335==1, 1,0)!=0)

TP_PAO1.Sputum.224 <- which(ifelse(PAO1_ess_model$PAO1.Sputum.224==PAO1_ess_model$Sputum_Model & PAO1_ess_model$PAO1.Sputum.224==1, 1,0)!=0)
TN_PAO1.Sputum.224 <- which(ifelse(PAO1_ess_model$PAO1.Sputum.224==PAO1_ess_model$Sputum_Model & PAO1_ess_model$PAO1.Sputum.224==0, 1,0)!=0)
FP_PAO1.Sputum.224 <- which(ifelse(PAO1_ess_model$PAO1.Sputum.224!=PAO1_ess_model$Sputum_Model & PAO1_ess_model$PAO1.Sputum.224==0, 1,0)!=0)
FN_PAO1.Sputum.224 <- which(ifelse(PAO1_ess_model$PAO1.Sputum.224!=PAO1_ess_model$Sputum_Model & PAO1_ess_model$PAO1.Sputum.224==1, 1,0)!=0)

TP_PAO1.Sputum.405 <- which(ifelse(PAO1_ess_model$PAO1.Sputum.405==PAO1_ess_model$Sputum_Model & PAO1_ess_model$PAO1.Sputum.405==1, 1,0)!=0)
TN_PAO1.Sputum.405 <- which(ifelse(PAO1_ess_model$PAO1.Sputum.405==PAO1_ess_model$Sputum_Model & PAO1_ess_model$PAO1.Sputum.405==0, 1,0)!=0)
FP_PAO1.Sputum.405 <- which(ifelse(PAO1_ess_model$PAO1.Sputum.405!=PAO1_ess_model$Sputum_Model & PAO1_ess_model$PAO1.Sputum.405==0, 1,0)!=0)
FN_PAO1.Sputum.405 <- which(ifelse(PAO1_ess_model$PAO1.Sputum.405!=PAO1_ess_model$Sputum_Model & PAO1_ess_model$PAO1.Sputum.405==1, 1,0)!=0)

TP_PAO1.Pyruvate.179 <- which(ifelse(PAO1_ess_model$PAO1.Pyruvate.179==PAO1_ess_model$Sputum_Model & PAO1_ess_model$PAO1.Pyruvate.179==1, 1,0)!=0)
TN_PAO1.Pyruvate.179 <- which(ifelse(PAO1_ess_model$PAO1.Pyruvate.179==PAO1_ess_model$Sputum_Model & PAO1_ess_model$PAO1.Pyruvate.179==0, 1,0)!=0)
FP_PAO1.Pyruvate.179 <- which(ifelse(PAO1_ess_model$PAO1.Pyruvate.179!=PAO1_ess_model$Sputum_Model & PAO1_ess_model$PAO1.Pyruvate.179==0, 1,0)!=0)
FN_PAO1.Pyruvate.179 <- which(ifelse(PAO1_ess_model$PAO1.Pyruvate.179!=PAO1_ess_model$Sputum_Model & PAO1_ess_model$PAO1.Pyruvate.179==1, 1,0)!=0)

TP_PAO1.Succinate.640 <- which(ifelse(PAO1_ess_model$PAO1.Succinate.640==PAO1_ess_model$Sputum_Model & PAO1_ess_model$PAO1.Succinate.640==1, 1,0)!=0)
TN_PAO1.Succinate.640 <- which(ifelse(PAO1_ess_model$PAO1.Succinate.640==PAO1_ess_model$Sputum_Model & PAO1_ess_model$PAO1.Succinate.640==0, 1,0)!=0)
FP_PAO1.Succinate.640 <- which(ifelse(PAO1_ess_model$PAO1.Succinate.640!=PAO1_ess_model$Sputum_Model & PAO1_ess_model$PAO1.Succinate.640==0, 1,0)!=0)
FN_PAO1.Succinate.640 <- which(ifelse(PAO1_ess_model$PAO1.Succinate.640!=PAO1_ess_model$Sputum_Model & PAO1_ess_model$PAO1.Succinate.640==1, 1,0)!=0)

PAO1_ess_stats <- data.frame(stat = c(length(TP_PAO1.LB.913),
                                      length(TN_PAO1.LB.913),
                                      length(FP_PAO1.LB.913),
                                      length(FN_PAO1.LB.913),
                                      length(TP_PAO1.LB.201),
                                      length(TN_PAO1.LB.201),
                                      length(FP_PAO1.LB.201),
                                      length(FN_PAO1.LB.201),
                                      length(TP_PAO1.LB.335),
                                      length(TN_PAO1.LB.335),
                                      length(FP_PAO1.LB.335),
                                      length(FN_PAO1.LB.335),
                                      length(TP_PAO1.Sputum.224),
                                      length(TN_PAO1.Sputum.224),
                                      length(FP_PAO1.Sputum.224),
                                      length(FN_PAO1.Sputum.224),
                                      length(TP_PAO1.Sputum.405),
                                      length(TN_PAO1.Sputum.405),
                                      length(FP_PAO1.Sputum.405),
                                      length(FN_PAO1.Sputum.405),
                                      length(TP_PAO1.Pyruvate.179),
                                      length(TN_PAO1.Pyruvate.179),
                                      length(FP_PAO1.Pyruvate.179),
                                      length(FN_PAO1.Pyruvate.179),
                                      length(TP_PAO1.Succinate.640),
                                      length(TN_PAO1.Succinate.640),
                                      length(FP_PAO1.Succinate.640),
                                      length(FN_PAO1.Succinate.640)),
                             group = c('PAO1.LB.913','PAO1.LB.913','PAO1.LB.913','PAO1.LB.913',
                                       'PAO1.LB.201','PAO1.LB.201','PAO1.LB.201','PAO1.LB.201',
                                       'PAO1.LB.335','PAO1.LB.335','PAO1.LB.335','PAO1.LB.335',
                                       'PAO1.Sputum.224','PAO1.Sputum.224','PAO1.Sputum.224','PAO1.Sputum.224',
                                       'PAO1.Sputum.405','PAO1.Sputum.405','PAO1.Sputum.405','PAO1.Sputum.405',
                                       'PAO1.Pyruvate.179','PAO1.Pyruvate.179','PAO1.Pyruvate.179','PAO1.Pyruvate.179',
                                       'PAO1.Succinate.640','PAO1.Succinate.640','PAO1.Succinate.640','PAO1.Succinate.640'),
                             type = c('TP','TN','FP','FN',
                                      'TP','TN','FP','FN',
                                      'TP','TN','FP','FN',
                                      'TP','TN','FP','FN',
                                      'TP','TN','FP','FN',
                                      'TP','TN','FP','FN',
                                      'TP','TN','FP','FN'))

# calculate statistics for PA14
TP_PA14.LB.1544 <- which(ifelse(PA14_ess_model$PA14.LB.1544==PA14_ess_model$LB_Model & PA14_ess_model$PA14.LB.1544==1, 1,0)!=0)
TN_PA14.LB.1544 <- which(ifelse(PA14_ess_model$PA14.LB.1544==PA14_ess_model$LB_Model & PA14_ess_model$PA14.LB.1544==0, 1,0)!=0)
FP_PA14.LB.1544 <- which(ifelse(PA14_ess_model$PA14.LB.1544!=PA14_ess_model$LB_Model & PA14_ess_model$PA14.LB.1544==0, 1,0)!=0)
FN_PA14.LB.1544 <- which(ifelse(PA14_ess_model$PA14.LB.1544!=PA14_ess_model$LB_Model & PA14_ess_model$PA14.LB.1544==1, 1,0)!=0)

TP_PA14.LB.634 <- which(ifelse(PA14_ess_model$PA14.LB.634==PA14_ess_model$LB_Model & PA14_ess_model$PA14.LB.634==1, 1,0)!=0)
TN_PA14.LB.634 <- which(ifelse(PA14_ess_model$PA14.LB.634==PA14_ess_model$LB_Model & PA14_ess_model$PA14.LB.634==0, 1,0)!=0)
FP_PA14.LB.634 <- which(ifelse(PA14_ess_model$PA14.LB.634!=PA14_ess_model$LB_Model & PA14_ess_model$PA14.LB.634==0, 1,0)!=0)
FN_PA14.LB.634 <- which(ifelse(PA14_ess_model$PA14.LB.634!=PA14_ess_model$LB_Model & PA14_ess_model$PA14.LB.634==1, 1,0)!=0)

TP_PA14.Sputum.510 <- which(ifelse(PA14_ess_model$PA14.Sputum.510==PA14_ess_model$Sputum_Model & PA14_ess_model$PA14.Sputum.510==1, 1,0)!=0)
TN_PA14.Sputum.510 <- which(ifelse(PA14_ess_model$PA14.Sputum.510==PA14_ess_model$Sputum_Model & PA14_ess_model$PA14.Sputum.510==0, 1,0)!=0)
FP_PA14.Sputum.510 <- which(ifelse(PA14_ess_model$PA14.Sputum.510!=PA14_ess_model$Sputum_Model & PA14_ess_model$PA14.Sputum.510==0, 1,0)!=0)
FN_PA14.Sputum.510 <- which(ifelse(PA14_ess_model$PA14.Sputum.510!=PA14_ess_model$Sputum_Model & PA14_ess_model$PA14.Sputum.510==1, 1,0)!=0)

PA14_ess_stats <- data.frame(stat = c(length(TP_PA14.LB.1544),
                                      length(TN_PA14.LB.1544),
                                      length(FP_PA14.LB.1544),
                                      length(FN_PA14.LB.1544),
                                      length(TP_PA14.LB.634),
                                      length(TN_PA14.LB.634),
                                      length(FP_PA14.LB.634),
                                      length(FN_PA14.LB.634),
                                      length(TP_PA14.Sputum.510),
                                      length(TN_PA14.Sputum.510),
                                      length(FP_PA14.Sputum.510),
                                      length(FN_PA14.Sputum.510)),
                             group = c('PA14.LB.1544','PA14.LB.1544','PA14.LB.1544','PA14.LB.1544',
                                       'PA14.LB.634','PA14.LB.634','PA14.LB.634','PA14.LB.634',
                                       'PA14.Sputum.510','PA14.Sputum.510','PA14.Sputum.510','PA14.Sputum.510'),
                             type = c('TP','TN','FP','FN',
                                      'TP','TN','FP','FN',
                                      'TP','TN','FP','FN'))

# combine PA14 statistics and PAO1 statistics into a single dataframe
combined_ess_stats <- rbind(PAO1_ess_stats, PA14_ess_stats)

# calculate accuracy 
acc_PAO1.LB.913 <- (PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.LB.913' & PAO1_ess_stats$type == 'TP'] + PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.LB.913' & PAO1_ess_stats$type == 'TN'])/sum(PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.LB.913'])*100
acc_PAO1.LB.201 <- (PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.LB.201' & PAO1_ess_stats$type == 'TP'] + PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.LB.201' & PAO1_ess_stats$type == 'TN'])/sum(PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.LB.201'])*100
acc_PAO1.LB.335 <- (PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.LB.335' & PAO1_ess_stats$type == 'TP'] + PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.LB.335' & PAO1_ess_stats$type == 'TN'])/sum(PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.LB.335'])*100
acc_PAO1.Sputum.224 <- (PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.Sputum.224' & PAO1_ess_stats$type == 'TP'] + PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.Sputum.224' & PAO1_ess_stats$type == 'TN'])/sum(PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.Sputum.224'])*100
acc_PAO1.Sputum.405 <- (PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.Sputum.405' & PAO1_ess_stats$type == 'TP'] + PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.Sputum.405' & PAO1_ess_stats$type == 'TN'])/sum(PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.Sputum.405'])*100
acc_PAO1.Pyruvate.179 <- (PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.Pyruvate.179' & PAO1_ess_stats$type == 'TP'] + PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.Pyruvate.179' & PAO1_ess_stats$type == 'TN'])/sum(PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.Pyruvate.179'])*100
acc_PAO1.Succinate.640 <- (PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.Succinate.640' & PAO1_ess_stats$type == 'TP'] + PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.Succinate.640' & PAO1_ess_stats$type == 'TN'])/sum(PAO1_ess_stats$stat[PAO1_ess_stats$group == 'PAO1.Succinate.640'])*100

acc_PA14.LB.1544 <- (PA14_ess_stats$stat[PA14_ess_stats$group == 'PA14.LB.1544' & PA14_ess_stats$type == 'TP'] + PA14_ess_stats$stat[PA14_ess_stats$group == 'PA14.LB.1544' & PA14_ess_stats$type == 'TN'])/sum(PA14_ess_stats$stat[PA14_ess_stats$group == 'PA14.LB.1544'])*100
acc_PA14.LB.634 <- (PA14_ess_stats$stat[PA14_ess_stats$group == 'PA14.LB.634' & PA14_ess_stats$type == 'TP'] + PA14_ess_stats$stat[PA14_ess_stats$group == 'PA14.LB.634' & PA14_ess_stats$type == 'TN'])/sum(PA14_ess_stats$stat[PA14_ess_stats$group == 'PA14.LB.634'])*100
acc_PA14.Sputum.510 <- (PA14_ess_stats$stat[PA14_ess_stats$group == 'PA14.Sputum.510' & PA14_ess_stats$type == 'TP'] + PA14_ess_stats$stat[PA14_ess_stats$group == 'PA14.Sputum.510' & PA14_ess_stats$type == 'TN'])/sum(PA14_ess_stats$stat[PA14_ess_stats$group == 'PA14.Sputum.510'])*100

combined_acc <- data.frame(acc_PAO1.LB.913,
                           acc_PAO1.LB.201,
                           acc_PAO1.LB.335,
                           acc_PAO1.Sputum.224,
                           acc_PAO1.Sputum.405,
                           acc_PAO1.Pyruvate.179,
                           acc_PAO1.Succinate.640,
                           acc_PA14.LB.1544,
                           acc_PA14.LB.634,
                           acc_PA14.Sputum.510)
colnames(combined_acc) <- c('PAO1.LB.913',
                            'PAO1.LB.201',
                            'PAO1.LB.335',
                            'PAO1.Sputum.224',
                            'PAO1.Sputum.405',
                            'PAO1.Pyruvate.179',
                            'PAO1.Succinate.640',
                            'PA14.LB.1544',
                            'PA14.LB.634',
                            'PA14.Sputum.510'
                            )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Determine consensus metabolic essential and non-essential genes #####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PA14_HCess_metab_LB <- PA14_ess_model$Genes[which(ifelse(PA14_ess_model$PA14.LB.1544==PA14_ess_model$PA14.LB.634 & PA14_ess_model$PA14.LB.1544==1, 1,0)!=0)]
PA14_HCnoness_metab_LB <- PA14_ess_model$Genes[which(ifelse(PA14_ess_model$PA14.LB.1544==PA14_ess_model$PA14.LB.634 & PA14_ess_model$PA14.LB.1544==0, 1,0)!=0)]

PAO1_HCess_metab_LB <- PAO1_ess_model$Genes[which(ifelse(PAO1_ess_model$PAO1.LB.913==PAO1_ess_model$PAO1.LB.335 & PAO1_ess_model$PAO1.LB.913==PAO1_ess_model$PAO1.LB.201 & PAO1_ess_model$PAO1.LB.913==1, 1,0)!=0)]
PAO1_HCnoness_metab_LB <- PAO1_ess_model$Genes[which(ifelse(PAO1_ess_model$PAO1.LB.913==PAO1_ess_model$PAO1.LB.335 & PAO1_ess_model$PAO1.LB.913==PAO1_ess_model$PAO1.LB.201 & PAO1_ess_model$PAO1.LB.913==0, 1,0)!=0)]

PAO1_HCess_metab_Sputum <- PAO1_ess_model$Genes[which(ifelse(PAO1_ess_model$PAO1.Sputum.224==PAO1_ess_model$PAO1.Sputum.405 & PAO1_ess_model$PAO1.Sputum.224==1, 1,0)!=0)]
PAO1_HCnoness_metab_Sputum <- PAO1_ess_model$Genes[which(ifelse(PAO1_ess_model$PAO1.Sputum.224==PAO1_ess_model$PAO1.Sputum.405 & PAO1_ess_model$PAO1.Sputum.224==0, 1,0)!=0)]

# # output consensus essential/nonessential genes to a .csv file. Note that these files were stitched
# # together to generate the supplementary dataset files.
# write.csv(PA14_HCess_metab_LB, "PA14_consensusEss_metab_LB.csv")
# write.csv(PA14_HCnoness_metab_LB, "PA14_consensusNonEss_metab_LB.csv")
# write.csv(PAO1_HCess_metab_LB, "PAO1_consensusEss_metab_LB.csv")
# write.csv(PAO1_HCnoness_metab_LB, "PAO1_consensusNonEss_metab_LB.csv")
# write.csv(PAO1_HCess_metab_Sputum, "PAO1_consensusEss_metab_Sputum.csv")
# write.csv(PAO1_HCnoness_metab_Sputum, "PAO1_consensusNonEss_metab_Sputum.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Figure 1 - overlap of datasets #####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create dendrogram for both PAO1
PAO1_vitro_ess_dendro <- t(PAO1_vitro_ess[-10])
colnames(PAO1_vitro_ess_dendro) <- PAO1_vitro_ess_dendro[1,]
PAO1_vitro_ess_dendro <- PAO1_vitro_ess_dendro[-1,]
hc_PAO1 <- hclust(dist(PAO1_vitro_ess_dendro, method = "binary")) # uses Jaccard distance
hcdata_PAO1 <- dendro_data(hc_PAO1, type="rectangle")

ggplot(segment(hcdata_PAO1)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data=label(hcdata_PAO1),
            aes(label=label, x=x, y=-1), angle = 90,hjust="right") +
  theme_dendro() +
  ylim(-2,1)

# create dendrogram for both PA14
PA14_vitro_ess_dendro <- t(PA14_vitro_ess[-7])
colnames(PA14_vitro_ess_dendro) <- PA14_vitro_ess_dendro[1,]
PA14_vitro_ess_dendro <- PA14_vitro_ess_dendro[-1,]
hc_PA14 <- hclust(dist(PA14_vitro_ess_dendro, method = "binary")) # uses Jaccard distance
hcdata_PA14 <- dendro_data(hc_PA14, type="rectangle")

ggplot(segment(hcdata_PA14)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data=label(hcdata_PA14),
            aes(label=label, x=x, y=-1), angle = 90,hjust="right") +
  theme_dendro() +
  ylim(-2,1)

# create UpsetR plot to assess overlap for PAO1
upset(PAO1_vitro_ess, sets = rev(c("PAO1.LB.913","PAO1.LB.201","PAO1.LB.335","PAO1.Sputum.224","PAO1.Sputum.405","PAO1.Pyruvate.179","PAO1.Succinate.640")),
      sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on",
      point.size = 2, line.size=0.7,
      keep.order = TRUE,
      mainbar.y.label = "Number of Genes",
      sets.x.label = "Essential Genes per Study",
      query.legend = "top",
      queries = list(list(query = intersects,
                          params = list("PAO1.LB.201","PAO1.LB.335","PAO1.Pyruvate.179","PAO1.Sputum.224","PAO1.Sputum.405","PAO1.Succinate.640","PAO1.LB.913"),
                          color = "orange", active = T,
                          query.name = "All datasets")))

# create UpsetR plot to assess overlap for PA14
upset(PA14_vitro_ess, sets = rev(c("PA14.LB.1544","PA14.LB.634","PA14.Sputum.510")), 
      sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on",
      point.size = 2, line.size=0.7,
      keep.order = TRUE,
      mainbar.y.label = "Number of Genes",
      sets.x.label = "Essential Genes per Study",
      query.legend = "top",
      queries = list(list(query = intersects,
                          params = list("PA14.LB.634","PA14.Sputum.510","PA14.LB.1544"),
                          color = "orange", active = T,
                          query.name = "All datasets")))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Figure 2A - overlap with model ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

target <- c('TN','FN','TP','FP')
combined_ess_stats$type <- reorder.factor(combined_ess_stats$type, new.order = target)
combined_ess_stats <- arrange(combined_ess_stats,type)

ggplot(combined_ess_stats,aes(x = group, y = stat, fill = type)) +
  geom_bar(stat="identity", width = 0.5) +
  ylab("Genes") +
  scale_x_discrete(breaks = c("PA14.LB.634",
                              "PA14.LB.1544",
                              "PA14.Sputum.510",
                              "PAO1.LB.201",
                              "PAO1.LB.335",
                              "PAO1.LB.913",
                              "PAO1.Pyruvate.179",
                              "PAO1.Sputum.224",
                              "PAO1.Sputum.405",
                              "PAO1.Succinate.640"),
                   labels = c("PA14.LB.634",
                              "PA14.LB.1544",
                              "PA14.Sputum.510",
                              "PAO1.LB.201",
                              "PAO1.LB.335",
                              "PAO1.LB.913",
                              "PAO1.Pyruvate.179",
                              "PAO1.Sputum.224",
                              "PAO1.Sputum.405",
                              "PAO1.Succinate.640"),
                   limits = c("PA14.LB.634",
                              "PA14.LB.1544",
                              "PA14.Sputum.510",
                              "PAO1.LB.201",
                              "PAO1.LB.335",
                              "PAO1.LB.913",
                              "PAO1.Sputum.224",
                              "PAO1.Sputum.405",
                              "PAO1.Pyruvate.179",
                              "PAO1.Succinate.640")) +
  scale_fill_manual(breaks = c("FP","TP","FN","TN"),
                    labels = c("mismatch: model essential, screen nonessential",
                               "match: both essential",
                               "mismatch: model nonessential, screen essential",
                               "match: both nonessential"),
                    values = c(myPalette[4],
                               myPalette[6],
                               myPalette[3],
                               myPalette[5])) +
  guides(fill = guide_legend(reverse = TRUE)) + 
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=10),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text = element_text(size=10, colour = "black"))

# stats on size difference between FP and FN

sigTesting_stats <- rbind(filter(PA14_ess_stats, type == "FP" | type == "FN"), filter(PAO1_ess_stats, type == "FP" | type =="FN"))

FP_norm <- shapiro.test(filter(sigTesting_stats,type == "FP")$stat)$p.value
FN_norm <- shapiro.test(filter(sigTesting_stats,type == "FN")$stat)$p.value

wilcox.test(filter(sigTesting_stats,type == "FP")$stat, filter(sigTesting_stats,type == "FN")$stat, paired = TRUE)$p.value
#pval = 0.009765625

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Figure 2B - consensus subsystem assignment ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TP_PA14_all <- intersect(TP_PA14.LB.634,TP_PA14.LB.1544)
FP_PA14_all <- intersect(FP_PA14.LB.634,FP_PA14.LB.1544)
TN_PA14_all <- intersect(TN_PA14.LB.634,TN_PA14.LB.1544)
FN_PA14_all <- intersect(FN_PA14.LB.634,FN_PA14.LB.1544)

TP_PAO1_all <- intersect(intersect(TP_PAO1.LB.201,TP_PAO1.LB.335),TP_PAO1.LB.913)
FP_PAO1_all <- intersect(intersect(FP_PAO1.LB.201,FP_PAO1.LB.335),FP_PAO1.LB.913)
TN_PAO1_all <- intersect(intersect(TN_PAO1.LB.201,TN_PAO1.LB.335),TN_PAO1.LB.913)
FN_PAO1_all <- intersect(intersect(FN_PAO1.LB.201,FN_PAO1.LB.335),FN_PAO1.LB.913)

TP_PA14_all <- as.data.frame(PA14_ess_model$Genes[TP_PA14_all])
FP_PA14_all <- as.data.frame(PA14_ess_model$Genes[FP_PA14_all])
TN_PA14_all <- as.data.frame(PA14_ess_model$Genes[TN_PA14_all])
FN_PA14_all <- as.data.frame(PA14_ess_model$Genes[FN_PA14_all])

TP_PAO1_all <- as.data.frame(PAO1_ess_model$Genes[TP_PAO1_all])
FP_PAO1_all <- as.data.frame(PAO1_ess_model$Genes[FP_PAO1_all])
TN_PAO1_all <- as.data.frame(PAO1_ess_model$Genes[TN_PAO1_all])
FN_PAO1_all <- as.data.frame(PAO1_ess_model$Genes[FN_PAO1_all])

# # output consensus TP,FP,TN,FN genes to a .csv file. Note that these files were stitched together and saved as a .xls file called "highConfidenceGenes_LB_PA14.xls" OR "highConfidenceGenes_LB_PAO1.xls" in order to run the Matlab script for subsystem assignment.
# write.csv(TP_PA14_all, "TP_PA14_all.csv")
# write.csv(FP_PA14_all, "FP_PA14_all.csv")
# write.csv(TN_PA14_all, "TN_PA14_all.csv")
# write.csv(FN_PA14_all, "FN_PA14_all.csv")
# 
# write.csv(TP_PAO1_all, "TP_PAO1_all.csv")
# write.csv(FP_PAO1_all, "FP_PAO1_all.csv")
# write.csv(TN_PAO1_all, "TN_PAO1_all.csv")
# write.csv(FN_PAO1_all, "FN_PAO1_all.csv")

# after collecting this gene information, feed these files into a MATLAB script to get subsystem information for each gene

# after getting subsystem information for each gene and the number of genes_per for each functional category from the MATLAB script, upload this file into R to create a bar graph

gene_per_HC_PA14 <- read.csv('highConfidenceSubSys_PA14.csv')

# remove Exchange from the dataframe
gene_per_HC_PA14 <- gene_per_HC_PA14[-c(5),]

# change the name of the functional categories
gene_per_HC_PA14$Functional_categories <- c('Amino Acid',
                                            'Secondary Metabolites',
                                            'Carbohydrate',
                                            'Energy',
                                            'Glycan',
                                            'Lipid',
                                            'Cofactors and Vitamins',
                                            'Other Amino Acids',
                                            'Terpenoids and Polyketides',
                                            'Nucleotide',
                                            'Other',
                                            'Transport',
                                            'Virulence Factor',
                                            'Xenobiotics')


# add a sum column and sort by this column
gene_per_HC_PA14$sum <- gene_per_HC_PA14$TP_PA14 + gene_per_HC_PA14$TN_PA14 + gene_per_HC_PA14$FN_PA14 + gene_per_HC_PA14$FP_PA14
gene_per_HC_PA14 <- gene_per_HC_PA14[order(gene_per_HC_PA14$sum),]
gene_per_HC_PA14$Functional_categories <- factor(gene_per_HC_PA14$Functional_categories, levels = unique(gene_per_HC_PA14$Functional_categories))

# create a matches dataset
gene_per_HC_PA14_match <- gene_per_HC_PA14[,c(1,2,3)]

# tidy the datasets 
gene_per_HC_PA14_match <- gene_per_HC_PA14_match %>%
  gather(stat, genes_per, TP_PA14, TN_PA14)

# change all the 0 instances to 1 for log plotting and perform log plotting
gene_per_HC_PA14_match[gene_per_HC_PA14_match == 0] <- 1

gene_per_HC_PA14_match$genes_per <- log10(gene_per_HC_PA14_match$genes_per)

# make the bar graph for PA14
ggplot(gene_per_HC_PA14_match, aes(x = Functional_categories, y = genes_per, fill = stat)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge(width=0.5)) +
  ylab('log(Genes)') +
  coord_flip() +
  scale_fill_manual(breaks = c("TP_PA14","TN_PA14"),
                    labels = c("match: both essential",
                               "match: both nonessential"),
                    values = c(myPalette[4],
                               myPalette[3])) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=10),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text = element_text(size=10, colour = "black"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Figure 2E - LB vs Sputum analysis ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# obtain the consensus Sputum essentials
PAO1_HCess_Sputum <- which(ifelse(PAO1_vitro_ess$PAO1.Sputum.224==PAO1_vitro_ess$PAO1.Sputum.405 & PAO1_vitro_ess$PAO1.Sputum.224==1, 1,0)!=0)

# obtain the consensus LB essentials
PAO1_HCess_LB <- which(ifelse(PAO1_vitro_ess$PAO1.LB.913==PAO1_vitro_ess$PAO1.LB.201 & PAO1_vitro_ess$PAO1.LB.913==PAO1_vitro_ess$PAO1.LB.335 & PAO1_vitro_ess$PAO1.LB.913==1, 1,0)!=0)

# find the essentials that are shared between the two screens
PAO1_HCess_shared <- intersect(PAO1_HCess_Sputum,PAO1_HCess_LB)

# find the essentials that are unique to either screen
PAO1_HCess_Sputum_unique <- setdiff(PAO1_HCess_Sputum,PAO1_HCess_LB)
PAO1_HCess_LB_unique <- setdiff(PAO1_HCess_LB,PAO1_HCess_Sputum)

# find the high-confidence essentials that are also predicted essential by the model in LB conditions but NOT Sputum conditions
PAO1_HCess_model_LB_unique <- PAO1_ess_model$Genes[which(ifelse(PAO1_ess_model$PAO1.LB.913==PAO1_ess_model$PAO1.LB.335 & PAO1_ess_model$PAO1.LB.913==PAO1_ess_model$PAO1.LB.201 & PAO1_ess_model$PAO1.LB.913==PAO1_ess_model$LB_Model & PAO1_ess_model$PAO1.LB.913==1 & PAO1_ess_model$Sputum_Model==0, 1,0)!=0)]

# find the high-confidence essentials that are also predicted essential by the model in Sputum conditions but NOT LB conditions
PAO1_HCess_model_Sputum_unique <- PAO1_ess_model$Genes[which(ifelse(PAO1_ess_model$PAO1.Sputum.224==PAO1_ess_model$PAO1.Sputum.405 & PAO1_ess_model$PAO1.Sputum.224==PAO1_ess_model$Sputum_Model & PAO1_ess_model$PAO1.Sputum.224==1 & PAO1_ess_model$LB_Model==0, 1,0)!=0)]

# # output high-confidence PAO1_HCess_model_LB_unique or PAO1_HCess_model_Sputum_unique to a .csv file. These genes were investigated further using the model and MATLAB to determine reasons for their unique essentiality statuts
# write.csv(PAO1_HCess_model_LB_unique,"PAO1_HCess_model_LB_unique.csv")
# write.csv(PAO1_HCess_model_Sputum_unique,"PAO1_HCess_model_Sputum_unique.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table 2 discrepancies analysis ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Find FPs that come up in all the PA14 LB screens
# these are genes that indicate points of curation for future iterations of the model
FP_PA14_all <- intersect(FP_PA14.LB.634,FP_PA14.LB.1544)
FP_PA14_all <- as.data.frame(PA14_ess_model$Genes[FP_PA14_all])

# Find FPs that come up in all the PAO1 LB screens
# these are genes that indicate points of curation for future iterations of the model
FP_PAO1_all <- intersect(intersect(FP_PAO1.LB.201,FP_PAO1.LB.335),FP_PAO1.LB.913)
FP_PAO1_all <- data.frame(PAO1_ess_model$Genes[FP_PAO1_all])

# # output false postive genes for PAO1 and PA14 screens performed on LB to a .csv file. These genes were investigated further using the model and MATLAB to determine reasons for their false essentiality status
# write.csv(FP_PA14_all,"FP_PA14_all.csv")
# write.csv(FP_PAO1_all,"FP_PAO1_all.csv")

# # Compare the predictions of the updated PAO1 model to the consensus PAO1 LB genes

# calculate true positives
length(intersect(PAO1_HCess_metab_LB, PAO1_updated_model_ess$Genes[which(PAO1_updated_model_ess$LB_Model==1)])) #7
#calculate true negatives
PAO1_updated_model_noness <- setdiff(PAO1_model$Genes,PAO1_updated_model_ess$Genes[which(PAO1_updated_model_ess$LB_Model==1)])
length(intersect(PAO1_HCnoness_metab_LB, PAO1_updated_model_noness)) #848
# calculate false positives
length(intersect(PAO1_updated_model_ess$Genes, PAO1_HCnoness_metab_LB)) #15
# calculate false negatives
length(intersect(PAO1_updated_model_noness, PAO1_HCess_metab_LB)) #8

# # Compare the predictions of the updated PAO1 model to the consensus PAO1 Sputum genes

# calculate true positives
length(intersect(PAO1_HCess_metab_Sputum, PAO1_updated_model_ess$Genes[which(PAO1_updated_model_ess$Sputum_Model==1)])) #24
#calculate true negatives
PAO1_updated_model_noness <- setdiff(PAO1_model$Genes,PAO1_updated_model_ess$Genes[which(PAO1_updated_model_ess$Sputum_Model==1)])
length(intersect(PAO1_HCnoness_metab_Sputum, PAO1_updated_model_noness)) #878
# calculate false positives
length(intersect(PAO1_updated_model_ess$Genes, PAO1_HCnoness_metab_Sputum)) #25
# calculate false negatives
length(intersect(PAO1_updated_model_noness, PAO1_HCess_metab_Sputum)) #43

# # Compare the predictions of the updated PA14 model to the consensus PA14 LB genes

# calculate true positives
length(intersect(PA14_HCess_metab_LB, PA14_updated_model_ess$Genes[which(PA14_updated_model_ess$LB_Model==1)])) #45
#calculate true negatives
PA14_updated_model_noness <- setdiff(PA14_model$Genes,PA14_updated_model_ess$Genes[which(PA14_updated_model_ess$LB_Model==1)])
length(intersect(PA14_HCnoness_metab_LB, PA14_updated_model_noness)) #781
# calculate false positives
length(intersect(PA14_updated_model_ess$Genes, PA14_HCnoness_metab_LB)) #19
# calculate false negatives
length(intersect(PA14_updated_model_noness, PA14_HCess_metab_LB)) #68

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Figure 3B - minimal media analysis ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in minMedia_stats
minMedia_stats <- read.csv('model_minimalMediaAnalysis_stats.csv')
minMedia_stats$count <- as.numeric(rownames(minMedia_stats))
minMedia_stats$group <- c(1)

# remove empty first and last rows
minMedia_stats <- minMedia_stats[-c(1,41,42),]

# Average with 95% CI
ggplot(data = minMedia_stats, aes(x = count, y = Ave, group = 1)) + 
  geom_errorbar(aes(ymin=Ave-Std, ymax=Ave+Std), width = 0.1, colour = "black") +
  geom_point(shape = 17, size = 2) +
  xlab("Number of minimal media conditions") + 
  ylab("Number of overlap essential genes") +
  theme_bw() +
  theme(axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text = element_text(size=10, colour = "black"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Figure 4B and 4C - LB media analysis ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in LBMedia_stats
LBMedia_stats <- read.csv('model_LBMediaAnalysis_stats.csv')
LBMedia_stats$count <- as.numeric(rownames(LBMedia_stats))
LBMedia_stats$group <- c(1)

# remove empty first and last rows
LBMedia_stats <- LBMedia_stats[-c(1,22,23,24),]

# read in LBMedia_essGenesLength
LBMedia_length <- read.csv('model_LBMediaAnalysis_essGenesLength.csv')

# remove empty first column
LBMedia_length <- LBMedia_length[,-c(1,22,23)]

# combine all data columns into one column
LBMedia_length <- LBMedia_length %>%
  gather(comparison, length, X2:X21)

# remove X's in front of comparison
LBMedia_length$comparison <- substring(LBMedia_length$comparison,2)

LBMedia_length$comparison <- as.numeric(LBMedia_length$comparison)

LBMedia_length$group <- c(1)

# Average with standard deviation
a <- ggplot(data = LBMedia_stats, aes(x = count, y = Ave, group = 1)) + 
  geom_errorbar(aes(ymin=Ave-Std, ymax=Ave+Std), width = 0.1, colour = "black") +
  geom_point() +
  geom_point(aes(x = count, y = Overlap, group = 1), shape = 17, size = 2) + 
  xlab("Number of LB media components") + 
  ylab("Number of essential genes") +
  theme_bw() +
  theme(axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text = element_text(size=10, colour = "black"))

# read inLBMediaAnalysis_overlapLists
samp1 <- read.csv('model_LBMediaAnalysis_overlapList_sample_1.csv')
samp2 <- read.csv('model_LBMediaAnalysis_overlapList_sample_2.csv')
samp3 <- read.csv('model_LBMediaAnalysis_overlapList_sample_3.csv')
samp4 <- read.csv('model_LBMediaAnalysis_overlapList_sample_4.csv')
samp5 <- read.csv('model_LBMediaAnalysis_overlapList_sample_5.csv')
samp6 <- read.csv('model_LBMediaAnalysis_overlapList_sample_6.csv')
samp7 <- read.csv('model_LBMediaAnalysis_overlapList_sample_7.csv')
samp8 <- read.csv('model_LBMediaAnalysis_overlapList_sample_8.csv')
samp9 <- read.csv('model_LBMediaAnalysis_overlapList_sample_9.csv')
samp10 <- read.csv('model_LBMediaAnalysis_overlapList_sample_10.csv')

samp1 <- samp1[c(-99,-100),]
samp2 <- samp2[c(-99,-100),]
samp3 <- samp3[c(-99,-100),]
samp4 <- samp4[c(-99,-100),]
samp5 <- samp5[c(-99,-100),]
samp6 <- samp6[c(-99,-100),]
samp7 <- samp7[c(-99,-100),]
samp8 <- samp8[c(-99,-100),]
samp9 <- samp9[c(-99,-100),]
samp10 <- samp10[c(-99,-100),]

samp1 <- samp1[,c(-1,-22,-23)]
samp2 <- samp2[,c(-1,-22,-23)]
samp3 <- samp3[,c(-1,-22,-23)]
samp4 <- samp4[,c(-1,-22,-23)]
samp5 <- samp5[,c(-1,-22,-23)]
samp6 <- samp6[,c(-1,-22,-23)]
samp7 <- samp7[,c(-1,-22,-23)]
samp8 <- samp8[,c(-1,-22,-23)]
samp9 <- samp9[,c(-1,-22,-23)]
samp10 <- samp10[,c(-1,-22,-23)]

min1 <- apply(samp1,2,min)
min2 <- apply(samp2,2,min)
min3 <- apply(samp3,2,min)
min4 <- apply(samp4,2,min)
min5 <- apply(samp5,2,min)
min6 <- apply(samp6,2,min)
min7 <- apply(samp7,2,min)
min8 <- apply(samp8,2,min)
min9 <- apply(samp9,2,min)
min10 <- apply(samp10,2,min)

myfunction <- function (val,df,minval){
  minFirst <- min(which(df[,val] == minval[val]))
  return(minFirst)
}

minFirst1 <- sapply(1:20,myfunction,df = samp1, minval = min1)
minFirst2 <- sapply(1:20,myfunction,df = samp2, minval = min2)
minFirst3 <- sapply(1:20,myfunction,df = samp3, minval = min3)
minFirst4 <- sapply(1:20,myfunction,df = samp4, minval = min4)
minFirst5 <- sapply(1:20,myfunction,df = samp5, minval = min5)
minFirst6 <- sapply(1:20,myfunction,df = samp6, minval = min6)
minFirst7 <- sapply(1:20,myfunction,df = samp7, minval = min7)
minFirst8 <- sapply(1:20,myfunction,df = samp8, minval = min8)
minFirst9 <- sapply(1:20,myfunction,df = samp9, minval = min9)
minFirst10 <- sapply(1:20,myfunction,df = samp10, minval = min10)

df_minFirst <- rbind(minFirst1,minFirst2,minFirst3,minFirst4,minFirst5,minFirst6,minFirst7,minFirst8,minFirst9,minFirst10)

df_minFirst_stats <- data.frame(Ave = apply(df_minFirst,2,mean))
df_minFirst_stats$Std <- apply(df_minFirst,2,sd)
df_minFirst_stats$count <- 2:21

b <- ggplot(data = df_minFirst_stats, aes(x = count, y = Ave, group = 1)) + 
  geom_errorbar(aes(ymin=Ave-Std, ymax=Ave+Std), width = 0.1, colour = "black") +
  geom_point(shape = 15) +
  xlab("Number of LB media components") + 
  ylab("Number of replicates to reach overlap") +
  theme_bw() +
  theme(axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text = element_text(size=10, colour = "black"))

multiplot(a,b, cols = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Figure S1- Re-analysis of PAO1 LB Tn-seq data ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#create Venn diagram to assess overlap of original PAO1 LB screens
counts <- function(dataset) {
  venn_counts <- PAO1_vitro_ess
  for (i in 1:length(dataset)) {
    venn_counts <- subset(venn_counts, venn_counts[dataset[i]] == T)
  }
  nrow(venn_counts)
}

plotDatasets <- function(a, ...) {
  grid.newpage()
  if (length(a) == 1) {
    out <- draw.single.venn(counts(a), ...)
  }
  if (length(a) == 2) {
    out <- draw.pairwise.venn(counts(a[1]), counts(a[2]), counts(a[1:2]), ...)
  }
  if (length(a) == 3) {
    out <- draw.triple.venn(counts(a[1]), counts(a[2]), counts(a[3]), counts(a[1:2]), 
                            counts(a[2:3]), counts(a[c(1, 3)]), counts(a), ...)
  }
  if (length(a) == 4) {
    out <- draw.quad.venn(counts(a[1]), counts(a[2]), counts(a[3]), counts(a[4]), 
                          counts(a[1:2]), counts(a[c(1, 3)]), counts(a[c(1, 4)]), counts(a[2:3]), 
                          counts(a[c(2, 4)]), counts(a[3:4]), counts(a[1:3]), counts(a[c(1, 2, 
                                                                                         4)]), counts(a[c(1, 3, 4)]), counts(a[2:4]), counts(a), ...)
  }
  if (!exists("out")) 
    out <- "Oops"
  return(out)
}

plotDatasets(c("PAO1.LB.335","PAO1.LB.201"),
             category=c("PAO1.LB.335", "PAO1.LB.201"),
             lty = "blank",
             fill = c("gray","gray"),
             alpha = 0.5)

#create UpsetR plot to assess overlap for original PAO1 LB screens
upset(PAO1_vitro_ess, sets = rev(c("PAO1.LB.335", "PAO1.LB.201")),
      sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on",
      point.size = 2, line.size=0.7,
      keep.order = TRUE,
      mainbar.y.label = "Number of Genes",
      sets.x.label = "Essential Genes per Study",
      query.legend = "top",
      queries = list(list(query = intersects,
                          params = list("PAO1.LB.335", "PAO1.LB.201"),
                          color = "orange", active = T,
                          query.name = "All datasets")))

# read in files for re-analyzed Tn-seq data for PAO1.LB.201 and PAO1.LB.335
df_335 <- read.table(file = 'PAO1.LB.335_reanalyzed.tsv', sep = '\t', header = TRUE)
df_201 <- read.table(file = 'PAO1.LB.201_reanalyzed.tsv', sep = '\t', header = TRUE)

# filter dataframes for essential genes
ess_335 <- filter(df_335,Essentiality=="Reduced" & padj<0.05 & Uncertainty < 0.1)
ess_201 <- filter(df_201,Essentiality=="Reduced" & padj<0.05 & Uncertainty < 0.1)

# create a dataframe that combines the two essential dataframes
GenesEss <- union(ess_335$id, ess_201$id)
essentials <- data.frame(Genes=GenesEss, PAO1.LB.335 = NA, PAO1.LB.201 = NA)
essentials$PAO1.LB.335[which(essentials$Genes %in% ess_335$id == TRUE)] <- 1
essentials$PAO1.LB.201[which(essentials$Genes %in% ess_201$id == TRUE)] <- 1
essentials$PAO1.LB.335[is.na(essentials$PAO1.LB.335) == TRUE] <- 0
essentials$PAO1.LB.201[is.na(essentials$PAO1.LB.201) == TRUE] <- 0

#create Venn diagram to assess overlap of re-analyzed PAO1 LB screens
counts <- function(dataset) {
  venn_counts <- essentials
  for (i in 1:length(dataset)) {
    venn_counts <- subset(venn_counts, venn_counts[dataset[i]] == T)
  }
  nrow(venn_counts)
}

plotDatasets <- function(a, ...) {
  grid.newpage()
  if (length(a) == 1) {
    out <- draw.single.venn(counts(a), ...)
  }
  if (length(a) == 2) {
    out <- draw.pairwise.venn(counts(a[1]), counts(a[2]), counts(a[1:2]), ...)
  }
  if (length(a) == 3) {
    out <- draw.triple.venn(counts(a[1]), counts(a[2]), counts(a[3]), counts(a[1:2]), 
                            counts(a[2:3]), counts(a[c(1, 3)]), counts(a), ...)
  }
  if (length(a) == 4) {
    out <- draw.quad.venn(counts(a[1]), counts(a[2]), counts(a[3]), counts(a[4]), 
                          counts(a[1:2]), counts(a[c(1, 3)]), counts(a[c(1, 4)]), counts(a[2:3]), 
                          counts(a[c(2, 4)]), counts(a[3:4]), counts(a[1:3]), counts(a[c(1, 2, 
                                                                                         4)]), counts(a[c(1, 3, 4)]), counts(a[2:4]), counts(a), ...)
  }
  if (!exists("out")) 
    out <- "Oops"
  return(out)
}

plotDatasets(c("PAO1.LB.335","PAO1.LB.201"),
             category=c("Re-analyzed PAO1.LB.335","Re-analyzed PAO1.LB.201"),
             lty = "blank",
             fill = c("gray","gray"),
             alpha = 0.5)

#create UpsetR plot to assess overlap for re-analyzed PAO1 LB screens
upset(essentials, sets = rev(c("PAO1.LB.335","PAO1.LB.201")),
      sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on",
      point.size = 2, line.size=0.7,
      keep.order = TRUE,
      mainbar.y.label = "Number of Genes",
      sets.x.label = "Essential Genes per Study",
      query.legend = "top",
      queries = list(list(query = intersects,
                          params = list("PAO1.LB.335","PAO1.LB.201"),
                          color = "orange", active = T,
                          query.name = "All datasets")))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Figure S2- UpsetR for LB ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#create Venn diagram to assess overlap of PAO1 LB screens
counts <- function(dataset) {
  venn_counts <- PAO1_vitro_ess
  for (i in 1:length(dataset)) {
    venn_counts <- subset(venn_counts, venn_counts[dataset[i]] == T)
  }
  nrow(venn_counts)
}

plotDatasets <- function(a, ...) {
  grid.newpage()
  if (length(a) == 1) {
    out <- draw.single.venn(counts(a), ...)
  }
  if (length(a) == 2) {
    out <- draw.pairwise.venn(counts(a[1]), counts(a[2]), counts(a[1:2]), ...)
  }
  if (length(a) == 3) {
    out <- draw.triple.venn(counts(a[1]), counts(a[2]), counts(a[3]), counts(a[1:2]), 
                            counts(a[2:3]), counts(a[c(1, 3)]), counts(a), ...)
  }
  if (length(a) == 4) {
    out <- draw.quad.venn(counts(a[1]), counts(a[2]), counts(a[3]), counts(a[4]), 
                          counts(a[1:2]), counts(a[c(1, 3)]), counts(a[c(1, 4)]), counts(a[2:3]), 
                          counts(a[c(2, 4)]), counts(a[3:4]), counts(a[1:3]), counts(a[c(1, 2, 
                                                                                         4)]), counts(a[c(1, 3, 4)]), counts(a[2:4]), counts(a), ...)
  }
  if (!exists("out")) 
    out <- "Oops"
  return(out)
}

plotDatasets(c("PAO1.LB.913","PAO1.LB.201","PAO1.LB.335"),
             category=c("PAO1.LB.913","PAO1.LB.201","PAO1.LB.335"),
             lty = "blank",
             fill = c("gray","gray","gray"),
             alpha = 0.5)

#create Venn diagram to assess overlap of PA14 LB screens
counts <- function(dataset) {
  venn_counts <- PA14_vitro_ess
  for (i in 1:length(dataset)) {
    venn_counts <- subset(venn_counts, venn_counts[dataset[i]] == T)
  }
  nrow(venn_counts)
}

plotDatasets <- function(a, ...) {
  grid.newpage()
  if (length(a) == 1) {
    out <- draw.single.venn(counts(a), ...)
  }
  if (length(a) == 2) {
    out <- draw.pairwise.venn(counts(a[1]), counts(a[2]), counts(a[1:2]), ...)
  }
  if (length(a) == 3) {
    out <- draw.triple.venn(counts(a[1]), counts(a[2]), counts(a[3]), counts(a[1:2]), 
                            counts(a[2:3]), counts(a[c(1, 3)]), counts(a), ...)
  }
  if (length(a) == 4) {
    out <- draw.quad.venn(counts(a[1]), counts(a[2]), counts(a[3]), counts(a[4]), 
                          counts(a[1:2]), counts(a[c(1, 3)]), counts(a[c(1, 4)]), counts(a[2:3]), 
                          counts(a[c(2, 4)]), counts(a[3:4]), counts(a[1:3]), counts(a[c(1, 2, 
                                                                                         4)]), counts(a[c(1, 3, 4)]), counts(a[2:4]), counts(a), ...)
  }
  if (!exists("out")) 
    out <- "Oops"
  return(out)
}

plotDatasets(c("PA14.LB.1544","PA14.LB.634"),
             category=c("PA14.LB.1544","PA14.LB.634"),
             lty = "blank",
             fill = c("gray","gray"),
             alpha = 0.5)

#create UpsetR plot to assess overlap for PAO1 LB screens
upset(PAO1_vitro_ess, sets = rev(c("PAO1.LB.913","PAO1.LB.201","PAO1.LB.335")),
      sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on",
      point.size = 2, line.size=0.7,
      keep.order = TRUE,
      mainbar.y.label = "Number of Genes",
      sets.x.label = "Essential Genes per Study",
      query.legend = "top",
      queries = list(list(query = intersects,
                          params = list("PAO1.LB.201","PAO1.LB.335","PAO1.LB.913"),
                          color = "orange", active = T,
                          query.name = "All datasets")))

#create UpsetR plot to assess overlap for PA14 LB screens
upset(PA14_vitro_ess, sets = rev(c("PA14.LB.1544","PA14.LB.634")), 
      sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on",
      point.size = 2, line.size=0.7,
      keep.order = TRUE,
      mainbar.y.label = "Number of Genes",
      sets.x.label = "Essential Genes per Study",
      query.legend = "top",
      queries = list(list(query = intersects,
                          params = list("PA14.LB.634","PA14.LB.1544"),
                          color = "orange", active = T,
                          query.name = "All datasets")))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Figure S3- KEGG functional categories for PAO1 ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# after getting subsystem information for each gene and the number of genes_per for each functional category from the MATLAB script, upload this file into R to create a bar graph

gene_per_HC_PAO1 <- read.csv('highConfidenceSubSys_PAO1.csv')

# remove Exchange from the dataframe
gene_per_HC_PAO1 <- gene_per_HC_PAO1[-c(5),]

# change the name of the functional categories
gene_per_HC_PAO1$Functional_categories <- c('Amino Acid',
                                            'Secondary Metabolites',
                                            'Carbohydrate',
                                            'Energy',
                                            'Glycan',
                                            'Lipid',
                                            'Cofactors and Vitamins',
                                            'Other Amino Acids',
                                            'Terpenoids and Polyketides',
                                            'Nucleotide',
                                            'Other',
                                            'Transport',
                                            'Virulence Factor',
                                            'Xenobiotics')

# add a sum column and sort by this column
gene_per_HC_PAO1$sum <- gene_per_HC_PAO1$TP_PAO1 + gene_per_HC_PAO1$TN_PAO1 + gene_per_HC_PAO1$FN_PAO1 + gene_per_HC_PAO1$FP_PAO1
gene_per_HC_PAO1 <- gene_per_HC_PAO1[order(gene_per_HC_PAO1$sum),]
gene_per_HC_PAO1$Functional_categories <- factor(gene_per_HC_PAO1$Functional_categories, levels = unique(gene_per_HC_PAO1$Functional_categories))

# create a matches dataset
gene_per_HC_PAO1_match <- gene_per_HC_PAO1[,c(1,2,3)]

# tidy the datasets 
gene_per_HC_PAO1_match <- gene_per_HC_PAO1_match %>%
  gather(stat, genes_per, TP_PAO1, TN_PAO1)

# change all the 0 instances to 1 for log plotting and perform log plotting
gene_per_HC_PAO1_match[gene_per_HC_PAO1_match == 0] <- 1

gene_per_HC_PAO1_match$genes_per <- log10(gene_per_HC_PAO1_match$genes_per)

# make the bar graph for PAO1
ggplot(gene_per_HC_PAO1_match, aes(x = Functional_categories, y = genes_per, fill = stat)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab('log(Genes)') +
  coord_flip() +
  scale_fill_manual(breaks = c("TP_PAO1","TN_PAO1"),
                    labels = c("match: both essential",
                               "match: both nonessential"),
                    values = c(myPalette[4],
                               myPalette[3])) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=10),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text = element_text(size=10, colour = "black"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table S2 - Model accuracy stats ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TableS2 <- t(combined_acc)
colnames(TableS2) <- c("% Accuracy")

# # output accuracy metrics to a .csv file. 
# write.csv(TableS2,"Table_S2.csv")