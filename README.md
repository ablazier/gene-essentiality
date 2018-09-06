# gene-essentiality

Code for "Reconciling high-throughput gene essentiality data with metabolic network reconstructions" by Blazier and Papin (2018). Currently in submission.

# Contents:
Code_S1.R - R code and data to generate paper figures

Code_S2.m- MATLAB code to perform gene essentiality predictions

Code_S3.m - MATLAB code to assign subsystems to high-confidence metabolic essential and non-essential genes (Figure 2B)

Code_S4.m - MATLAB code to perform Sputum vs LB essentiality comparison (Figure 2E)

Code_S5.m - MATLAB code to perform minimal media analysis (Figure 3)

Code_S6.m - MATLAB code to perform LB media analysis (Figure 4)

Code_S7.m - MATLAB code to perform gene essentiality predictions on updated models (Table S3)

Dataset_S1.csv - PAO1 candidate essential genes for in vitro screens. Called by Code_S1.

Dataset_S2.csv - PA14 candidate essential genes for in vitro screens. Called by Code_S1.

Dataset_S3.csv - PAO1 model predicted essential genes for in silico screens. Called by Code_S1.

Dataset_S4.csv - PA14 model predicted essential genes for in silico screens. Called by Code_S1.

Dataset_S10.csv - PAO1 candidate essential genes for in vitro screens for the updated PAO1 model. Called by Code_S1.

Dataset_S11.csv - PA14 candidate essential genes for in vitro screens for the updated PA14 model. Called by Code_S1.

PAO1.LB.201_reanalyzed.tsv - Results from re-analysis of PAO1.LB.201 dataset. Called by Code_S1.

PAO1.LB.335_reanalyzed.tsv - Results from re-analysis of PAO1.LB.335 dataset. Called by Code_S1.

addExchangeRxn_JB.m - Function to add an exchange reaction to the models. Called by Code_S5.m*

broadSubsys.m - Function to evaluate the functional subsystems of the models. Called by Code_S3.*

changeMedia_SEED.m - Function to change the media conditions of the models. Called by Code_S2, Code_S3, Code_S4, Code_S5, Code_S6, Code_S7.*

highConfidenceGenes_LB_PA14.xlsx - PA14 high-confidence metabolic essential/non-essential genes compared to model predictions. Called by Code_S3.

highConfidenceGenes_LB_PAO1.xlsx - PAO1 high-confidence metabolic essential/non-essential genes compared to model predictions. Called by Code_S3.

highConfidenceSubSys_PA14.csv - Subsystem assignment for PA14 high-confidence metabolic essential/non-essential genes compared to model predictions. Called by Code_S1.

highConfidenceSubSys_PAO1.csv - Subsystem assignment for PAO1 high-confidence metabolic essential/non-essential genes compared to model predictions. Called by Code_S1.

model_LBMediaAnalysis_essGenesLength.csv - Number of essential genes for 100 simulations across 21 LB formulations. Called by Code_S1.

model_LBMediaAnalysis_overlapList_sample1.csv - One run assessing the number the shared essential genes over 100 simulations across 21 LB formulations. Called by Code_S1.

model_LBMediaAnalysis_overlapList_sample2.csv - One run assessing the number the shared essential genes over 100 simulations across 21 LB formulations. Called by Code_S1.

model_LBMediaAnalysis_overlapList_sample3.csv - One run assessing the number the shared essential genes over 100 simulations across 21 LB formulations. Called by Code_S1.

model_LBMediaAnalysis_overlapList_sample4.csv - One run assessing the number the shared essential genes over 100 simulations across 21 LB formulations. Called by Code_S1.

model_LBMediaAnalysis_overlapList_sample5.csv - One run assessing the number the shared essential genes over 100 simulations across 21 LB formulations. Called by Code_S1.

model_LBMediaAnalysis_overlapList_sample6.csv - One run assessing the number the shared essential genes over 100 simulations across 21 LB formulations. Called by Code_S1.

model_LBMediaAnalysis_overlapList_sample7.csv - One run assessing the number the shared essential genes over 100 simulations across 21 LB formulations. Called by Code_S1.

model_LBMediaAnalysis_overlapList_sample8.csv - One run assessing the number the shared essential genes over 100 simulations across 21 LB formulations. Called by Code_S1.

model_LBMediaAnalysis_overlapList_sample9.csv - One run assessing the number the shared essential genes over 100 simulations across 21 LB formulations. Called by Code_S1.

model_LBMediaAnalysis_overlapList_sample10.csv - One run assessing the number the shared essential genes over 100 simulations across 21 LB formulations. Called by Code_S1.

model_LBMediaAnalysis_stats.csv - Statistical information for essential genes over 100 simulations across 21 LB formulations. Called by Code_S1.

model_minimalMediaAnalysis_stats.csv - Statistical information for essential genes over 500 simulations across 40 combinations of minimal media conditions. Called by Code_S1.

modelGenes_PA14.csv - PA14 genes in the model. Called by Code_S1.

modelGenes_PAO1.csv - PAO1 genes in the model. Called by Code_S1.

multiplot.R - Function to plot multiple panels. Called by Code_S1. *

pa14.mat - Workspace containing the PA14 model (Bartell, Blazier et al., 2017). The model can also be downloaded from: http://bme.virginia.edu/csbl/Downloads1.html. Called by Code_S2, Code_S3, Code_S5, and Code_S6.

pao1.mat - Workspace containing the PAO1 model (Bartell, Blazier et al., 2017). The model can also be downloaded from: http://bme.virginia.edu/csbl/Downloads1.html. Called by Code_S2, Code_S3, and Code_S4.

pa14_updated.mat - Workspace containing the updated PA14 model. Called by Code_S7.

pao1.mat - Workspace containing the PAO1 model (Bartell, Blazier et al., 2017). The model can also be downloaded from: http://bme.virginia.edu/csbl/Downloads1.html. Called by Code_S2, Code_S3, and Code_S4.

pao1_updated.mat - Workspace containing the updated PAO1 model. Called by Code_S7.

*Note. Not all functions were originally written by Blazier and Papin.

# Software information:
MATLAB R2016a

COBRA Toolbox 2.0.5

Gurobi 6.5 solver

optGPSampler1.1

R 3.3.3
# R-code instructions:
Download the entire supplementary information contents into a folder
The following packages will need to be installed:

tidyverse

ggdendro

gdata

UpSetR

VennDiagram

In R-studio, create a new R-project within the above folder
# MATLAB-code instructions:
Download the entire supplementary information contents into a folder

In MATLAB, change the working directory to the above folder
