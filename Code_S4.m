%% set-up the environment
clear all
close all
clc
initCobraToolbox()
changeCobraSolver('gurobi5');

%% load models

PAO1=load('PAO1.mat');
PAO1 = PAO1.PAO1;
PAO1 = changeMedia_SEED(PAO1, 1,'');
PAO1 = changeObjective(PAO1,'PAO1_Biomass',1);
PAO1_sol = optimizeCbModel(PAO1) % confirm value is ~15.7298

PAO1_LB = changeMedia_SEED(PAO1, 1,''); % create an LB environment model
PAO1_LB = changeObjective(PAO1_LB,'PAO1_Biomass',1);
PAO1_LB_sol = optimizeCbModel(PAO1_LB)

PAO1_scfm = changeMedia_SEED(PAO1, 2,''); % create a Sputum environment model
PAO1_scfm = changeObjective(PAO1_scfm,'PAO1_Biomass',1);
PAO1_scfm_sol = optimizeCbModel(PAO1_scfm)

%% Determine biomass components that cannot be synthesized by sputum 
% unique essential genes

[PAO1_3527, hasEffect] = deleteModelGenes(PAO1_scfm,'PA3527');
[missingMets_3527, presentMets] = biomassPrecursorCheck(PAO1_3527);

[PAO1_3050, hasEffect] = deleteModelGenes(PAO1_scfm,'PA3050');
[missingMets_3050, presentMets] = biomassPrecursorCheck(PAO1_3050);

[PAO1_2876, hasEffect] = deleteModelGenes(PAO1_scfm,'PA2876');
[missingMets_2876, presentMets] = biomassPrecursorCheck(PAO1_2876);

[PAO1_0402, hasEffect] = deleteModelGenes(PAO1_scfm,'PA0402');
[missingMets_0402, presentMets] = biomassPrecursorCheck(PAO1_0402);

missingMets_3527
missingMets_3050
missingMets_2876
missingMets_0402

%missing metabolites for all four knockouts:
% cpd00046
% cpd00091
% cpd00298
% cpd00206
%all four of these are involved in pyrimidine metabolism

%% flux sampling

% Note, in order to run this next chunk, I needed to change the directory 
% to optGpSampler_1.1_Matlab

flux_LB = optGpSampler(PAO1_LB,[],3000,1,1,'gurobi',0);
flux_scfm = optGpSampler(PAO1_scfm,[],3000,1,1,'gurobi',0);

% pyrimidine metabolism reactions
rxn_IDs = {
    'rxn00414';
    'rxn01018';
    'rxn01465';
    'rxn01360';
    'rxn01362';
    'rxn00710';
    'rxn00711';
    'rxn00717';
    'rxn00117';
    'rxn00119';
    'rxn00708';
    'rxn00412';
    'rxn00409';
    'rxn00364';
    'rxn00363';
    'rxn06076';
    'rxn01673';
    'rxn01219';
    'rxn01218';
    'rxn01672';
    'rxn01678';
    'rxn01517';
    'rxn01521';
    'rxn01520';
    'rxn01512';
    'rxn01513';
    'rxn01145'};

fluxDist_LB = zeros(length(rxn_IDs),3000);
fluxDist_scfm = zeros(length(rxn_IDs),3000);

fluxMean_LB = zeros(length(rxn_IDs),1);
fluxMean_scfm = zeros(length(rxn_IDs),1);

fluxMedian_LB = zeros(length(rxn_IDs),1);
fluxMedian_scfm = zeros(length(rxn_IDs),1);

for i = 1:length(rxn_IDs)
    fluxDist_LB(i,:) = flux_LB.points(findRxnIDs(PAO1_LB,rxn_IDs(i)),:);
    fluxDist_scfm(i,:) = flux_scfm.points(findRxnIDs(PAO1_LB,rxn_IDs(i)),:);
    
    fluxMean_LB(i) = mean(fluxDist_LB(i,:));
    fluxMean_scfm(i) = mean(fluxDist_scfm(i,:));
    
    fluxMedian_LB(i) = median(fluxDist_LB(i,:));
    fluxMedian_scfm(i) = median(fluxDist_scfm(i,:));
end  

%% output results

% Note, remember to change directory back to the folder with the other
% supplementary information files

xlswrite('SPUTUMvsLB_fluxSampling.xlsx', cellstr('Reaction IDs'), 'Sheet1','A1');
xlswrite('SPUTUMvsLB_fluxSampling.xlsx', rxn_IDs, 'Sheet1','A2');
xlswrite('SPUTUMvsLB_fluxSampling.xlsx', cellstr('Flux Mean: LB'), 'Sheet1','B1');
xlswrite('SPUTUMvsLB_fluxSampling.xlsx', fluxMean_LB, 'Sheet1','B2');
xlswrite('SPUTUMvsLB_fluxSampling.xlsx', cellstr('Flux Mean: Sputum'), 'Sheet1','C1');
xlswrite('SPUTUMvsLB_fluxSampling.xlsx', fluxMean_scfm, 'Sheet1','C2');
xlswrite('SPUTUMvsLB_fluxSampling.xlsx', cellstr('Flux Median: LB'), 'Sheet1','D1');
xlswrite('SPUTUMvsLB_fluxSampling.xlsx', fluxMedian_LB, 'Sheet1','D2');
xlswrite('SPUTUMvsLB_fluxSampling.xlsx', cellstr('Flux Median: Sputum'), 'Sheet1','E1');
xlswrite('SPUTUMvsLB_fluxSampling.xlsx', fluxMedian_scfm, 'Sheet1','E2');