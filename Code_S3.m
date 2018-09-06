%% set-up the environment
clear all
close all
clc
initCobraToolbox()
changeCobraSolver('gurobi5');

%% load models

PA14=load('pa14.mat');
PA14 = PA14.PA14;
PA14 = changeMedia_SEED(PA14, 1,'');
PA14 = changeObjective(PA14,'PA14_Biomass',1);
PA14_sol = optimizeCbModel(PA14) % confirm value is ~15.7298

PAO1=load('pao1.mat');
PAO1 = PAO1.PAO1;
PAO1 = changeMedia_SEED(PAO1, 1,'');
PAO1 = changeObjective(PAO1,'PAO1_Biomass',1);
PAO1_sol = optimizeCbModel(PAO1) % confirm value is ~15.7298

%% load high-confidence genes

% load in high confidence genes
[num,PA14_highConfidence,raw] = xlsread('highConfidenceGenes_LB_PA14.xlsx',1);
[num,PAO1_highConfidence,raw] = xlsread('highConfidenceGenes_LB_PAO1.xlsx',1);

% separate out by stat
TP_PA14 = PA14_highConfidence(2:end,1);
TN_PA14 = PA14_highConfidence(2:end,2);
FP_PA14 = PA14_highConfidence(2:end,3);
FN_PA14 = PA14_highConfidence(2:end,4);

TP_PAO1 = PAO1_highConfidence(2:end,1);
TN_PAO1 = PAO1_highConfidence(2:end,2);
FP_PAO1 = PAO1_highConfidence(2:end,3);
FN_PAO1 = PAO1_highConfidence(2:end,4);

% remove empty values
TP_PA14 = TP_PA14(~cellfun('isempty',TP_PA14));
TN_PA14 = TN_PA14(~cellfun('isempty',TN_PA14));
FP_PA14 = FP_PA14(~cellfun('isempty',FP_PA14));
FN_PA14 = FN_PA14(~cellfun('isempty',FN_PA14));

TP_PAO1 = TP_PAO1(~cellfun('isempty',TP_PAO1));
TN_PAO1 = TN_PAO1(~cellfun('isempty',TN_PAO1));
FP_PAO1 = FP_PAO1(~cellfun('isempty',FP_PAO1));
FN_PAO1 = FN_PAO1(~cellfun('isempty',FN_PAO1));

%% genes_per calculation

% re-assign subsystems so that they match the broadly categorized KEGG
% subsystems
[pao1_modnew,pao1_subsys,pao1_matches,pao1_missing,pao1_multiple] = broadSubsys(PAO1,3);
[pa14_modnew,pa14_subsys,pa14_matches,pa14_missing,pa14_multiple] = broadSubsys(PA14,3);
% determine the number of unique subsystems (i.e., 'functional categories'
% there are (should result in 15)
unique_subsys = unique(pao1_subsys);

% initialize variables
genes_per_TP_PA14 = zeros(length(unique_subsys),1);
genes_per_TN_PA14 = zeros(length(unique_subsys),1);
genes_per_FP_PA14 = zeros(length(unique_subsys),1);
genes_per_FN_PA14 = zeros(length(unique_subsys),1);

genes_per_TP_PAO1 = zeros(length(unique_subsys),1);
genes_per_TN_PAO1 = zeros(length(unique_subsys),1);
genes_per_FP_PAO1 = zeros(length(unique_subsys),1);
genes_per_FN_PAO1 = zeros(length(unique_subsys),1);

for i = 1:length(pa14_modnew.genes)
     activeRxns_sparse = pa14_modnew.rxnGeneMat(:,i);
     activeRxns = find(activeRxns_sparse);
     activeRxns_subSysBroad = pa14_modnew.subSysBroad(activeRxns(1));
     if isempty(find(strcmp(TP_PA14,pa14_modnew.genes(i)))) ~= 1
         genes_per_TP_PA14(find(strcmp(unique_subsys,activeRxns_subSysBroad))) = genes_per_TP_PA14(find(strcmp(unique_subsys,activeRxns_subSysBroad))) + 1;
     elseif isempty(find(strcmp(TN_PA14,pa14_modnew.genes(i)))) ~= 1
         genes_per_TN_PA14(find(strcmp(unique_subsys,activeRxns_subSysBroad))) = genes_per_TN_PA14(find(strcmp(unique_subsys,activeRxns_subSysBroad))) + 1;
     elseif isempty(find(strcmp(FP_PA14,pa14_modnew.genes(i)))) ~= 1
         genes_per_FP_PA14(find(strcmp(unique_subsys,activeRxns_subSysBroad))) = genes_per_FP_PA14(find(strcmp(unique_subsys,activeRxns_subSysBroad))) + 1;
     else isempty(find(strcmp(FN_PA14,pa14_modnew.genes(i)))) ~= 1
         genes_per_FN_PA14(find(strcmp(unique_subsys,activeRxns_subSysBroad))) = genes_per_FN_PA14(find(strcmp(unique_subsys,activeRxns_subSysBroad))) + 1;
     end
end

for i = 1:length(pao1_modnew.genes)
     activeRxns_sparse = pao1_modnew.rxnGeneMat(:,i);
     activeRxns = find(activeRxns_sparse);
     activeRxns_subSysBroad = pao1_modnew.subSysBroad(activeRxns(1));
     if isempty(find(strcmp(TP_PAO1,pao1_modnew.genes(i)))) ~= 1
         genes_per_TP_PAO1(find(strcmp(unique_subsys,activeRxns_subSysBroad))) = genes_per_TP_PAO1(find(strcmp(unique_subsys,activeRxns_subSysBroad))) + 1;
     elseif isempty(find(strcmp(TN_PAO1,pao1_modnew.genes(i)))) ~= 1
         genes_per_TN_PAO1(find(strcmp(unique_subsys,activeRxns_subSysBroad))) = genes_per_TN_PAO1(find(strcmp(unique_subsys,activeRxns_subSysBroad))) + 1;
     elseif isempty(find(strcmp(FP_PAO1,pao1_modnew.genes(i)))) ~= 1
         genes_per_FP_PAO1(find(strcmp(unique_subsys,activeRxns_subSysBroad))) = genes_per_FP_PAO1(find(strcmp(unique_subsys,activeRxns_subSysBroad))) + 1;
     else isempty(find(strcmp(FN_PAO1,pao1_modnew.genes(i)))) ~= 1
         genes_per_FN_PAO1(find(strcmp(unique_subsys,activeRxns_subSysBroad))) = genes_per_FN_PAO1(find(strcmp(unique_subsys,activeRxns_subSysBroad))) + 1;
     end
end

%% output the Subsys results

xlswrite('highConfidenceSubSys_PAO1.xlsx', unique_subsys, 'Sheet1','A1');
xlswrite('highConfidenceSubSys_PAO1.xlsx', genes_per_TP_PAO1, 'Sheet1','B1');
xlswrite('highConfidenceSubSys_PAO1.xlsx', genes_per_TN_PAO1, 'Sheet1','C1');
xlswrite('highConfidenceSubSys_PAO1.xlsx', genes_per_FP_PAO1, 'Sheet1','D1');
xlswrite('highConfidenceSubSys_PAO1.xlsx', genes_per_FN_PAO1, 'Sheet1','E1');

xlswrite('highConfidenceSubSys_PA14_BHIupdate.xlsx', unique_subsys, 'Sheet1','A1');
xlswrite('highConfidenceSubSys_PA14_BHIupdate.xlsx', genes_per_TP_PA14, 'Sheet1','B1');
xlswrite('highConfidenceSubSys_PA14_BHIupdate.xlsx', genes_per_TN_PA14, 'Sheet1','C1');
xlswrite('highConfidenceSubSys_PA14_BHIupdate.xlsx', genes_per_FP_PA14, 'Sheet1','D1');
xlswrite('highConfidenceSubSys_PA14_BHIupdate.xlsx', genes_per_FN_PA14, 'Sheet1','E1');

%% biomass precursor check

essMissingMets_TP_PAO1 = cell(59,length(TP_PAO1)); %59 is the number of metabolites in biomass reaction?
essMissingMetsNames_TP_PAO1 = cell(59,length(TP_PAO1)); %59 is the number of metabolites in biomass reaction?

for i = 1:length(TP_PAO1)
    clear missingMets_names
    PAO1_deleted = deleteModelGenes(PAO1,TP_PAO1(i));
    [missingMets,presentMets] = biomassPrecursorCheck(PAO1_deleted);
    essMissingMets_TP_PAO1(1:length(missingMets),i) = missingMets;
    for j = 1:length(missingMets)
        metIDX = find(strcmp(missingMets(j),PAO1.mets));
        missingMets_names(j) = PAO1.metNames(metIDX);
    end
    missingMets_names = missingMets_names';
    essMissingMetsNames_TP_PAO1(1:length(missingMets_names),i) = missingMets_names;
end

essMissingMets_TP_PA14 = cell(59,length(TP_PA14)); %59 is the number of metabolites in biomass reaction?
essMissingMetsNames_TP_PA14 = cell(59,length(TP_PA14)); %59 is the number of metabolites in biomass reaction?

for i = 1:length(TP_PA14)
    clear missingMets_names
    PA14_deleted = deleteModelGenes(PA14,TP_PA14(i));
    [missingMets,presentMets] = biomassPrecursorCheck(PA14_deleted);
    essMissingMets_TP_PA14(1:length(missingMets),i) = missingMets;
    for j = 1:length(missingMets)
        metIDX = find(strcmp(missingMets(j),PA14.mets));
        missingMets_names(j) = PA14.metNames(metIDX);
    end
    missingMets_names = missingMets_names';
    essMissingMetsNames_TP_PA14(1:length(missingMets_names),i) = missingMets_names;
end

xlswrite('pa14biomass.xlsx',essMissingMetsNames_TP_PA14);
