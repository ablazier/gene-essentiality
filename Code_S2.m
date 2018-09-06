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

%% LB essentiality

PA14 = changeMedia_SEED(PA14, 1,'');
PA14_sol = optimizeCbModel(PA14)

PAO1 = changeMedia_SEED(PAO1, 1,'');
PAO1_sol = optimizeCbModel(PAO1)

essIDX_PA14 = [];
essGenes_PA14 = [];

for i = 1:length(PA14.genes)
    PA14_deleted = deleteModelGenes(PA14,PA14.genes(i));
    sol = optimizeCbModel(PA14_deleted);
    if sol.f < 0.0001
        essIDX_PA14 = [essIDX_PA14,i];
        essGenes_PA14 = [essGenes_PA14,PA14.genes(i)];
    end
end

essMissingMets_PA14 = cell(59,length(essGenes_PA14)); %59 is the number of metabolites in biomass reaction?

for i = 1:length(essGenes_PA14)
    PA14_deleted = deleteModelGenes(PA14,essGenes_PA14(i));
    [missingMets,presentMets] = biomassPrecursorCheck(PA14_deleted);
    essMissingMets_PA14(1:length(missingMets),i) = missingMets;
end

essIDX_PAO1 = [];
essGenes_PAO1 = [];

for i = 1:length(PAO1.genes)
    PAO1_deleted = deleteModelGenes(PAO1,PAO1.genes(i));
    sol = optimizeCbModel(PAO1_deleted);
    if sol.f < 0.0001
        essIDX_PAO1 = [essIDX_PAO1,i];
        essGenes_PAO1 = [essGenes_PAO1,PAO1.genes(i)];
    end
end

essMissingMets_PAO1 = cell(59,length(essGenes_PAO1)); %59 is the number of metabolites in biomass reaction?

for i = 1:length(essGenes_PAO1)
    PAO1_deleted = deleteModelGenes(PAO1,essGenes_PAO1(i));
    [missingMets,presentMets] = biomassPrecursorCheck(PAO1_deleted);
    essMissingMets_PAO1(1:length(missingMets),i) = missingMets;
end

% output the results
% these results will be used to generate Dataset_S3 and Dataset_S4
% these results will be used in determining functional explanations for
% essentiality
%xlswrite('model_biomassPrecursorCheck_PA14.xlsx',essGenes_PA14,'LB_results','B1');
%xlswrite('model_biomassPrecursorCheck_PA14.xlsx',essMissingMets_PA14,'LB_results','B2');

%xlswrite('model_biomassPrecursorCheck_PAO1.xlsx',essGenes_PAO1,'LB_results','B1');
%xlswrite('model_biomassPrecursorCheck_PAO1.xlsx',essMissingMets_PAO1,'LB_results','B2');

%% Sputum essentiality

PA14 = changeMedia_SEED(PA14, 2,'');
PA14_sol = optimizeCbModel(PA14)

PAO1 = changeMedia_SEED(PAO1, 2,'');
PAO1_sol = optimizeCbModel(PAO1)

essIDX_PA14 = [];
essGenes_PA14 = [];

for i = 1:length(PA14.genes)
    PA14_deleted = deleteModelGenes(PA14,PA14.genes(i));
    sol = optimizeCbModel(PA14_deleted);
    if sol.f < 0.0001
        essIDX_PA14 = [essIDX_PA14,i];
        essGenes_PA14 = [essGenes_PA14,PA14.genes(i)];
    end
end

essMissingMets_PA14 = cell(59,length(essGenes_PA14)); %59 is the number of metabolites in biomass reaction?

for i = 1:length(essGenes_PA14)
    PA14_deleted = deleteModelGenes(PA14,essGenes_PA14(i));
    [missingMets,presentMets] = biomassPrecursorCheck(PA14_deleted);
    essMissingMets_PA14(1:length(missingMets),i) = missingMets;
end

essIDX_PAO1 = [];
essGenes_PAO1 = [];

for i = 1:length(PAO1.genes)
    PAO1_deleted = deleteModelGenes(PAO1,PAO1.genes(i));
    sol = optimizeCbModel(PAO1_deleted);
    if sol.f < 0.0001
        essIDX_PAO1 = [essIDX_PAO1,i];
        essGenes_PAO1 = [essGenes_PAO1,PAO1.genes(i)];
    end
end

essMissingMets_PAO1 = cell(59,length(essGenes_PAO1)); %59 is the number of metabolites in biomass reaction?

for i = 1:length(essGenes_PAO1)
    PAO1_deleted = deleteModelGenes(PAO1,essGenes_PAO1(i));
    [missingMets,presentMets] = biomassPrecursorCheck(PAO1_deleted);
    essMissingMets_PAO1(1:length(missingMets),i) = missingMets;
end

% output the results
% these results will be used to generate Dataset_S3 and Dataset_S4
% these results will be used in determining functional explanations for
% essentiality
%xlswrite('model_biomassPrecursorCheck_PA14.xlsx',essGenes_PA14,'SCFM_results','B1');
%xlswrite('model_biomassPrecursorCheck_PA14.xlsx',essMissingMets_PA14,'SCFM_results','B2');

%xlswrite('model_biomassPrecursorCheck_PAO1.xlsx',essGenes_PAO1,'SCFM_results','B1');
%xlswrite('model_biomassPrecursorCheck_PAO1.xlsx',essMissingMets_PAO1,'SCFM_results','B2');

%% Pyruvate essentiality

PAO1 = changeMedia_SEED(PAO1, 3,'EX_cpd00020(e)');
PAO1_sol = optimizeCbModel(PAO1)

essIDX_PAO1 = [];
essGenes_PAO1 = [];

for i = 1:length(PAO1.genes)
    PAO1_deleted = deleteModelGenes(PAO1,PAO1.genes(i));
    sol = optimizeCbModel(PAO1_deleted);
    if sol.f < 0.0001
        essIDX_PAO1 = [essIDX_PAO1,i];
        essGenes_PAO1 = [essGenes_PAO1,PAO1.genes(i)];
    end
end

essMissingMets_PAO1 = cell(59,length(essGenes_PAO1)); %59 is the number of metabolites in biomass reaction?

for i = 1:length(essGenes_PAO1)
    PAO1_deleted = deleteModelGenes(PAO1,essGenes_PAO1(i));
    [missingMets,presentMets] = biomassPrecursorCheck(PAO1_deleted);
    essMissingMets_PAO1(1:length(missingMets),i) = missingMets;
end

% output the results
% these results will be used to generate Dataset_S3
% these results will be used in determining functional explanations for
% essentiality
%xlswrite('model_biomassPrecursorCheck_PAO1.xlsx',essGenes_PAO1,'Pyruvate_results','B1');
%xlswrite('model_biomassPrecursorCheck_PAO1.xlsx',essMissingMets_PAO1,'Pyruvate_results','B2');

%% Succinate essentiality

PAO1 = changeMedia_SEED(PAO1, 3,'EX_cpd00036(e)');
PAO1_sol = optimizeCbModel(PAO1)

essIDX_PAO1 = [];
essGenes_PAO1 = [];

for i = 1:length(PAO1.genes)
    PAO1_deleted = deleteModelGenes(PAO1,PAO1.genes(i));
    sol = optimizeCbModel(PAO1_deleted);
    if sol.f < 0.0001
        essIDX_PAO1 = [essIDX_PAO1,i];
        essGenes_PAO1 = [essGenes_PAO1,PAO1.genes(i)];
    end
end

essMissingMets_PAO1 = cell(59,length(essGenes_PAO1)); %59 is the number of metabolites in biomass reaction?

for i = 1:length(essGenes_PAO1)
    PAO1_deleted = deleteModelGenes(PAO1,essGenes_PAO1(i));
    [missingMets,presentMets] = biomassPrecursorCheck(PAO1_deleted);
    essMissingMets_PAO1(1:length(missingMets),i) = missingMets;
end

% output the results
% these results will be used to generate Dataset_S3
% these results will be used in determining functional explanations for
% essentiality
%xlswrite('model_biomassPrecursorCheck_PAO1.xlsx',essGenes_PAO1,'Succinate_results','B1');
%xlswrite('model_biomassPrecursorCheck_PAO1.xlsx',essMissingMets_PAO1,'Succinate_results','B2');