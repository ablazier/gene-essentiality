%% set-up the environment
clear all
close all
clc
initCobraToolbox()
changeCobraSolver('gurobi5');

%% load models

PAO1=load('pao1_updated.mat');
PAO1 = PAO1.PAO1;
PAO1 = changeMedia_SEED(PAO1, 1,'');
PAO1 = changeObjective(PAO1,'PAO1_Biomass',1);
PAO1_sol = optimizeCbModel(PAO1) % confirm value is ~15.7299

PA14=load('pa14_updated.mat');
PA14 = PA14.PA14;
PA14 = changeMedia_SEED(PA14, 1,'');
PA14 = changeObjective(PA14,'PA14_Biomass',1);
PA14_sol = optimizeCbModel(PA14) % confirm value is ~15.7298

%% LB essentiality

PAO1 = changeMedia_SEED(PAO1, 1,'');
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

PA14 = changeMedia_SEED(PA14, 1,'');
PA14_sol = optimizeCbModel(PA14)

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

% output the results
% these results will be used to generate Dataset_SXXX
%xlswrite('model_essentiality_PAO1_updated.xlsx',essGenes_PAO1,'LB_results','B1');
%xlswrite('model_essentiality_PA14_updated.xlsx',essGenes_PA14,'LB_results','B1');

%% Sputum essentiality

PAO1 = changeMedia_SEED(PAO1, 2,'');
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

% output the results
% these results will be used to generate Dataset_SXXX
%xlswrite('model_essentiality_PAO1_updated.xlsx',essGenes_PAO1,'Sputum_results','B1');