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

%% set exchanges

%Finds the indices of the exchange reactions in the model
exchangerxns = [];
for rxn = 1:length(PA14.rxns);
    exchangerxns(rxn) = strncmp('EX_',PA14.rxns{rxn},3);
end
exchangeindices = find(exchangerxns); 

PA14_LB = PA14;

PA14_LB.lb(exchangeindices) = zeros(size(exchangeindices));
PA14_LB.ub(exchangeindices) = 1000*ones(size(exchangeindices));
PA14_LB.lb(findRxnIDs(PA14_LB,{'EX_cpd00007(e)'})) = -20;  %-1000; %Limit aerobic growth to 20 mmol/gDW/hr of O2 (Added by Phil, EX_O2 removed from openexchanges)

openexchanges = {
'EX_cpd00001(e)'    %H2O
'EX_cpd00009(e)'    %Phosphate
'EX_cpd00011(e)'    %CO2
'EX_cpd00021(e)'    %Fe2+
'EX_cpd00034(e)'    %Zn2+
'EX_cpd00048(e)'    %Sulfate
'EX_cpd00058(e)'    %Cu2+
'EX_cpd00205(e)'    %K+
'EX_cpd00254(e)'    %Mg
'EX_cpd00971(e)'    %Na+
'EX_cpd01012(e)'    %Cd2+
'EX_cpd00067(e)'    %H+
'EX_cpd00528(e)' % added SEED Nitrogen (changed by JAB, 01/14/14, not in older LB media formulations)
};

%Limited nutrients in LB
LBexchanges = {
'EX_cpd00023(e)'    %L-Glutamate
'EX_cpd00027(e)'    %D-Glucose
'EX_cpd00033(e)'    %Glycine
'EX_cpd00035(e)'    %L-Alanine
'EX_cpd00039(e)'    %L-Lysine
'EX_cpd00041(e)'    %L-Aspartate
'EX_cpd00051(e)'    %L-Arginine
'EX_cpd00054(e)'    %L-Serine
'EX_cpd00060(e)'    %L-Methionine
'EX_cpd00065(e)'    %Tryptophan
'EX_cpd00066(e)'    %L-Phenylalanine
'EX_cpd00069(e)'    %L-Tyrosine
'EX_cpd00084(e)'    %L-Cysteine
'EX_cpd00107(e)'    %L-Leucine
'EX_cpd00119(e)'    %L-Histidine
'EX_cpd00129(e)'    %L-Proline
'EX_cpd00156(e)'    %L-Valine
'EX_cpd00161(e)'    %L-Threonine
'EX_cpd00305(e)'    %Thiamin
'EX_cpd00322(e)'    %L-Isoleucine
'EX_cpd00092(e)'    %Uracil
'EX_cpd00307(e)'    %Cytosine
'EX_cpd03091(e)'    %5'-Deoxyadenosine
};

%% initiate loop

essGenes = zeros(length(PA14_LB.genes), 100);
essGenes_length = zeros(length(LBexchanges),100);
overlapList = zeros(1,length(LBexchanges));
overlapList_sample = zeros(length(LBexchanges),100);
aveEssGenesLength = zeros(1,length(LBexchanges));
stdEssGenesLength = zeros(1,length(LBexchanges));
medianEssGenesLength = zeros(1,length(LBexchanges));
maxEssGenesLength = zeros(1,length(LBexchanges));
minEssGenesLength = zeros(1,length(LBexchanges));

for i = 2:(length(LBexchanges)-1)
    clearvars overlap essGenes

    for j = 1:100
        modelsample = PA14_LB;
        
        randNum = datasample(1:length(LBexchanges),i,'Replace',false);

        LBexchanges_sample = LBexchanges(randNum);
        modelsample.lb(find(ismember(PA14.rxns,openexchanges))) = -1000*ones(size(openexchanges));
        modelsample.lb(find(ismember(PA14.rxns,LBexchanges_sample))) = -10*ones(size(LBexchanges_sample));
        modelsample.ub(find(ismember(PA14.rxns,LBexchanges_sample))) = 10*ones(size(LBexchanges_sample));
        
        % there are a lot of instances of no growth when there are few
        % available limited exchanges
        sol = optimizeCbModel(modelsample);
        
        while sol.f < 0.0001
            modelsample = PA14_LB;
        
        randNum = datasample(1:length(LBexchanges),i,'Replace',false);

        LBexchanges_sample = LBexchanges(randNum);
        modelsample.lb(find(ismember(PA14.rxns,openexchanges))) = -1000*ones(size(openexchanges));
        modelsample.lb(find(ismember(PA14.rxns,LBexchanges_sample))) = -10*ones(size(LBexchanges_sample));
        modelsample.ub(find(ismember(PA14.rxns,LBexchanges_sample))) = 10*ones(size(LBexchanges_sample));
        
        % there are a lot of instances of no growth when there are few
        % available limited exchanges
        sol = optimizeCbModel(modelsample);
        
        end
            
        % perform gene essentiality predictions
        if sol.f > 0.0001
            for k = 1:length(modelsample.genes)
                modelsample_deleted = deleteModelGenes(modelsample,modelsample.genes(k));
                sol = optimizeCbModel(modelsample_deleted);
                if sol.f < 0.0001
                    essGenes(k,j) = k;
                end
            end
        else
            essGenes(:,j) = 0;
        end
    end
    
    % find overlap in essential genes
    overlap = intersect(find(essGenes(:,1)),find(essGenes(:,2)));
    overlapList_sample(i,1) = length(overlap);
    for m = 3:100
        overlap = intersect(overlap, find(essGenes(:,m)));
        overlapList_sample(i,m-2) = length(overlap);
    end
    
    overlapList(i) = length(overlap);
    
    % find length of essential gene lists
    for n = 1:100
        essGenes_length(i,n) = length(find(essGenes(:,n)));
    end
    
    aveEssGenesLength(i) = mean(essGenes_length(i,:));
    stdEssGenesLength(i) = std(essGenes_length(i,:));
    medianEssGenesLength(i) = median(essGenes_length(i,:));
    maxEssGenesLength(i) = max(essGenes_length(i,:));
    minEssGenesLength(i) = min(essGenes_length(i,:));
    
    count = i

end

%% output 

xlswrite('model_LBMediaAnalysis_overlapList_sample_180119run.xlsx',overlapList_sample');
