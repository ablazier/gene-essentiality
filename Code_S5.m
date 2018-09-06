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

%% Perform minimal media essentiality screens

carbs = {'cpd01949';'cpd00361';'cpd00162';'cpd00294';'cpd01242';'cpd00599';'cpd00361';'cpd00136';'cpd00029';'cpd00142';'cpd00182';'cpd00211';'cpd00331';'cpd01502';'cpd00137';'cpd00266';'cpd00080';'cpd00117';'cpd00158';'cpd00082';'cpd00072';'cpd00280';'cpd00222';'cpd00089';'cpd00079';'cpd00164';'cpd00386';'cpd00314';'cpd00138';'cpd00105';'cpd00609';'cpd00550';'cpd00588';'cpd00666';'cpd00154';'cpd00047';'cpd00106';'cpd00100';'cpd00033';'cpd00155';'cpd00139';'cpd11589';'cpd11592';'cpd11588';'cpd00040';'cpd00851';'cpd00246';'cpd00380';'cpd00035';'cpd11585';'cpd00224';'cpd00051';'cpd00132';'cpd00041';'cpd00023';'cpd00053';'cpd00119';'cpd00227';'cpd00322';'cpd00159';'cpd00107';'cpd00039';'cpd00130';'cpd00060';'cpd00064';'cpd00066';'cpd00129';'cpd00054';'cpd00666';'cpd00161';'cpd00156';'cpd00308';'cpd00179';'cpd01262';'cpd03320';'cpd00121';'cpd00666';'cpd00652';'cpd00477';'cpd00489';'cpd00141';'cpd00118';'cpd00020';'cpd00036';'cpd00184';'cpd00794';'cpd00249';'cpd00581';'cpd00027';'cpd00094';'cpd00024';'cpd00797';'cpd00281'};
carbs_names = {'2,3-Butanediol';'2,3-Butanone';'2-Aminoethanol';'2-Deoxy Adenosine';'2-Deoxy-D-Ribose';'2-Hydroxy Benzoic Acid';'3-Hydroxy 2-Butanone';'4-Hydroxy Benzoic Acid';'Acetic Acid';'Acetoacetic Acid';'Adenosine';'Butyric Acid';'cis-Aconitic Acid';'Citraconic Acid';'Citric Acid';'D,L-Carnitine';['D,L-?-Glycerol- ' char(10) 'Phosphate'];'D-Alanine';'D-Cellobiose';'D-Fructose';'D-Fructose-6-Phosphate';'D-Galacturonic Acid';'D-Gluconic Acid';'D-Glucose-1-Phosphate';'D-Glucose-6-Phosphate';'D-Glucuronic Acid';'D-Malic Acid';'D-Mannitol';'D-Mannose';'D-Ribose';'D-Saccharic Acid';'D-Serine';'D-Sorbitol';'D-Tartaric Acid';'D-Xylose';'Formic Acid';'Fumaric Acid';'Glycerol';'Glycine';'Glycogen';'Glycolic Acid';'Glycyl-L-Aspartic Acid';'Glycyl-L-Glutamic Acid';'Glycyl-L-Proline';'Glyoxylic Acid';'Hydroxy-L-Proline';'Inosine';'Itaconic Acid';'L-Alanine';'L-Alanyl-Glycine';'L-Arabinose';'L-Arginine';'L-Asparagine';'L-Aspartic Acid';'L-Glutamic Acid';'L-Glutamine';'L-Histidine';'L-Homoserine';'L-Isoleucine';'L-Lactic Acid';'L-Leucine';'L-Lysine';'L-Malic Acid';'L-Methionine';'L-Ornithine';'L-Phenylalanine';'L-Proline';'L-Serine';'L-Tartaric Acid';'L-Threonine';'L-Valine';'Malonic Acid';'Maltose';'Maltotriose';'m-Hydroxy Phenyl Acetic Acid';'m-Inositol';'m-Tartaric Acid';'Mucic Acid';'N-Acetyl-L-Glutamic Acid';'p-Hydroxy Phenyl Acetic Acid';'Propionic Acid';'Putrescine';'Pyruvic Acid';'Succinic Acid';'Thymidine';'D-Trehalose';'Uridine';'Urocanic Acid';'?-D-Glucose';'?-Keto-Butyric Acid';'?-Keto-Glutaric Acid';'?-Hydroxy Butyric Acid';'?-Amino Butyric Acid'};

% perform gene essentiality analysis on different minimal media
essIDX_min = zeros(length(PA14.genes), length(carbs));

for i=1:length(carbs)
    [PA14_min,rxn_min]=addExchangeRxn_JB(PA14,strcat(carbs(i),'[e]'),-10,1000); %check if transporter is present
    PA14_min = changeMedia_SEED(PA14_min, 3,rxn_min);
    PA14_min.c=zeros(size(PA14_min.c));
    PA14_min=changeObjective(PA14_min,'PA14_Biomass',1);
    sol = optimizeCbModel(PA14_min);
    
    if sol.f > 0.0001
        for j = 1:length(PA14_min.genes)
            PA14_min_deleted = deleteModelGenes(PA14_min,PA14_min.genes(j));
            sol = optimizeCbModel(PA14_min_deleted);
            if sol.f < 0.0001
                essIDX_min(j,i) = j;
            end
        end
    end
end

% determine on which carbon sources the model grew and essential genes
% were predicted also obtain essential gene lists for each carbon source on
% which the model grew
carbs_grow = [];

for i = 1:length(carbs)
    carbs_essIDX = find(essIDX_min(:,i));
    if isempty(carbs_essIDX) == 0
        carbs_grow = [carbs_grow,carbs(i)];
    end
end

carbs_essGenes = repmat({''},length(PA14.genes),length(carbs_grow));
carbs_count = 0;

for i = 1:length(carbs)
    carbs_essIDX = find(essIDX_min(:,i));
    if isempty(carbs_essIDX) == 0
        carbs_count = carbs_count+1;
        carbs_essGenes(1:length(carbs_essIDX),carbs_count) = PA14.genes(carbs_essIDX);
    end
end

%% Perform sampling of essential gene lists

overlapList = zeros(length(carbs_grow)-1,500);
aveOverlap = zeros(1,length(carbs_grow));
stdOverlap = zeros(1,length(carbs_grow));
medianOverlap = zeros(1,length(carbs_grow));
maxOverlap = zeros(1,length(carbs_grow));
minOverlap = zeros(1,length(carbs_grow));

% run the for-loop for pairwise duplicates all the way through sets of 41
% carbon sources
for i = 2:(length(carbs_grow)-1)
        for j = 1:500
            randNum = randi(length(carbs_grow),1,i);
            randEss = carbs_essGenes(:, randNum);
            overlap = intersect(randEss(:,1),randEss(:,2));
            if i > 2
                for k = 3:i
                    overlap = intersect(overlap,randEss(:,k));
                end
            end
            overlapList(i,j) = length(overlap);
        end
        
        aveOverlap(i) = mean(overlapList(i,:));
        stdOverlap(i) = std(overlapList(i,:));
        medianOverlap(i) = median(overlapList(i,:));
        maxOverlap(i) = max(overlapList(i,:));
        minOverlap(i) = min(overlapList(i,:));
end

%% calculate core and accessory essentials

% determine union of essentials and core essentials
allEss = union(carbs_essGenes(:,1),carbs_essGenes(:,2));
coreEss = intersect(carbs_essGenes(:,1),carbs_essGenes(:,2));
for i = 3:length(carbs_grow)
    allEss = union(allEss,carbs_essGenes(:,i));
    coreEss = intersect(coreEss,carbs_essGenes(:,i));
end

% calculate core essentials
coreNum = length(coreEss)/length(allEss)*100

% calculate accessory essentials
accessoryNum = (length(allEss)-length(coreEss))/length(allEss)*100

%% output data

xlswrite('model_minimalMediaAnalysis_stats.xlsx',aveOverlap');
xlswrite('model_minimalMediaAnalysis_stats.xlsx',stdOverlap',1,'B1');
xlswrite('model_minimalMediaAnalysis_stats.xlsx',medianOverlap',1,'C1');
xlswrite('model_minimalMediaAnalysis_stats.xlsx',maxOverlap',1,'D1');
xlswrite('model_minimalMediaAnalysis_stats.xlsx',minOverlap',1,'E1');

xlswrite('model_minimalMediaAnalysis_overlap.xlsx',overlapList');

