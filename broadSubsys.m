%Jennifer Bartell 12/2012
%categorizing specific subsystem assignments into broader (KEGG) subsystems

function [modnew,subsys,matches,missing,multiple]=broadSubsys(model,syssource)

%clear all; close all; clc
%load('paoLB.mat');

%model=paoLB;

if syssource==1
    sys=model.subSystems;
elseif syssource==2
    sys=model.genepathways(:,3);
elseif syssource==3
    sys=model.KEGGSubsys;
end

% sys={'Glycolysis / Gluconeogenesis'
% 'Citrate cycle(TCA cycle)'
% 'Pentose phosphate pathway'
% 'Inositol metabolism'
% 'Pentose and glucuronate interconversions'
% 'Fructose and mannose metabolism'
% 'Galactose metabolism'
% 'Ascorbate and aldarate metabolism'
% 'Fatty acid biosynthesis'
% 'Fatty acid elongation in mitochondria'
% 'Fatty acid metabolism'
% 'Synthesis and degradation of ketone bodies'
% 'Biosynthesis of steroids'
% 'Ubiquinone biosynthesis'
% 'Urea cycle and metabolism of amino groups'
% 'Purine metabolism'
% 'Pyrimidine metabolism'
% 'Glutamate metabolism'
% 'Alanine and aspartate metabolism'
% 'Tetracycline biosynthesis'
% 'Glycine, serine and threonine metabolism'
% 'Methionine metabolism'
% 'Cysteine metabolism'
% 'Valine,leucine and isoleucine degradation'
% 'Valine,leucine and isoleucine biosynthesis'
% 'Lysine biosynthesis'
% 'Lysine degradation'
% 'Arginine and proline metabolism'
% 'Histidine metabolism'
% 'Tyrosine metabolism'
% 'Phenylalanine metabolism'
% 'gamma-Hexachlorocyclohexane degradation'
% 'Benzoate degradation via hydroxylation'
% 'Bisphenol A degradation'
% 'Tryptophan metabolism'
% 'Phenylalanine tyrosine and tryptophan biosynthesis'
% 'Novobiocin biosynthesis'
% 'Benzoxazinone biosynthesis'
% 'beta-Alanine metabolism'
% 'Taurine and hypotaurine metabolism'
% 'Aminophosphonate metabolism'
% 'Selenoamino acid metabolism'
% 'Cyanoamino acid metabolism'
% 'D-Glutamine and D-glutamate metabolism'
% 'D-Alanine metabolism'
% 'Glutathione metabolism'
% 'Starch and sucrose metabolism'
% 'Nucleotide sugars metabolism'
% 'Streptomycin biosynthesis'
% 'Polyketide sugar unit biosynthesis'
% 'Aminosugars metabolism'
% 'Lipopolysaccharide biosynthesis'
% 'Peptidoglycan biosynthesis'
% 'Glycerolipid metabolism'
% 'Inositol phosphate metabolism'
% 'Glycerophospholipid metabolism'
% 'Pyruvate metabolism'
% 'Biphenyl degradation'
% 'Toluene and xylene degradation'
% '"24-Dichlorobenzoate degradation"'
% '1- and 2-Methylnaphthalene degradation'
% 'Naphthalene and anthracene degradation'
% '"14-Dichlorobenzene degradation"'
% 'Fluorene degradation'
% 'Carbazole degradation'
% 'Glyoxylate and dicarboxylate metabolism'
% 'Benzoate degradation via CoA ligation'
% 'Propanoate metabolism'
% 'Ethylbenzene degradation'
% 'Styrene degradation'
% 'Butanoate metabolism'
% 'One carbon pool by folate'
% 'Methane metabolism'
% 'Carbon fixation in photosynthetic organisms'
% 'Reductive carboxylate cycle(CO2 fixation)'
% 'Thiamine metabolism'
% 'Riboflavin metabolism'
% 'Vitamin B6 metabolism'
% 'Nicotinate and nicotinamide metabolism'
% 'Pantothenate and CoA biosynthesis'
% 'Biotin metabolism'
% 'Folate biosynthesis'
% 'Atrazine degradation'
% 'Porphyrin and chlorophyll metabolism'
% 'Terpenoid biosynthesis'
% 'Nitrogen metabolism'
% 'Sulfur metabolism'
% 'Caprolactam degradation'
% 'Alkaloid biosynthesis I'
% 'Alkaloid biosynthesis II'
% 'Aminoacyl-tRNA biosynthesis'
% 'Biosynthesis of siderophore group nonribosomal peptides'
% 'Biosynthesis of vancomycin group antibiotics'
% 'Phosphatidylinositol signaling system'};

%systosort=unique(paoLB.subSystems(~ismember(paoLB.subSystems,sys)))

%family=subsysKEGGgroup(syscheck);

    for chk=1:length(sys)

        %residual from ToBiN model
        %if ismember(model.rxns(chk),'IR10359')
        %    sys(chk)={'Energy metabolism'};
        %end
        
        if ismember(sys(chk),'Biosynthesis of terpenoids and steroids')
            sys(chk)={'Metabolism of Terpenoids and Polyketides'};
        end
        %organization: 1st name in each row is broad subsystem category, following list
        %is all minor subsystem names filed under that category
        c(1).s={'Carbohydrate Metabolism','Inositol metabolism','Glycolysis / Gluconeogenesis','Citrate cycle (TCA cycle)','Pentose phosphate pathway','Pentose and glucuronate interconversions','Fructose and mannose metabolism','Galactose metabolism','Ascorbate and aldarate metabolism','Starch and sucrose metabolism','Amino sugar and nucleotide sugar metabolism','Pyruvate metabolism','Glyoxylate and dicarboxylate metabolism','Propanoate metabolism','Butanoate metabolism','C5-Branched dibasic acid metabolism','Inositol phosphate metabolism','Insulin signaling pathway',};
        c(2).s={'Energy Metabolism','Oxidative phosphorylation','Photosynthesis','Photosynthesis - antenna proteins','Carbon fixation in photosynthetic organisms','Carbon fixation pathways in prokaryotes','Carbon fixation pathways in prokaryotes','Methane metabolism','Nitrogen metabolism','Sulfur metabolism','Reductive carboxylate cycle (CO2 fixation)'};
        c(3).s={'Lipid Metabolism','Biosynthesis of steroids','Fatty acid biosynthesis','Fatty acid elongation in mitochondria','Fatty acid metabolism','Synthesis and degradation of ketone bodies','Steroid biosynthesis','Primary bile acid biosynthesis','Secondary bile acid biosynthesis','Steroid hormone biosynthesis','Glycerolipid synthesis','Glycerolipid metabolism','Glycerophospholipid metabolism','Ether lipid metabolism','Sphingolipid metabolism','Arachidonic acid metabolism','Linoleic acid metabolism','alpha-Linolenic acid metabolism','Biosynthesis of unsaturated fatty acids'};
        c(4).s={'Nucleotide Metabolism','Purine metabolism','Pyrimidine metabolism'};
        c(5).s={'Amino Acid Metabolism','Glutamate metabolism','Alanine, aspartate and glutamate metabolism','Alanine, aspartate, and glutamate metabolism','Alanine and aspartate metabolism','Glycine, serine and threonine metabolism','Glycine, serine, and threonine metabolism','Methionine metabolism','Cysteine metabolism','Cysteine and methionine metabolism','Valine, leucine and isoleucine degradation','Valine, leucine, and isoleucine degradation','Valine, leucine and isoleucine biosynthesis','Valine, leucine, and isoleucine biosynthesis','Lysine biosynthesis','Lysine degradation','Arginine and proline metabolism','Urea cycle and metabolism of amino groups','Histidine degradation','Histidine metabolism','Tyrosine metabolism','Pheynylalanine metabolism','Phenylalanine metabolism','Tryptophan metabolism','Phenylalanine, tyrosine and tryptophan biosynthesis'};
        c(6).s={'Metabolism of Other Amino Acids','beta-Alanine metabolism','Taurine and hypotaurine metabolism','Aminophosphonate metabolism','Phosphonate and phosphinate metabolism','Selenocompound metabolism','Cyanoamino acid metabolism','D-Glutamine and D-glutamate metabolism','D-Arginine and D-ornithine metabolism','','D-Alanine metabolism','Glutathione metabolism','Selenoamino acid metabolism'};
        c(7).s={'Glycan Biosynthesis and Metabolism','N-Glycan biosynthesis','Various types of N-glycan biosynthesis','Mucin type O-Glycan biosynthesis','Mucin type O-Glycan biosynthesis','Other types of O-glycan biosynthesis','Glycosaminoglycan biosynthesis - chondroitin sulfate','Glycosaminoglycan biosynthesis - heparan sulfate','Glycosaminoglycan biosynthesis - keratan sulfate','Glycosaminoglycan degradation','Glycosylphosphatidylinositol(GPI)-anchor biosynthesis','Glycosphingolipid biosynthesis - lacto and neolacto series','Glycosphingolipid biosynthesis - globo series','Glycosphingolipid biosynthesis - ganglio series','Lipopolysaccharide biosynthesis','Lipopolysaccharide biosynthesis (PAO)','Lipopolysaccharide biosynthesis (PS)','Peptidoglycan biosynthesis','Other glycan degradation'};
        c(8).s={'Metabolism of Cofactors and Vitamins','Folate metabolism','Thiamine metabolism','Riboflavin metabolism','Vitamin B6 metabolism','Nicotinate and nicotinamide metabolism','Pantothenate and CoA biosynthesis','Biotin metabolism','Lipoic acid metabolism','Folate biosynthesis','One carbon pool by folate','Retinol metabolism','Porphyrin and chlorophyll metabolism','Ubiquinone biosynthesis','Ubiquinone metabolism','Ubiquinone and other terpenoid-quinone biosynthesis','Ubiquinone and other terpenoid-quinone biosynthesis '};
        c(9).s={'Metabolism of Terpenoids and Polyketides','Terpenoid biosynthesis','Terpenoid backbone biosynthesis','Monoterpenoid biosynthesis','Sesquiterpenoid biosynthesis','Diterpenoid biosynthesis','Carotenoid biosynthesis','Brassinosteroid biosynthesis','Insect hormone biosynthesis','Zeatin biosynthesis','Limonene and pinene degradation','Geraniol degradation','Type I polyketide structures','Biosynthesis of 12-, 14- and 16-membered macrolides','Biosynthesis of ansamycins','Biosynthesis of type II polyketide backbone','Biosynthesis of type II polyketide products','Tetracycline biosynthesis','Polyketide sugar unit biosynthesis','Nonribosomal peptide structures','Biosynthesis of vancomycin group antibiotics';};
        c(10).s={'Biosynthesis of Other Secondary Metabolites','Phenylpropanoid biosynthesis','Stilbenoid, diarylheptanoid and gingerol biosynthesis','Flavonoid biosynthesis','Flavone and flavonol biosynthesis','Anthocyanin biosynthesis','Isoflavonoid biosynthesis','Indole alkaloid biosynthesis','Isoquinoline alkaloid biosynthesis','Tropane, piperidine and pyridine alkaloid biosynthesis','Acridone alkaloid biosynthesis','Caffeine metabolism','Betalain biosynthesis','Glucosinolate biosynthesis','Benzoxazinoid biosynthesis','Penicillin and cephalosporin biosynthesis','beta-Lactam resistance','Streptomycin biosynthesis','Butirosin and neomycin biosynthesis','Clavulanic acid biosynthesis','Puromycin biosynthesis','Novobiocin biosynthesis';};
        c(11).s={'Xenobiotics Biodegradation and Metabolism','Benzoate degradation via CoA ligation','Naphthalene and anthracene degradation','Benzoate degradation','Benzoate degradation via hydroxylation','1,2-Dichloroethane degradation','2,4-Dichlorobenzoate degradation','Aminobenzoate degradation','Fluorobenzoate degradation','Chloroalkane and chloroalkene degradation','gamma-Hexachlorocyclohexane degradation','Chlorocyclohexane and chlorobenzene degradation','Toluene degradation','Toluene and xylene degradation','Xylene degradation','Nitrotoluene degradation','Ethylbenzene degradation','Styrene degradation','Atrazine degradation','Caprolactam degradation','DDT degradation','Bisphenol degradation','Dioxin degradation','Naphthalene degradation','Polycyclic aromatic hydrocarbon degradation','Metabolism of xenobiotics by cytochrome P450','Drug metabolism - cytochrome P450','Drug metabolism - other enzymes';};
        c(12).s={'Transport','ABC transporters'};
        %c(12).s={'Nada'};
        c(13).s={'Virulence Factors','Biosynthesis of siderophore group nonribosomal peptides','Phenazine biosynthesis','Quorum sensing, AHL','Quorum sensing, PQS','Rhamnolipid biosynthesis','Flagellar assembly','Quorum sensing','Rhamnolipid synthesis','Cepacian synthesis','Alginate synthesis'};
        c(14).s={'Exchange'};
        %c(14).s={'Nada'};
        c(15).s={'Demand','demand',};
        tfall=[];
        % -1 is no matches found
        % -2 is multiple matches found

        %sysparts=regexp(sys(chk),'split','|');

        for cats=1:length(c)
            tf=sum(strcmpi(sys(chk),c(cats).s));

            if tf>0; tf1=1; else tf1=0; end

            tfall=[tfall;tf1 tf];

            sysBroad(cats)=c(cats).s(1);
            
        end
        sysBroad(16)={'Other'};
        sysBroad(17)={'Multiple Subsystems'};
        
        
        if sum(tfall(:,1))==1
            indt=find(tfall(:,1)==1);
            matches(chk,:)=[chk,indt];
            subSysBroad(chk)=sysBroad(indt);
        elseif sum(tfall(:,1))<1
            matches(chk,:)=[chk,16];
            subSysBroad(chk)={'Other'};
        elseif sum(tfall(:,1))>1
            matches(chk,:)=[chk,17];
            subSysBroad(chk)={'Multiple Subsystems'};
        end
        
        if isempty(sys{chk})
            subSysBroad(chk)={'Other'};
        end
    end
    
    subsys=subSysBroad';
    modnew=model;
    modnew.subSysBroad=subSysBroad';
    modnew.subSystems=subSysBroad';
    modnew.subSysOld=model.subSystems;
    missing=find(ismember(subSysBroad,'Other'));
    multiple=find(ismember(subSysBroad,'Multiple Subsystems'));
end