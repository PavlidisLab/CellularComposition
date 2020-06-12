% in this code I am examining the brain-specific modules
% represented in TAN , TSN, GTEx, GTEXCRC and how much of it is due
% to the cell type correction. 
% 0. Inputs: TAN, TSN, GTEX BIN net, GTEX CRC BIN NET, SC BIN NETS. SC
% expresssion (so I can do DE analysis between the subtypes)-
% marker genes. 
% PHASE one
% 1. functions enriched in TAN and TSN 
% 2. functions enriched GTEx , GTExCRC
% 3. Give me the brain specific ones from the union of these two, with
% list of genes. 
% PHASE two
% 4. give me the contribution of the marker genes to these modules
% 5. Give me the presentation of any of these functions in the SC
% datasets (a binary result)
% 6. Give me the count of genes, in any of these modules,
% differentially expressed in any cell type (major cell types?)
% 10. Getting the modules which do not have overlap 

% TODO for the new CTC, get this one! 
book= load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'correctedBinNets_logCorrected_jusResiduals.mat'])

% 0. loading stuff  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 

% affy stuff
load('~/data/general/GPL570GemmaMapNEW.mat')
load('~/data/general/tissueExpGenes/brainExpGenes0.8.mat')
load(['~/networks/tissues/brain/binaryNet_FDR5e-' ...
      '5_0.8Expr_Ind0.10.mat'])
tan = binNet;
clear binNet;

FDR = '0012'
FC = '3'
load(sprintf(['~/resultsAndFigures/firstProject/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))
tsn = finalTable(2).wholeNet;

% for GTEx 
load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat') 
load('~/data/GTEx/Brain_Cortex_expGenes.mat')
load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex.mat'])
load('~/resultsAndFigures/secondProject/GTExRegression/correctedBinNets_logCorrected.mat')

% >>>>> marker genes
% affy markers
load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes.mat'])
load(['~/data/cellTypeVarianceFiles/' ...
      'markersFull.mat'])
load(['~/data/cellTypeVarianceFiles/' ...
      'markerSet.mat'])

% GTEx markers
load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex.mat'])

% The GO stuff
load('~/data/general/GOdata_GTExV6_nonComp.mat')
GTExGO = GOdata;
clear GOdata;
load(['~/data/general/GOdata_GPL570_07_nonComp.mat']);
affyGO = GOdata;
clear GOdata

% 1,2. functions enriched GTEx , GTExCRC TAN, TSN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for the GTEx: 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

wholeFMat = GTExGO.matP;
% terms with <=200 and >= 5 genes
sib = sum(wholeFMat);
inTerms = ((sib >=5) + (sib<=200)) == 2;

GTExExpGenes = GTExFiveNets.nets(1).expGenes;

fMat = wholeFMat(GTExExpGenes, inTerms);

inTermNames = GTExGO.GOTerms(inTerms);
inTermsGOID = GTExGO.GOID(inTerms);
rawGeneCount = sib(inTerms);
expGeneCounts = sum(fMat);
linkCounts = zeros(size(rawGeneCount));
potLinkCounts = zeros(size(rawGeneCount));
ps = size(inTermNames);

% to save: total count of raw genes for that term, count of genes in the
% network, count of links, count of potential links, function term,
% function ID, function p, .fdr  
myNet = GTExFiveNets.nets(1).net01;
totalLinkCount = sum(myNet(:))
totalPotLinkCount = size(myNet,1) * (size(myNet, 1) -1) /2
for i = 1:length(inTermNames)
    i
    thisTermGenes = logical(fMat(:, i));
    smallNet = myNet(thisTermGenes, thisTermGenes);
    linkCounts(i) = sum(smallNet(:));
    potLinkCounts(i) = size(smallNet, 1) * (size(smallNet, 1) -1)/2;
    ps(i) = 1 - hygecdf(linkCounts(i) - 1,  totalPotLinkCount, ...
                        totalLinkCount, potLinkCounts(i));
end

gtexBinNetBlood.inTermNames = inTermNames;
gtexBinNetBlood.inTermsGOID = inTermsGOID;
gtexBinNetBlood.rawGeneCount = rawGeneCount;
gtexBinNetBlood.expGeneCounts = expGeneCounts;
gtexBinNetBlood.linkCounts = linkCounts;
gtexBinNetBlood.potLinkCounts = potLinkCounts;
gtexBinNetBlood.ps = ps;

save('~/resultsAndFigures/secondProject/moduleRep_GTExBlood.mat' ...
     , 'gtexBinNetBlood')

load('~/resultsAndFigures/secondProject/moduleRepresentation/moduleRep_GTExBlood.mat')

% select the terms 
wholeFMat = GTExGO.matP;
% terms with <=200 and >= 5 genes
sib = sum(wholeFMat);
inTerms = ((sib >=5) + (sib<=200)) == 2;

GTExExpGenes = GTExFiveNets.nets(2).expGenes;

fMat = wholeFMat(GTExExpGenes, inTerms);

inTermNames = GTExGO.GOTerms(inTerms);
inTermsGOID = GTExGO.GOID(inTerms);
rawGeneCount = sib(inTerms);
expGeneCounts = sum(fMat);
linkCounts = zeros(size(rawGeneCount));
potLinkCounts = zeros(size(rawGeneCount));
ps = zeros(size(inTermNames));

myNet = GTExFiveNets.nets(3).net01;
totalLinkCount = sum(myNet(:))
totalPotLinkCount = size(myNet,1) * (size(myNet, 1) -1) /2
for i = 1:length(inTermNames)
    i
    thisTermGenes = logical(fMat(:, i));
    smallNet = myNet(thisTermGenes, thisTermGenes);
    linkCounts(i) = sum(smallNet(:));
    potLinkCounts(i) = size(smallNet, 1) * (size(smallNet, 1) -1)/2;
    ps(i) = 1 - hygecdf(linkCounts(i) - 1,  totalPotLinkCount, ...
                        totalLinkCount, potLinkCounts(i));
end

gtexBinNetLiver.inTermNames = inTermNames;
gtexBinNetLiver.inTermsGOID = inTermsGOID;
gtexBinNetLiver.rawGeneCount = rawGeneCount;
gtexBinNetLiver.expGeneCounts = expGeneCounts;
gtexBinNetLiver.linkCounts = linkCounts;
gtexBinNetLiver.potLinkCounts = potLinkCounts;
gtexBinNetLiver.ps = ps;

save('~/resultsAndFigures/secondProject/moduleRep_GTExLiver.mat' ...
     , 'gtexBinNetLiver')

load('~/resultsAndFigures/secondProject/moduleRepresentation/moduleRep_GTExLiver.mat')


% >>>>>
% select the terms 
wholeFMat = GTExGO.matP;
% terms with <=200 and >= 5 genes
sib = sum(wholeFMat);
inTerms = ((sib >=5) + (sib<=200)) == 2;

GTExExpGenes = GTExFiveNets.nets(2).expGenes;

fMat = wholeFMat(GTExExpGenes, inTerms);

inTermNames = GTExGO.GOTerms(inTerms);
inTermsGOID = GTExGO.GOID(inTerms);
rawGeneCount = sib(inTerms);
expGeneCounts = sum(fMat);
linkCounts = zeros(size(rawGeneCount));
potLinkCounts = zeros(size(rawGeneCount));
ps = size(inTermNames);

% to save: total count of raw genes for that term, count of genes in the
% network, count of links, count of potential links, function term,
% function ID, function p, .fdr  
myNet = GTExFiveNets.nets(2).net01;
totalLinkCount = sum(myNet(:))
totalPotLinkCount = size(myNet,1) * (size(myNet, 1) -1) /2
for i = 1:length(inTermNames)
    i
    thisTermGenes = logical(fMat(:, i));
    smallNet = myNet(thisTermGenes, thisTermGenes);
    linkCounts(i) = sum(smallNet(:));
    potLinkCounts(i) = size(smallNet, 1) * (size(smallNet, 1) -1)/2;
    ps(i) = 1 - hygecdf(linkCounts(i) - 1,  totalPotLinkCount, ...
                        totalLinkCount, potLinkCounts(i));
end

gtexBinNet.inTermNames = inTermNames;
gtexBinNet.inTermsGOID = inTermsGOID;
gtexBinNet.rawGeneCount = rawGeneCount;
gtexBinNet.expGeneCounts = expGeneCounts;
gtexBinNet.linkCounts = linkCounts;
gtexBinNet.potLinkCounts = potLinkCounts;
gtexBinNet.ps = ps;

save('~/resultsAndFigures/secondProject/moduleRep_GTExBrainCortex.mat' ...
     , 'gtexBinNet')
load('~/resultsAndFigures/secondProject/moduleRepresentation/moduleRep_GTExBrainCortex.mat')

gtexBinNetIns = (gtexBinNet.ps * length(gtexBinNet.ps)) < .1;

inTermNames = GTExGO.GOTerms(inTerms);
inTermsGOID = GTExGO.GOID(inTerms);
rawGeneCount = sib(inTerms);
expGeneCounts = sum(fMat);
linkCounts = zeros(size(rawGeneCount));
potLinkCounts = zeros(size(rawGeneCount));
ps = size(inTermNames);

myNet = ctc.net01; 
totalLinkCount = sum(myNet(:))
totalPotLinkCount = size(myNet,1) * (size(myNet, 1) -1) /2
for i = 1:length(inTermNames)
    i
    thisTermGenes = logical(fMat(:, i));
    smallNet = myNet(thisTermGenes, thisTermGenes);
    linkCounts(i) = sum(smallNet(:));
    potLinkCounts(i) = size(smallNet, 1) * (size(smallNet, 1) -1)/2;
    ps(i) = 1 - hygecdf(linkCounts(i) - 1,  totalPotLinkCount, ...
                        totalLinkCount, potLinkCounts(i));
end

gtexBinNetCTC.inTermNames = inTermNames;
gtexBinNetCTC.inTermsGOID = inTermsGOID;
gtexBinNetCTC.rawGeneCount = rawGeneCount;
gtexBinNetCTC.expGeneCounts = expGeneCounts;
gtexBinNetCTC.linkCounts = linkCounts;
gtexBinNetCTC.potLinkCounts = potLinkCounts;
gtexBinNetCTC.ps = ps;

save('~/resultsAndFigures/secondProject/moduleRep_GTExBrainCortex_CTC_redo.mat' ...
     , 'gtexBinNetCTC')

save('~/resultsAndFigures/secondProject/moduleRep_GTExBrainCortex_CTC_logCorrected_residues.mat' ...
     , 'gtexBinNetCTC')
load('~/resultsAndFigures/secondProject/moduleRepresentation/moduleRep_GTExBrainCortex_CTC_logCorrected.mat')

% >>> priliminary check: yes, some neuronal functions are missing
% as expected >> the question is if those genes are expressed in
% the SC data and the modules are not represented. Also, I could
% examine the enrichment of AffyTAN and TSN for the brain specific  modules
gtexBinNetCTCIns = (gtexBinNetCTC.ps * length(gtexBinNetCTC.ps)) < ...
    .1;

kado = gtexBinNetCTCIns - gtexBinNetIns;
diffTerms = inTermNames(kado == 1);

fileName = 'diffTerms_temp_ctc.csv'
file = ['~/resultsAndFigures/secondProject/' fileName]
fid = fopen(file, 'w')

% printing the Jaccard Sim
for i = 1:length(diffTerms)
    fprintf(fid, ['%s\n'], diffTerms{i})
end
fclose(fid)

% for affy TAN and TSN
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

load(['~/networks/tissues/brain/binaryNet_FDR5e-' ...
      '5_0.8Expr_Ind0.10.mat'])
tan = binNet;

wholeFMat = affyGO.matP;
% terms with <=200 and >= 5 genes
sib = sum(wholeFMat);
inTerms = ((sib >=5) + (sib<=200)) == 2;

load('~/data/general/tissueExpGenes/brainExpGenes0.8.mat')
affyExpGenes = expGenesInd;

fMat = wholeFMat(expGenesInd, inTerms);

inTermNames = affyGO.GOTerms(inTerms);
inTermsGOID = affyGO.GOID(inTerms);
rawGeneCount = sib(inTerms);
expGeneCounts = sum(fMat);
linkCounts = zeros(size(rawGeneCount));
potLinkCounts = zeros(size(rawGeneCount));
ps = size(inTermNames);

% to save: total count of raw genes for that term, count of genes in the
% network, count of links, count of potential links, function term,
% function ID, function p, .fdr  
myNet = tan(expGenesInd, expGenesInd);
totalLinkCount = sum(myNet(:))
totalPotLinkCount = sum(affyExpGenes) * (sum(affyExpGenes) -1) /2
for i = 1:length(inTermNames)
    i
    thisTermGenes = logical(fMat(:, i));
    smallNet = myNet(thisTermGenes, thisTermGenes);
    linkCounts(i) = sum(smallNet(:));
    potLinkCounts(i) = size(smallNet, 1) * (size(smallNet, 1) -1)/2;
    ps(i) = 1 - hygecdf(linkCounts(i) - 1,  totalPotLinkCount, ...
                        totalLinkCount, potLinkCounts(i));
end

tanFen.inTermNames = inTermNames;
tanFen.inTermsGOID = inTermsGOID;
tanFen.rawGeneCount = rawGeneCount;
tanFen.expGeneCounts = expGeneCounts;
tanFen.linkCounts = linkCounts;
tanFen.potLinkCounts = potLinkCounts;
tanFen.ps = ps;

save('~/resultsAndFigures/secondProject/moduleRep_AffyTAN.mat' ...
     , 'tanFen')
load('~/resultsAndFigures/secondProject/moduleRepresentation/moduleRep_AffyTAN.mat')

wholeFMat = affyGO.matP;
% terms with <=200 and >= 5 genes
sib = sum(wholeFMat);
inTerms = ((sib >=5) + (sib<=200)) == 2;

affyExpGenes = expGenesInd;

fMat = wholeFMat(expGenesInd, inTerms);

inTermNames = affyGO.GOTerms(inTerms);
inTermsGOID = affyGO.GOID(inTerms);
rawGeneCount = sib(inTerms);
expGeneCounts = sum(fMat);
linkCounts = zeros(size(rawGeneCount));
potLinkCounts = zeros(size(rawGeneCount));
ps = size(inTermNames);


% to save: total count of raw genes for that term, count of genes in the
% network, count of links, count of potential links, function term,
% function ID, function p, .fdr  
myNet = tsn(expGenesInd, expGenesInd);
totalLinkCount = sum(myNet(:))
totalPotLinkCount = sum(affyExpGenes) * (sum(affyExpGenes) -1) /2
for i = 1:length(inTermNames)
    i
    thisTermGenes = logical(fMat(:, i));
    smallNet = myNet(thisTermGenes, thisTermGenes);
    linkCounts(i) = sum(smallNet(:));
    potLinkCounts(i) = size(smallNet, 1) * (size(smallNet, 1) -1)/2;
    ps(i) = 1 - hygecdf(linkCounts(i) - 1,  totalPotLinkCount, ...
                        totalLinkCount, potLinkCounts(i));
end

tsnsFen.inTermNames = inTermNames;
tsnFen.inTermsGOID = inTermsGOID;
tsnFen.rawGeneCount = rawGeneCount;
tsnFen.expGeneCounts = expGeneCounts;
tsnFen.linkCounts = linkCounts;
tsnFen.potLinkCounts = potLinkCounts;
tsnFen.ps = ps;

save('~/resultsAndFigures/secondProject/moduleRepresentation/moduleRep_AffyTSN.mat' ...
     , 'tsnFen')
load('~/resultsAndFigures/secondProject/moduleRepresentation/moduleRep_AffyTSN.mat')

fdr = ps * length(ps);
printFs = inTermNames(fdr < .1);

fileName = 'enrichedTerms.csv'
file = ['~/resultsAndFigures/secondProject/' fileName]
fid = fopen(file, 'w')

for i = 1:length(printFs)
    fprintf(fid, ['%s\n'], printFs{i})
end
fclose(fid)

% for TAN blood 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['~/networks/tissues/blood/binaryNet_FDR5e-' ...
      '5_0.8Expr_Ind0.10.mat'])
tanBlood = binNet;
load('~/data/general/tissueExpGenes/bloodExpGenes0.8.mat')

wholeFMat = affyGO.matP;
% terms with <=200 and >= 5 genes
sib = sum(wholeFMat);
inTerms = ((sib >=5) + (sib<=200)) == 2;

affyExpGenes = expGenesInd;

fMat = wholeFMat(expGenesInd, inTerms);

inTermNames = affyGO.GOTerms(inTerms);
inTermsGOID = affyGO.GOID(inTerms);
rawGeneCount = sib(inTerms);
expGeneCounts = sum(fMat);
linkCounts = zeros(size(rawGeneCount));
potLinkCounts = zeros(size(rawGeneCount));
ps = size(inTermNames);

% to save: total count of raw genes for that term, count of genes in the
% network, count of links, count of potential links, function term,
% function ID, function p, .fdr  
myNet = tanBlood(expGenesInd, expGenesInd);
totalLinkCount = sum(myNet(:))
totalPotLinkCount = sum(affyExpGenes) * (sum(affyExpGenes) -1) /2
for i = 1:length(inTermNames)
    i
    thisTermGenes = logical(fMat(:, i));
    smallNet = myNet(thisTermGenes, thisTermGenes);
    linkCounts(i) = sum(smallNet(:));
    potLinkCounts(i) = size(smallNet, 1) * (size(smallNet, 1) -1)/2;
    ps(i) = 1 - hygecdf(linkCounts(i) - 1,  totalPotLinkCount, ...
                        totalLinkCount, potLinkCounts(i));
end

tanBloodFen.inTermNames = inTermNames;
tanBloodFen.inTermsGOID = inTermsGOID;
tanBloodFen.rawGeneCount = rawGeneCount;
tanBloodFen.expGeneCounts = expGeneCounts;
tanBloodFen.linkCounts = linkCounts;
tanBloodFen.potLinkCounts = potLinkCounts;
tanBloodFen.ps = ps;

save('~/resultsAndFigures/secondProject/moduleRep_AffyTANBlood.mat' ...
     , 'tanBloodFen')
load(['~/resultsAndFigures/secondProject/' ...
      'moduleRep_AffyTANBlood.mat'])

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
load(['~/networks/tissues/liver/binaryNet_FDR5e-' ...
      '5_0.8Expr_Ind0.10.mat'])
tanLiver = binNet;
load('~/data/general/tissueExpGenes/liverExpGenes0.8.mat')

wholeFMat = affyGO.matP;
% terms with <=200 and >= 5 genes
sib = sum(wholeFMat);
inTerms = ((sib >=5) + (sib<=200)) == 2;

affyExpGenes = expGenesInd;

fMat = wholeFMat(expGenesInd, inTerms);

inTermNames = affyGO.GOTerms(inTerms);
inTermsGOID = affyGO.GOID(inTerms);
rawGeneCount = sib(inTerms);
expGeneCounts = sum(fMat);
linkCounts = zeros(size(rawGeneCount));
potLinkCounts = zeros(size(rawGeneCount));
ps = size(inTermNames);

% to save: total count of raw genes for that term, count of genes in the
% network, count of links, count of potential links, function term,
% function ID, function p, .fdr  
myNet = tanLiver(expGenesInd, expGenesInd);
totalLinkCount = sum(myNet(:))
totalPotLinkCount = sum(affyExpGenes) * (sum(affyExpGenes) -1) /2
for i = 1:length(inTermNames)
    i
    thisTermGenes = logical(fMat(:, i));
    smallNet = myNet(thisTermGenes, thisTermGenes);
    linkCounts(i) = sum(smallNet(:));
    potLinkCounts(i) = size(smallNet, 1) * (size(smallNet, 1) -1)/2;
    ps(i) = 1 - hygecdf(linkCounts(i) - 1,  totalPotLinkCount, ...
                        totalLinkCount, potLinkCounts(i));
end

tanLiverFen.inTermNames = inTermNames;
tanLiverFen.inTermsGOID = inTermsGOID;
tanLiverFen.rawGeneCount = rawGeneCount;
tanLiverFen.expGeneCounts = expGeneCounts;
tanLiverFen.linkCounts = linkCounts;
tanLiverFen.potLinkCounts = potLinkCounts;
tanLiverFen.ps = ps;

save('~/resultsAndFigures/secondProject/moduleRep_AffyTANLiver.mat' ...
     , 'tanLiverFen')
load('~/resultsAndFigures/secondProject/moduleRep_AffyTANLiver.mat')

% putting together the functions: tanFen, tsnFen, gtexBinNet, gtexBinNetCTC
% getting the ID and termNames for each result 

% tan
fdr = .05;

thisPs = tanFen.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < fdr;

Ids = tanFen.inTermsGOID(selected);
terms = tanFen.inTermNames(selected);

% tsn
thisPs = tsnFen.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < fdr;
sum(selected) 

Ids = [Ids tsnFen.inTermsGOID(selected)];
terms = [terms tsnFen.inTermNames(selected)];

% gtexBinNet
thisPs = gtexBinNet.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < fdr;
sum(selected) 

Ids = [Ids gtexBinNet.inTermsGOID(selected)];
terms = [terms gtexBinNet.inTermNames(selected)];

% gtexBinNetCTC
thisPs = gtexBinNetCTC.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < fdr;
sum(selected) 

Ids = [Ids gtexBinNetCTC.inTermsGOID(selected)];
terms = [terms gtexBinNetCTC.inTermNames(selected)];

% gtexBinNetBlood
thisPs = gtexBinNetBlood.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < fdr;
sum(selected) 

Ids = [Ids gtexBinNetBlood.inTermsGOID(selected)];
terms = [terms gtexBinNetBlood.inTermNames(selected)];

% gtexBinNetLiver
thisPs = gtexBinNetLiver.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < fdr;
sum(selected) 

Ids = [Ids gtexBinNetBlood.inTermsGOID(selected)];
terms = [terms gtexBinNetBlood.inTermNames(selected)];

% tanliver
thisPs = tanLiverFen.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < fdr;
sum(selected) 

Ids = [Ids tanLiverFen.inTermsGOID(selected)];
terms = [terms tanLiverFen.inTermNames(selected)];

% tanblood
thisPs = tanBloodFen.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < fdr;
sum(selected) 

Ids = [Ids tanBloodFen.inTermsGOID(selected)];
terms = [terms tanBloodFen.inTermNames(selected)];

% >> marking the terms based on presence in the networks
[uniqueTermIDs, inds] = (unique(Ids));
uniqueTermNames = terms(inds);
funPresence = zeros(8, length(uniqueTermNames));

% tan
thisPs = tanFen.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < fdr;

thisIds = tanFen.inTermsGOID(selected);
[a, b] = ismember(thisIds, uniqueTermIDs);
funPresence(1, b) = 1;

% tsn
thisPs = tsnFen.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < fdr;

thisIds = tsnFen.inTermsGOID(selected);
[a, b] = ismember(thisIds, uniqueTermIDs);
funPresence(2, b) = 1;

% gtexBinNet
thisPs = gtexBinNet.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < fdr;

thisIds = gtexBinNet.inTermsGOID(selected);
[a, b] = ismember(thisIds, uniqueTermIDs);
funPresence(3, b) = 1;

% gtexBinNetCTC
thisPs = gtexBinNetCTC.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < fdr;

thisIds = gtexBinNetCTC.inTermsGOID(selected);
[a, b] = ismember(thisIds, uniqueTermIDs);
funPresence(4, b) = 1;

% gtexBinNetBlood
thisPs = gtexBinNetBlood.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < fdr;

thisIds = gtexBinNetBlood.inTermsGOID(selected);
[a, b] = ismember(thisIds, uniqueTermIDs);
funPresence(5, b) = 1;

% gtexBinNetLiver
thisPs = gtexBinNetLiver.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < fdr;

thisIds = gtexBinNetLiver.inTermsGOID(selected);
[a, b] = ismember(thisIds, uniqueTermIDs);
funPresence(6, b) = 1;

% tan liver
thisPs = tanLiverFen.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < fdr;

thisIds = tanLiverFen.inTermsGOID(selected);
[a, b] = ismember(thisIds, uniqueTermIDs);
funPresence(7, b) = 1;

% tanblood
thisPs = tanBloodFen.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < fdr;

thisIds = tanBloodFen.inTermsGOID(selected);
[a, b] = ismember(thisIds, uniqueTermIDs);
funPresence(8, b) = 1;

funcRes.funPresence = funPresence;
funcRes.uniqueTermIDs = uniqueTermIDs;
funcRes.uniqueTermNames = uniqueTermNames;

save('~/resultsAndFigures/secondProject/moduleRepresentation/funcRes_withAffyBloodLiver_FDR05.mat', ...
     'funcRes')

load('~/resultsAndFigures/secondProject/moduleRepresentation/funcRes_withAffyBloodLiver_FDR05.mat')

load(['~/resultsAndFigures/secondProject/moduleRepresentation/' ...
      'funcRes.mat'])

sumPresence = sum(funcRes.funPresence(1:4, :));
filterFor = sum(funcRes.funPresence(5:6, :)) > 0;
nonSpecificBrain = ((sumPresence > 0) + filterFor) == 2;

brainAll = sumPresence > 0;
brainFiltered = (brainAll - filterFor) > 0;

selectedFunPresence = funcRes.funPresence(1:4, brainFiltered);
selectedTermNames = funcRes.uniqueTermNames(brainFiltered);
selectedTermIDs = funcRes.uniqueTermIDs(brainFiltered);
selectedSumPresence = sumPresence(brainFiltered);

fileName = 'enrichedTerms_TanTsnGTExGTEXCTC_GTExBloodGTExLiver.tsv'
file = ['~/resultsAndFigures/secondProject/' fileName]
fid = fopen(file, 'w')

fprintf(fid, ['ID\tterm\ttan\ttsn\tgtex\tgtexCTC\t', ...
              'sum\n'])
for i = 1:length(brainFiltered)
    myTerm = selectedTermNames{i};
    sib = strtok(myTerm, ',');
    fprintf(fid, ['%d\t%s\t%d\t%d\t%d\t%d\t%d\n'], selectedTermIDs(i), ...
            sib, selectedFunPresence(1, i), selectedFunPresence(2, ...
                                                      i), selectedFunPresence(3, ...
                                                      i), selectedFunPresence(4, ...
                                                      i), selectedSumPresence(i));
end
fclose(fid)

% marking the brain specific terms
% printing the terms, ids and sum
fileName = 'enrichedTerms_allNets.tsv'
file = ['~/resultsAndFigures/secondProject/' fileName]
fid = fopen(file, 'w')

fprintf(fid, ['ID\tterm\ttan\ttsn\tgtex\tgtexCTC\t', ...
              'sum\n'])
for i = 1:length(funPresence)
    myTerm = uniqueTermNames{i};
    sib = strtok(myTerm, ',');
    fprintf(fid, ['%d\t%s\t%d\t%d\t%d\t%d\t%d\n'], uniqueTermIDs(i), ...
            sib, funPresence(1, i), funPresence(2, ...
                                                      i), funPresence(3, ...
                                                      i), funPresence(4, ...
                                                      i), sumPresence(i));
end
fclose(fid)

% TODO - DONE. 1. get the affy GO terms updated (some terms seem to be
% old/new) - 2. get the blood gtex as negative (it is ok if
% negative is not brain specific, since we are going to check the
% expression levels for batches of genes  - no harm.) 

% TODO as a background info, I want the count of genes expressed,
% count of links and the actual pvalues for each of these
% terms. This reference will be used when we are checking the DE
% results 

% PHASE two
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO extracting the sets of genes

% TODO differential expression for the single cell data for group
% of genes with enriched functional terms: 
% for each functional term enriched in the brain networks (or even
% ctc), fetch the genes from GTEx FMAT - Do an ftest to see if
% those genes are expressed in a cell type specific manner in the
% SC data - if the count of DE genes is more than expected by
% random, look it closer. 
[a, b] = ismember(funcRes.uniqueTermIDs, GTExGO.GOID);
fMat = GTExGO.matP(:, b(a));
geneCount = sum(fMat);
hist(geneCount)

myGenes = GTExGO.geneSymbols(logical(fMat(:, 1)));
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% From here, I should extract a combination of bulk tissue brain
% specific functional terms. I should examine (represent) their
% presence in the  two GTEx networks (like a heatmap)

% in the next step, I should examine the expression level of these
% sets of genes in the single cell data  >> are these sets of genes
% expressed differentially in any of the SC data? 

% the following is copy pasted, fix it and things. 
dataFolder = ['~/data/brainSingleCell/']
load('~/data/brainSingleCell/filDataSet_intronAndExon_V4.mat')
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exonAndIntron_V4.mat'])
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exonAndIntron_V4_clusterLabels.mat'])

% the following is copy pasted, fix it and things. 
dataFolder = ['~/data/brainSingleCell/']
load('~/data/brainSingleCell/filDataSet_exon_V4.mat')
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exon_V4.mat'])
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exon_V4_clusterLabels.mat'])


load('~/resultsAndFigures/secondProject/moduleRepresentation/funcRes.mat')

expMat = filDataSet.expMat(:, clusterMeta.inCells);
% normalizing for the million count
sumExp = sum(expMat);
milFac = sumExp ./ 1000000;
normExp = zeros(size(expMat));
sampleCount = size(expMat, 2);
for i = 1:sampleCount
    normExp(:, i) = expMat(:, i)./milFac(i);
end

% get the quantiles for the median network and things 
% I am looking at 31 clusters > those with > 100 samples
inClusterCount = 31;
gCount = size(expMat, 1);
netNames = cell(1, 31);
medMat = zeros(gCount, inClusterCount);
qs = zeros(100, 31);
for i = 1:inClusterCount
    i
    netNames{i} = clusterMeta.sortedClusterNames(i);
    % getting the exp  mat
    [a, b] = ismember(clusterMeta.clusters, ...
                      clusterMeta.sortedClusterNames(i));
    thisExpMat = log2(normExp(:, a)+1);
    thisSCount = sum(a);
    thisExpGenes = logical(sum(thisExpMat > 0,2) > thisSCount/5);
    
    thisExpGenesIDs = find(thisExpGenes);
    for j = 1:length(thisExpGenesIDs)
        myrow = thisExpGenesIDs(j);
        inSamples = thisExpMat(myrow, :) > 0;
        medMat(myrow, i) = median(thisExpMat(myrow, inSamples));
    end
    
    %    medMat(:, i) = median(thisExpMat');
    
    % sib = corr(thisExpMat');
    % upCorr = sib(logical(triu(ones(size(sib)), 1)));
    % qs(:, i) = quantile(upCorr, [0:.01:.99]);
end
% >> medmat with no correction for expression
medMat(:, [11 17 18]) = []; % removing the non neural clusters 

tic
[sib, p2] = corr(medMat');
toc

upP2 = p2(logical(triu(ones(size(p2)),1)));
upSib = sib(logical(triu(ones(size(sib)),1)));
fdr = upP2 * sum(upSib > 0);
fdr = upP2 * size(upP2,1);

temp1 = p2 < (.1/(sum(upSib >0)));
temp2 = sib > 0;
medNet2 = (temp1 + temp2) == 2;
upMedNet2 = triu(medNet2, 1);

tic
[sib, p1] = corr(medMat');
toc

upP1 = p1(logical(triu(ones(size(p1)),1)));
fdr = upP1 * size(upP1,1);

medNet = p1 < (.1/(length(upP1)));
upMedNet = triu(medNet, 1);

% > medmat with correction for expression
tic
[sibc, p1c] = corr(medMat');
toc

upP1c = p1c(logical(triu(ones(size(p1c)),1)));
fdrc = upP1c * size(upP1c,1);

medNetc = p1c < (.1/(length(upP1c)));
upMedNetc = triu(medNetc);

medNets.unCorrected = upMedNet;
medNets.corrected = upMedNetc;
medNets.nonNeuronsRemoved = upMedNet2;
medNets.oneSided = upMedNet2;
medNets.oneSidedNNR = upMedNet2;
medNets.geneSymbols = filDataSet.geneSyms;

heatmap(corr(medMat))

save(['~/resultsAndFigures/secondProject/moduleRepresentation/' ...
      'medianNetworks_31Clusters.mat'], 'medNets')

save(['~/resultsAndFigures/secondProject/moduleRepresentation/' ...
      'SCclusterQuantiles.mat'], 'qs')

load(['~/resultsAndFigures/secondProject/moduleRepresentation/' ...
      'SCclusterQuantiles.mat'])

fCount = length(funcRes.uniqueTermIDs);
sib = corr(medMat(:, 1:31)');
upCorr = sib(logical(triu(ones(size(sib)), 1)));
allMedMatqs = quantile(upCorr, [0:.01:.99]);
medLinkCounts90 = zeros(fCount, 1);
medLinkCounts95 = zeros(fCount, 1);
medLinkCounts99 = zeros(fCount, 1);
scFGeneCount = zeros(1, fCount);
for j = 1:fCount 
    j
    % getting the genes 
    myGenes = GTExGO.geneSymbols(logical(fMat(:, j)));
    [a, b] = ismember(myGenes, filDataSet.geneSyms);
    scFGeneCount(j) = sum(a);
    smallMedMat = medMat(b(a), :);
    
    if(sum(a) >= 3)
        sib = corr(smallMedMat');
        upCorr = sib(logical(triu(ones(size(sib)), 1)));
        medLinkCounts90(j) = sum(upCorr > allMedMatqs(90));
        medLinkCounts95(j) = sum(upCorr > allMedMatqs(95));
        medLinkCounts99(j) = sum(upCorr > allMedMatqs(99));
    end
end

mrRes.medLinkCounts90 = medLinkCounts90;
mrRes.medLinkCounts95 = medLinkCounts95;
mrRes.medLinkCounts99 = medLinkCounts99;

fCount = length(uniqueTermIDs);
[a, b] = ismember(uniqueTermIDs, GTExGO.GOID);
geneCount = sum(fMat);
linkCounts90 = zeros(fCount, 31);
linkCounts95 = zeros(fCount, 31);
linkCounts99 = zeros(fCount, 31);
for i = 1:inClusterCount
    i
    netNames{i} = clusterMeta.sortedClusterNames(i);
    % getting the exp  mat
    [a, b] = ismember(clusterMeta.clusters, ...
                      clusterMeta.sortedClusterNames(i));
    thisExpMat = log2(normExp(:, a)+1);
    
    medMat(:, i) = median(thisExpMat');
    

    for j = 1:fCount 
        % getting the genes 
        j
        myGenes = GTExGO.geneSymbols(logical(fMat(:, j)));
        [a, b] = ismember(myGenes, filDataSet.geneSyms);
        smallExpMat = thisExpMat(b(a), :);
        
        if(sum(a) >= 3)
        sib = corr(smallExpMat');
        upCorr = sib(logical(triu(ones(size(sib)), 1)));
        linkCounts90(j, i) = sum(upCorr > qs(90, i));
        linkCounts95(j, i) = sum(upCorr > qs(95, i));
        linkCounts99(j, i) = sum(upCorr > qs(99, i));
        end
    end
end

mrRes.linkCounts90 = linkCounts90;
mrRes.linkCounts95 = linkCounts95;
mrRes.linkCounts99 = linkCounts99;

save(['~/resultsAndFigures/secondProject/moduleRepresentation/' ...
      'mrRes.mat'], 'mrRes')

load(['~/resultsAndFigures/secondProject/moduleRepresentation/' ...
      'mrRes.mat'])
load(['~/resultsAndFigures/secondProject/moduleRepresentation/' ...
      'funcRes.mat'])
load('~/resultsAndFigures/secondProject/moduleRep_GTExBrainCortex.mat')
%[a,b] = ismember(uniqueTermIDs(nonSpecificBrain) ,
%gtexBinNet.inTermsGOID);

brainFiltered = ((sum(funcRes.funPresence(1:4, :)) > 0) - ...
    sum(funcRes.funPresence(5:6, :))) > 0;

nonSpecificBrain = (((sum(funcRes.funPresence(5:6, :)) > 0) + ...
                    (sum(funcRes.funPresence(1:4, :)) > 0))) == 2;

brainFilteredIDs = find(brainFiltered);
[a,b] = ismember(funcRes.uniqueTermIDs(brainFilteredIDs) , ...
                 gtexBinNet.inTermsGOID);
brainFilteredIDs = brainFilteredIDs(a);

otherLinkCounts_med = mrRes.medLinkCounts99(nonSpecificBrain);
brainLinkCounts_med = mrRes.medLinkCounts99(brainFilteredIDs);

otherLinkCounts_net1 = mrRes.linkCounts99(nonSpecificBrain,9);
brainLinkCounts_net1 = mrRes.linkCounts99(brainFilteredIDs,9);

% for non specific functions, count of links in the scNet1 to the
% count of links in the median thing
ratioOther = otherLinkCounts_net1 ./ (otherLinkCounts_med+1); 
% for general function I don't
% expect to see an imbalance. 
% for specific brain functions, count of links in the net1 to the
% count of links in the median: we expect that count of links in
% the brain sc be lower, due to cell type specific in the median thing
ratioBrain = brainLinkCounts_net1 ./(brainLinkCounts_med+1);

ratioOther = otherLinkCounts_med./(otherLinkCounts_net1+1);
ratioBrain = (brainLinkCounts_med)./ ...
    (brainLinkCounts_net1+1);

% getting the gtex link count
[a,b] = ismember(funcRes.uniqueTermIDs(nonSpecificBrain) , ...
                 gtexBinNet.inTermsGOID);
otherLinkCounts_bulk = gtexBinNet.linkCounts(b(a))';
otherLinkCounts_ctc = gtexBinNetCTC.linkCounts(b(a))';
[a,b] = ismember(funcRes.uniqueTermIDs(brainFilteredIDs) , ...
                 gtexBinNet.inTermsGOID);
brainLinkCounts_bulk = gtexBinNet.linkCounts(b(a))';
brainLinkCounts_ctc = gtexBinNetCTC.linkCounts(b(a))';

ratioOtherBulk = otherLinkCounts_bulk ./ (otherLinkCounts_net1+1);
ratioBrainBulk = brainLinkCounts_bulk ./ (brainLinkCounts_net1+1);

% which term has the most loss?
lossRatioBrain = (brainLinkCounts_bulk - brainLinkCounts_net1) ./ ...
    max(brainLinkCounts_bulk, brainLinkCounts_net1);
myTerms = funcRes.uniqueTermNames(brainFilteredIDs);
[a, b] = sort(lossRatioBrain, 'descend');
sortedTerms = myTerms(b);
sorted_bulk = brainLinkCounts_bulk(b);

halva=[a, sorted_bulk];

% correlation is a sanity check, it is hard to interpret : it is high
kado = [otherLinkCounts_med otherLinkCounts_net1 ...
        otherLinkCounts_bulk otherLinkCounts_ctc];
sib = corr(kado, 'type', 'Spearman')
h = figure
heatmap({'med', 'net1', 'bulk', 'ctc'}, {'med', 'net1', 'bulk', ...
                    'ctc'}, sib);
title('other')
book = [brainLinkCounts_med brainLinkCounts_net1 brainLinkCounts_bulk, ...
       brainLinkCounts_ctc];
sib = corr(book, 'type', 'Spearman')
h = figure
title('brain')
heatmap({'med', 'net1', 'bulk', 'ctc'}, {'med', 'net1', 'bulk', ...
                    'ctc'}, sib);
title('brain')


% however, for specific brain functions, it might mean sth if the
% order is different between bulk (combination of signal),
% med(tissue specific signal), net1 (tissue based signal) and the
% ctc. The question is, do we have "some" brain functions who
% change order dramatically when we use med?  

% >>> sorting the counts and checking which is which!
myTerms = funcRes.uniqueTermNames(brainFilteredIDs);
[a_med, b_med] = sort(brainLinkCounts_med, 'descend');
[a_bulk, b_bulk] = sort(brainLinkCounts_bulk, 'descend');
[a_net1, b_net1] = sort(brainLinkCounts_net1, 'descend');
[a_ctc, b_ctc] = sort(brainLinkCounts_ctc, 'descend');

[ao_med, bo_med] = ismember(b_bulk, b_med);
[ao_net1, bo_net1] = ismember(b_bulk, b_net1);
[ao_ctc, bo_ctc] = ismember(b_bulk, b_ctc);

halva = [bo_bulk, bo_net1, bo_ctc];

% >>>> plotting the lc for the brain functions
myTerms = funcRes.uniqueTermNames(brainFilteredIDs);
[a_med, b_med] = sort(brainLinkCounts_med./potLink', 'descend');
[a_bulk, b_bulk] = sort(brainLinkCounts_bulk./potLink', 'descend');
[a_net1, b_net1] = sort(brainLinkCounts_net1./potLink', 'descend');
[a_ctc, b_ctc] = sort(brainLinkCounts_ctc./potLink', 'descend');

sortedMyTerms = myTerms(b_bulk);
sort_bulk = a_bulk;
sort_net1 = brainLinkCounts_net1(b_bulk)./(potLink(b_bulk)');
sort_med = brainLinkCounts_med(b_bulk)./(potLink(b_bulk)');
sort_ctc =  brainLinkCounts_ctc(b_bulk)./(potLink(b_bulk)');

book = [sort_bulk sort_med sort_net1 sort_ctc];
heatmap((book'))

% >>>> plotting the lc for the other functions
myTerms = funcRes.uniqueTermNames(nonSpecificBrain);
[a_med, b_med] = sort(otherLinkCounts_med ./potLink', 'descend');
[a_bulk, b_bulk] = sort(otherLinkCounts_bulk./potLink', 'descend');
[a_net1, b_net1] = sort(otherLinkCounts_net1./potLink', 'descend');
[a_ctc, b_ctc] = sort(otherLinkCounts_ctc./potLink', 'descend');

sortedMyTerms = myTerms(b_bulk);
sort_bulk = a_bulk;
sort_net1 = otherLinkCounts_net1(b_bulk)./potLink(b_bulk)';
sort_med = otherLinkCounts_med(b_bulk)./potLink(b_bulk)';
sort_ctc =  otherLinkCounts_ctc(b_bulk)./potLink(b_bulk)';

book = [sort_bulk sort_med sort_net1 sort_ctc];
corr(book(1:500, :))
heatmap((book(1:400, :)'))
% GOOD, check the genes for the interesting functions

myTerms = funcRes.uniqueTermNames(nonSpecificBrain);
[a, b] = sort(otherLinkCounts_bulk, 'descend');
sortedMyTerms = myTerms(b);
sort_med = a;
sort_net1 = otherLinkCounts_net1(b);
sort_bulk = otherLinkCounts_med(b);

kado = [sort_med sort_net1 sort_bulk];

% TODO review the preivous list and see where you are. I think we
% can have a result/discussion on the function representation - for
% interesing brain terms (in each of the networks -31) check the
% genes and give a discussion. The fact that there IS a
% correlation, suggests that the terms are not random, but mark
% them with the significance and check  the correlation for
% that. remember, we can only discuss the terms when there IS
% significance > check the gene expression and etc for interesting
% brain terms, where we have drop between the cell
% types. Basically, for every fucntion that is different from bulk
% to CTC, do you see differential expression patterns between
% subtypes? The median should select the links in favor of the cell
% type proportion induced links, exactly the oppositite of the
% direction we go from bulk to CTC. So, median and CTC should
% complement each other on the "meaningful" coexpression, while CTC
% should be more on the direction of net1  

% I need to normalize the function gene count. Just the median is
% correlated wit the observed correlation, but so is SC, which
% shows that we have both cell type effect and bulk tissue effect. 

% compare them to the link count in the bulk tissue, is the med
% effect similar to the effect of bulk tissue? (pleae don't tell me
% it is more...)

% after correction for the cell type, brain specific terms lose
% almost 40% of their representation in the network (are there any
% terms still specific?) and non specific terms gain 15% extra
% links. (compare sum(book) and sum(kado)). I tried the median
% network to capture potential cell type specific terms. However,
% for the med terms, I should first identify if I find any term
% enriched, either for the brain-specific terms or for the general
% terms >> I know that some of the more clear brain specific terms
% are enriched in the median network (is there a clear differential
% expression? they might be non-markers, but have clear DE)
% Therefore, 
% TODO 1: check the significance of the brain terms in median
% netowrk. Are any of them caused by clear differential expression? 
% TODO 2: compare this to the significance of the general terms in
% the median network. IF it works, we can say that enrichment of
% these terms in the bulk tissue are highly affected by the cell
% type proportion: however, if it is cell type proportion, how come
% the difference is not present in the sc and bulk comparison?
% TODO 3: the small difference between SC count of links and bulk
% count of links suggests little effect of cell type proportion in
% the enrichment of functional terms, however, maybe some functions
% are more affected : I did a check of the terms and it is really
% hard to say. I don't think we can say anything here. Maybe there
% is a one function or two, but it is hard to tell. 
% TODO 4: check that cell type specific clusters. 
% REMINDER We are investigating if the representation of the any
% functional modules is affected by the cell type proportion.
% TODO Our brain specific links were highly reproducible in GTEx data:
% alternatively, maybe we should check their presence in the median
% netowrk versus the individual networks! 
% TODO also representation of functions for each of the SC
% nets. brain and general. 

% CONCLUSIONS1 so far: it would be really nice if the med network
% had more links for brain related functions, but it has < 1/4
% links for the brain related functions and < 1/3 links for the
% general functions. So far I don't know how to interpret that. 
% CONCLUSION2 CTC drops many links from brain specific functions
% and gains some links from the other functions. This result was
% also present in the preliminary results and protein-complex
% analysis. 
% CONCLUSION3 Even after the drop, ctc has more links than the med
% network. I know that med network has highest count of links for
% the most relevant brain functions: CHECK THEM FOR differential
% expression! But the bare conclusion is that, sc has 85% of the
% links in the bulk (both specific and non-specific). So basically,
% what is happening is that SC just has lower rep of functions, all
% the functions. 

% TODO1 draft check:
clear
dataFolder = ['~/data/brainSingleCell/']
load('~/data/brainSingleCell/filDataSet_intronAndExon_V4.mat')
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exonAndIntron_V4.mat'])
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exonAndIntron_V4_clusterLabels.mat'])

load(['~/resultsAndFigures/secondProject/moduleRepresentation/' ...
      'funcRes.mat'])

% have the functions 
load('~/data/general/GOdata_GTExV6_nonComp.mat')
GTExGO = GOdata;

wholeFMat = GTExGO.matP;
% terms with <=200 and >= 5 genes
sib = sum(wholeFMat);
inTerms = ((sib >=5) + (sib<=200)) == 2;

GTExExpGenes = GTExFiveNets.nets(2).expGenes;
fMat = wholeFMat(GTExExpGenes, inTerms);

load(['~/resultsAndFigures/secondProject/moduleRepresentation/' ...
      'moduleRep_GTExBrainCortex_CTC.mat']) % gtexBinNetCTC
load(['~/resultsAndFigures/secondProject/' ...
      'moduleRep_GTExBrainCortex.mat']) % gtexBinNet
load(['~/resultsAndFigures/secondProject/' ...
      'moduleRep_GTExBlood.mat']) % gtexBinNetBlood
load(['~/resultsAndFigures/secondProject/moduleRepresentation/' ...
      'moduleRep_GTExLiver.mat']) % gtexBinNetLiver
load('~/resultsAndFigures/secondProject/moduleRep_AffyTSN.mat') % tsnFen
load('~/resultsAndFigures/secondProject/moduleRep_AffyTAN.mat') % tanFen

load(['~/resultsAndFigures/secondProject/moduleRepresentation/' ...
      'medianNetworks_31Clusters.mat'])

load(['~/networks/allenBrainSC/binNets_exon_V4_fdr_net1.mat'])

% >> marking the terms based on presence in the networks
uniqueTermNames = funcRes.uniqueTermNames;
uniqueTermIDs = funcRes.uniqueTermIDs;
funPresence = zeros(6, length(uniqueTermNames));

% tan
thisPs = tanFen.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < .1;

thisIds = tanFen.inTermsGOID(selected);
[a, b] = ismember(thisIds, uniqueTermIDs);
funPresence(1, b) = 1;

% tsn
thisPs = tsnFen.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < .1;

thisIds = tsnFen.inTermsGOID(selected);
[a, b] = ismember(thisIds, uniqueTermIDs);
funPresence(2, b) = 1;

% gtexBinNet
thisPs = gtexBinNet.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < .1;

thisIds = gtexBinNet.inTermsGOID(selected);
[a, b] = ismember(thisIds, uniqueTermIDs);
funPresence(3, b) = 1;

% gtexBinNetCTC
thisPs = gtexBinNetCTC.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < .1;

thisIds = gtexBinNetCTC.inTermsGOID(selected);
[a, b] = ismember(thisIds, uniqueTermIDs);
funPresence(4, b) = 1;

% gtexBinNetBlood
thisPs = gtexBinNetBlood.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < .1;

thisIds = gtexBinNetBlood.inTermsGOID(selected);
[a, b] = ismember(thisIds, uniqueTermIDs);
funPresence(5, b) = 1;

% gtexBinNetLiver
thisPs = gtexBinNetLiver.ps;
thisFDR = thisPs * length(thisPs);
selected = thisFDR < .1;

thisIds = gtexBinNetBlood.inTermsGOID(selected);
[a, b] = ismember(thisIds, uniqueTermIDs);
funPresence(6, b) = 1;

myNet = medNets.corrected;
myNet = medNets.nonNeuronsRemoved;
totalLinkCount = sum(myNet(:))
% getting FMAT: selecte the GTEx genes
[a, b] = ismember(GTExGO.geneSymbols, medNets.geneSymbols);
myNet = myNet(b(a), b(a));
tempFMat = GTExGO.matP(a, :);
% getting FMAT: selecte the terms
[a, b] = ismember(GTExGO.GOID, funcRes.uniqueTermIDs);
inTermNames = GTExGO.GOTerms(a);
inTermsGOID = GTExGO.GOID(a);
fMat = tempFMat(:, a);
totalPotLinkCount = size(myNet,1) * (size(myNet,1)-1)/2
for i = 1:size(fMat, 2)
    i
    thisTermGenes = logical(fMat(:, i));
    smallNet = myNet(thisTermGenes, thisTermGenes);
    expGeneCounts(i) = size(smallNet, 1);
    linkCounts(i) = sum(smallNet(:));
    potLinkCounts(i) = size(smallNet, 1) * (size(smallNet, 1) -1)/2;
    ps(i) = 1 - hygecdf(linkCounts(i) - 1,  totalPotLinkCount, ...
                        totalLinkCount, potLinkCounts(i));
end

medNetFenNNR.inTermNames = inTermNames;
medNetFenNNR.inTermsGOID = inTermsGOID;
medNetFenNNR.rawGeneCount = rawGeneCount;
medNetFenNNR.expGeneCounts = expGeneCounts;
medNetFenNNR.linkCounts = linkCounts;
medNetFenNNR.potLinkCounts = potLinkCounts;
medNetFenNNR.ps = ps;

save(['~/resultsAndFigures/secondProject/' ...
      'moduleRep_medinaCorrectedNNR.mat'], 'medNetFenNNR')

load(['~/resultsAndFigures/secondProject/' ...
      'moduleRep_medinaCorrectedNNR.mat'])


load('~/networks/allenBrainSC/binNets_exon_V4.mat')

load(['~/networks/allenBrainSC/binNets_exon_V4_fdr_net1.mat'])
% now for the SC 01 
myNet = net.net2;
totalLinkCount = sum(myNet(:))
% getting FMAT: selecte the GTEx genes
[a, b] = ismember(GTExGO.geneSymbols, filDataSet.geneSyms(net.expGenes));
myNet = myNet(b(a), b(a));
tempFMat = GTExGO.matP(a, :);
% getting FMAT: selecte the terms
[a, b] = ismember(GTExGO.GOID, funcRes.uniqueTermIDs);
inTermNames = GTExGO.GOTerms(a);
inTermsGOID = GTExGO.GOID(a);
fMat = tempFMat(:, a);
totalPotLinkCount = size(myNet,1) * (size(myNet,1)-1)/2
for i = 1:size(fMat, 2)
    i
    thisTermGenes = logical(fMat(:, i));
    smallNet = myNet(thisTermGenes, thisTermGenes);
    expGeneCounts(i) = size(smallNet, 1);
    linkCounts(i) = sum(smallNet(:));
    potLinkCounts(i) = size(smallNet, 1) * (size(smallNet, 1) -1)/2;
    ps(i) = 1 - hygecdf(linkCounts(i) - 1,  totalPotLinkCount, ...
                        totalLinkCount, potLinkCounts(i));
end

netSC01.inTermNames = inTermNames;
netSC01.inTermsGOID = inTermsGOID;
netSC01.rawGeneCount = rawGeneCount;
netSC01.expGeneCounts = expGeneCounts;
netSC01.linkCounts = linkCounts;
netSC01.potLinkCounts = potLinkCounts;
netSC01.ps = ps;

save(['~/resultsAndFigures/secondProject/' ...
      'moduleRep_netSC01_exon_V4.mat'], 'netSC01')

load(['~/resultsAndFigures/secondProject/' ...
      'moduleRep_netSC01_exon_V4.mat'])

load('~/data/brainSingleCell/filDataSet_exon_V4.mat')
% >>> getting all enriched terms in all the SC data
for i = 1:69
    i
    load(sprintf(['~/networks/allenBrainSC/allFiveNets/' ...
                  'binNets_Exon_V4_allFive_net%d.mat'], i)) 
    myNet = net.net005;
    totalLinkCount = sum(myNet(:))
    % getting FMAT: selecte the GTEx genes
    [a, b] = ismember(GTExGO.geneSymbols, filDataSet.geneSyms(net.expGenes));
    myNet = myNet(b(a), b(a));
    tempFMat = GTExGO.matP(a, :);
    % getting FMAT: selecte the terms
    [a, b] = ismember(GTExGO.GOID, funcRes.uniqueTermIDs);
    inTermNames = GTExGO.GOTerms(a);
    inTermsGOID = GTExGO.GOID(a);
    fMat = tempFMat(:, a);
    totalPotLinkCount = size(myNet,1) * (size(myNet,1)-1)/2
    ps = ones(1, size(fMat, 2)) * -1;
    for j = 1:size(fMat, 2)
        thisTermGenes = logical(fMat(:, j));
        if sum(thisTermGenes) >=  5
            smallNet = myNet(thisTermGenes, thisTermGenes);
            expGeneCounts(j) = size(smallNet, 1);
            linkCounts(j) = sum(smallNet(:));
            potLinkCounts(j) = size(smallNet, 1) * (size(smallNet, 1) -1)/2;
            ps(j) = 1 - hygecdf(linkCounts(j) - 1,  totalPotLinkCount, ...
                                totalLinkCount, potLinkCounts(j));
        end
    end

    netSC(i).inTermNames = inTermNames;
    netSC(i).inTermsGOID = inTermsGOID;
    netSC(i).rawGeneCount = rawGeneCount;
    netSC(i).expGeneCounts = expGeneCounts;
    netSC(i).linkCounts = linkCounts;
    netSC(i).potLinkCounts = potLinkCounts;
    netSC(i).ps = ps;
end

save(['~/resultsAndFigures/secondProject/' ...
      'moduleRep_netSCAll_exon_V4.mat'], 'netSC')

% selecting brain and non-brain functions using funcres. (that is the
% totalF you will be looking at)
% I removedCTC

% > getting the terms that are in GTEx, it matters
[inGTExTerms, b] = ismember(funcRes.uniqueTermIDs, ...
                  gtexBinNetBlood.inTermsGOID);
funPresence = funcRes.funPresence(:,inGTExTerms);
uniqueTermIDs = funcRes.uniqueTermIDs(inGTExTerms);
uniqueTermNames = funcRes.uniqueTermNames(inGTExTerms);

brainFiltered = ((sum(funPresence(1:3, :)) > 0) - ...
    sum(funPresence(5:6, :))) > 0;
brainIDs = uniqueTermIDs(brainFiltered);
brainNames = uniqueTermNames(brainFiltered);

nonSpecificBrain = (((sum(funPresence(5:6, :)) > 0) + ...
                    (sum(funPresence(1:3, :)) > 0))) == 2;
otherIDs = uniqueTermIDs(nonSpecificBrain);
otherNames = uniqueTermNames(nonSpecificBrain);

% what happens to the brain_bulk brain and nonspecific terms
brainTermPs = zeros(5, sum(brainFiltered));
nonBrainTermPs = zeros(5, sum(nonSpecificBrain));

[a, b] = ismember(uniqueTermIDs(brainFiltered), ...
                  tsnFen.inTermsGOID);
brainTermPs(1, :) = tsnFen.ps(b);

[a, b] = ismember(uniqueTermIDs(brainFiltered), ...
                  gtexBinNet.inTermsGOID);
brainTermPs(2, :) = gtexBinNet.ps(b);
brainTermPs(3, :) = gtexBinNetCTC.ps(b);
netSCps = netSC01.ps(inGTExTerms);
medps = medNetFen.ps(inGTExTerms);
mednnrps = medNetFenNNR.ps(inGTExTerms);
brainTermPs(4, :) = netSCps(brainFiltered);
brainTermPs(5, :) = medps(brainFiltered);
brainTermPs(6, :) = mednnrps(brainFiltered);

netSClc = netSC01.linkCounts(inGTExTerms);

plotPs = brainTermPs * 2475;
medSelected = (plotPs(6, :) < .1); % the terms got by the median
                                   % network
medNames = brainNames(plotPs(6, :) < .1);
plotPs(plotPs > .1)= nan;

h = figure
%heatmap(plotPs)
heatmap(medNames, {'tsn', 'bulk', 'CTC', 'SC', 'medNNR'}, plotPs([1 2 3 4 6] , ...
                                                  medSelected))
h = figure
heatmap(1:length(medNames), {'tsn', 'bulk', 'CTC', 'SC', 'medNNR'}, plotPs([1 2 3 4 6] , ...
                                                  medSelected))

% odds of having overlap between CTC and SC
(18/57) * (9/57)

% p of success using bioncdf (since I am suspecting the terms have
% overlap) 
1-binocdf(8, 57, .0499) % this is the p for such observation, if I
                        % get the likelihood ratio it is 2.8 -
                        % however, just thefact that NNR HAS some
                        % functions enriched is sth (potential
                        % test: NNR with removing 3 sets at a time,
                        % how does that perform?)

% we clearly lose more Brain specific terms when going from bulk to
% CTC. We also lose more brain terms when going from medNNR to SC -
% is there any similarity between what we lose between the two?
% either in bulk to CTC or from medNNR to SC? >>> the above p-value
% is that similarity Is there such similarity for the non specific
% brain ?

sum(~isnan(plotPs)')

% also the terms lost by the cell type correction  
% also make the median network without the non-neurons... I think
% these will go away. It is hard to distinguish the neuron types 
[a, b] = ismember(uniqueTermIDs(nonSpecificBrain), ...
                  tsnFen.inTermsGOID);
brainTermPs(1, :) = tsnFen.ps(b);

[a, b] = ismember(uniqueTermIDs(nonSpecificBrain), gtexBinNet.inTermsGOID);
nonBrainTermPs(2, :) = gtexBinNet.ps(b);
nonBrainTermPs(3, :) = gtexBinNetCTC.ps(b);
nonBrainTermPs(4, :) = netSCps(nonSpecificBrain);
nonBrainTermPs(5, :) = medps(nonSpecificBrain);
nonBrainTermPs(6, :) = mednnrps(nonSpecificBrain);

plotPs = nonBrainTermPs * 2475;
plotPs(plotPs > .1)= nan;
heatmap(plotPs)

medSelected = (plotPs(6, :) < .1); % the terms got by the median
                                   % network
medNames = otherNames(plotPs(6, :) < .1);
plotPs(plotPs > .1)= nan;

h = figure
heatmap(plotPs)
heatmap(medNames, {'tsn', 'bulk', 'CTC', 'SC', 'medNNR'}, plotPs([1 2 3 4 6] , ...
                                                  medSelected))
 
(13/133) * (55/ 133)

1 - binocdf(9 , 113, .0404) % still unlikely  p < .05

% for the terms: count of genes, overlaps and expression in each
% dataset. Also, investigate the terms that overlap between tsn, SC
% and CTC. I still don't know why they are in tsn 

fMat = GTExGO.matP;
[a, b] = ismember(uniqueTermIDs(brainFiltered), GTExGO.GOID);
brainFMat = fMat(GTExFiveNets.nets(2).expGenes, b(a));

medSelectedFMat = brainFMat(:, medSelected);
presentGenes =  sum(medSelectedFMat, 2) > 0;
GTExBrainExpGenes = GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(2).expGenes);
presentGenesSyms = GTExBrainExpGenes(presentGenes);

geneMedSFMat = medSelectedFMat(presentGenes, :);
geneOverlap = geneMedSFMat' * geneMedSFMat;
geneSum = sum(geneMedSFMat);
geneOverlapJC = zeros(57, 57);
for i = 1:57
    for j = i:57
        geneOverlapJC(i,j) = geneOverlap(i, j) / (geneSum(i) + ...
                                                  geneSum(j) - geneOverlap(i,j));
                geneOverlapJC(j,i) = geneOverlap(i, j) / (geneSum(i) ...
                                                          + geneSum(j)- geneOverlap(i,j));
    end
end

h = figure
heatmap([1:57], medNames, geneOverlapJC) % most of the functions presenet in all of
                       % them share some overlap with each
                       % other. the functions which lose term from
                       % median to CTC or SC are less likely to
                       % share  terms >> show me the waves of
                       % expression in CTC for these genes that are
                       % silenced between bulk and CTC: the
                       % correlated genes, and if they are markers,
                       % and how they go down after the CTC. I can
                       % show it for some of them: also, show me
                       % the median pattern that caused their
                       % correlation. What happens in TSN? can you
                       % explain it? My only guess is that these
                       % links are not "brain specific". They are
                       % present in brain, but not necessarily
                       % brain specific: do we see the same effect
                       % for the common processes? Just explain to
                       % me the bulk tissue process 

% load the ctc and bulk and med 
load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes.mat'])

load('~/data/GTEx/Brain_Cortex_expGenes.mat')

load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex.mat'])

bulkExp = (log2(dataSet.mat+1));% - result.regOut;
ctcExpTemp = bulkExp - result.regOut;
signs = -1*(ctcExpTemp < 0) + (ctcExpTemp>=0);
ctcLog = (log2(abs(ctcExpTemp+1)));
ctcExp = signs .* ctcLog;

logReg = (log2(abs(result.regOut+1)));
signs = -1*(result.regOut < 0) + (result.regOut>=0);
regExp = logReg .* signs;

gtexGeneSyms = dataSet.genes;

% medMat: for the median matrix 

netNames{1} = clusterMeta.sortedClusterNames(1);
% getting the exp  mat
[a, b] = ismember(clusterMeta.clusters, ...
                  clusterMeta.sortedClusterNames(1));
scExpMat = log2(normExp(:, a)+1);
scGeneSyms = filDataSet.geneSyms;
% for sc
[a, b] = ismember(thisGeneSyms, scGeneSyms);
fexSC = nan(length(thisGeneSyms), size(scExpMat, 2));
fexSC(a, :) = scExpMat(b(a), :);

% for medNet
[a, b] = ismember(thisGeneSyms, scGeneSyms);
fexMed = nan(length(thisGeneSyms), 28);
fexMed(a, :) = medMat(b(a), :);
heatmap(fexMed)

book = nansum(fexSC > 0);

% Now we have the genes: I want, for every function, do the
% clustering and see how they change. Take common functions first,
% the regulation of neurotransmitters: 
f = 1;
thisGeneSyms = fGenes(logical(testFMat(:, f)));
sum(ismember(thisGeneSyms, veryCorrGenesSyms))
f = f + 1
impGenes = [14, 38, 41, 50];
h = figure
hold all
plot((fexGTEx(14, :)))
plot((fexGTEx(38, :)))
plot((fexGTEx(41, :)))
plot((fexGTEx(50, :)))

plot(fexReg(14,:))

h = figure
hold all
plot((fexCTC(14, :)))
plot((fexCTC(38, :)))
plot((fexCTC(41, :)))
plot((fexCTC(50, :)))

% for GTEx
[a, b] = ismember(thisGeneSyms, gtexGeneSyms);
fexGTEx = bulkExp(b, :);
fexReg = result.regOut(b, :);

fexCTC = ctcExp(b, :);
h = figure
heatmap(zscore(fexGTEx')')
kado = zscore(fexGTEx')';
sib = corr(fexGTEx');
net = sib > .68;
heatmap(net + snGTEx)

h = figure
clustergram(fexGTEx', 'Cluster', 'row', 'ColumnPDist', 'correlation', ...
            'Linkage', 'average', 'Colormap', 'bone')

y = pdist(fexGTEx, 'correlation');
z = linkage(y);
h = figure
dendrogram(z, 0)
t = cluster(z, 'cutoff', 1.15);
hist(t, [1:max(t)])

% find the group of 15 gnees (it was cluster 8) that are
% correlated. these genes are fully correlated, no exception. in
% GTEx, find the group of genes that are also correlated wit
% hthem. there are 500 of those genes, fully correlated. I call
% them the veryCorrGenes
% veryCorrGens15 =

%   1x15 cell array

%   Columns 1 through 6

%     {'PPFIA3'}    {'BRSK1'}    {'STX1A'}    {'DAGLA'}    {'STX1B'}    {'SYT3'}

%   Columns 7 through 12

%     {'RAB3A'}    {'SYT7'}    {'SYN2'}    {'VAMP2'}    {'UNC13A'}    {'RIMS3'}

%   Columns 13 through 15

%     {'SNAP47'}    {'DOC2A'}    {'PNKD'}

veryCorrGenes27 = thisGeneSyms((t == 20));

[a, b] = ismember(veryCorrGenes15, gtexGenes)
[a, b] = ismember(veryCorrGenes27, gtexGenes)
[a, b] = ismember(veryCorrGenes15_2, gtexGenes)
% finding the genes correlated with these:
sib = GTExFiveNets.nets(2).net01;
fullSib = sib + sib';

kado = fullSib(b, :);
halva = sum(kado);
extendedVeryCorr = gtexGenes(halva >=19);

veryCorrGenesSyms = gtexGenes(veryCorrGenes);
save('~/resultsAndFigures/secondProject/veryCorrGenesSyms_GTEx.mat', ...
     'veryCorrGenesSyms')

[sa , sb] = sort(t);
h = figure
heatmap(fexGTEx(sb,:))
h = figure
heatmap(sb, sb,(sib(sb, sb) + 0))
h = figure
heatmap((sib(sb, sb) + 0))
h = figure
heatmap(sb, sb, snGTEx(sb, sb)+0)

myInds = find(t == 1)
h = figure
i = 1
plot(zscore(fexGTEx(myInds(i), :)), 'k')
hold on
i = i + 1
hold off

i = 1
plot(zscore(fexReg(myInds(i), :)), 'r')
hold on
i = i + 1


sb(60:65)
h = figure
heatmap(zscore(fexGTEx(sb(60:65),:)')')

h = figure
heatmap(zscore(fexCTC(sb(60:65),:)')')
sibctc = corr(fexCTC');
h = figure
heatmap(sibctc(sb, sb))
h = figure
heatmap(snCTC(sb, sb) + 0)

h = figure
heatmap(snGTEx(sb, sb) + 0)


inds = [5 6 20 21 11 30 24 2 27 9 22 15 17 1 7 13 19 26 25 8 16 4 ...
        12 10 23 14 18 3 28];
h = figure
heatmap(snGTEx + 0)

% for CTC
fexCTC = ctcExp(b(a), :);
h = figure
heatmap(fexCTC)

% for med 
[a, b] = ismember(thisGeneSyms, scGeneSyms);
fexMed = medMat(b(a), :);
h = figure
heatmap(fexMed)

% for sc
fexSC = scExpMat(b(a), :);


% for SC
[a, b] = ismember(thisGeneSyms, scGenes);
if sum(a) < length(a)
    snSC = zeros(size(tempNet));
    snSC(a==0, a==0) = nan;
    snSC(a > 0, a>0) = net.net2(b(a), b(a));
else
    snSC = net.net2(b, b);
end

% for medNet
[a, b] = ismember(thisGeneSyms, medGenes);
if sum(a) < length(a)
    snMed = zeros(size(tempNet));
    snMed(a==0, a==0) = nan;
    snMed(a > 0, a>0) = medNets.nonNeuronsRemoved(b(a), b(a));
else
    snMed = medNets.nonNeuronsRemoved(b, b);
end


testFMat = logical(geneMedSFMat);

% I want the wave and I want to know if they are markers 

thisF = 1;
% give me the link count in each network. Give me the expressed
% gene count in each network.  

h = figure
heatmap(log2(theseTermsLC(2:5,:)))

h = figure
heatmap(bulkExp(testFMat(:, thisF), :))
h = figure
heatmap(ctcExp(testFMat(:, thisF), :))

brainTermLC = zeros(6, sum(brainFiltered));
[a, b] = ismember(uniqueTermIDs(brainFiltered), ...
                  tsnFen.inTermsGOID);
brainTermLC(1, :) = tsnFen.linkCounts(b);

[a, b] = ismember(uniqueTermIDs(brainFiltered), ...
                  gtexBinNet.inTermsGOID);
brainTermLC(2, :) = gtexBinNet.linkCounts(b(a));
brainTermLC(3, :) = gtexBinNetCTC.linkCounts(b(a));
netSCps = netSC01.linkCounts(inGTExTerms);
medps = medNetFen.linkCounts(inGTExTerms);
mednnrps = medNetFenNNR.linkCounts(inGTExTerms);
brainTermLC(4, :) = netSCps(brainFiltered);
brainTermLC(5, :) = medps(brainFiltered);
brainTermLC(6, :) = mednnrps(brainFiltered);

theseTermsLC = brainTermLC(:, medSelected);


nonBrainTermLC = zeros(6, sum(nonSpecificBrain));
[a, b] = ismember(uniqueTermIDs(nonSpecificBrain), ...
                  tsnFen.inTermsGOID);
nonBrainTermLC(1, :) = tsnFen.linkCounts(b);

[a, b] = ismember(uniqueTermIDs(nonSpecificBrain), ...
                  gtexBinNet.inTermsGOID);
nonBrainTermLC(2, :) = gtexBinNet.linkCounts(b);
nonBrainTermLC(3, :) = gtexBinNetCTC.linkCounts(b);
netSCps = netSC01.linkCounts(inGTExTerms);
medps = medNetFen.linkCounts(inGTExTerms);
mednnrps = medNetFenNNR.linkCounts(inGTExTerms);
nonBrainTermLC(4, :) = netSCps(nonSpecificBrain);
nonBrainTermLC(5, :) = medps(nonSpecificBrain);
nonBrainTermLC(6, :) = mednnrps(nonSpecificBrain);

% give me the links in each network: we get the genes, we get the
% small network, we get the links, we get the expressed genes
GTExFiveNets.nets(2).net01 % got the genes
ctc.net01 % got the genes
gtexGenes = GTExFiveNets.uniqueGeneSyms((GTExFiveNets.nets(2).expGenes));

medNets.nonNeuronsRemoved % got the genes
medGenes = medNets.geneSymbols;
net.net2 % genes are in scBinNets.geneSyms
scGenes = medNets.geneSymbols(net.expGenes);

% testFMat is the matrix of genes and functions
fTerms = medNames;
fGenes = presentGenesSyms; 
fCount = size(testFMat, 1);
for f = 1:fCount
    % now the links and subNetwroks 
    thisGeneSyms = fGenes(logical(testFMat(:, f)));
    
    % for GTEx
    [a, b] = ismember(thisGeneSyms, gtexGenes);
    snGTEx = GTExFiveNets.nets(2).net01(b, b);
    tempNet = snGTEx;
    
    % for CTC
    [a, b] = ismember(thisGeneSyms, gtexGenes);
    snCTC = ctc.net01(b, b);
    
    % for SC
    [a, b] = ismember(thisGeneSyms, scGenes);
    if sum(a) < length(a)
        snSC = zeros(size(tempNet));
        snSC(a==0, a==0) = nan;
        snSC(a > 0, a>0) = net.net2(b(a), b(a));
    else
        snSC = net.net2(b, b);
    end
    
    % for medNet
    [a, b] = ismember(thisGeneSyms, medGenes);
    if sum(a) < length(a)
        snMed = zeros(size(tempNet));
        snMed(a==0, a==0) = nan;
        snMed(a > 0, a>0) = medNets.nonNeuronsRemoved(b(a), b(a));
    else
        snMed = medNets.nonNeuronsRemoved(b, b);
    end
end
h = figure
heatmap(snGTEx + 0)
h = figure
heatmap(snCTC + 0)
h = figure
heatmap(snSC + 0)
h = figure
heatmap(snMed + 0)

% get the common Net:

book = snMed + snSC;
halva = (~isnan(sum(book)));
ins = find(halva);

snsGTEx = snGTEx(ins, ins);
snsCTC = snCTC(ins, ins);
snsSC = snSC(ins, ins);
snsMed = snMed(ins, ins);

mymap = [1 1 1; .9 .9 0; 0 .3 .8; 0 .7 0]

plotMat = zeros(size(snsGTEx));
plotMat(snsGTEx >0) = 1;
plotMat(snsCTC > 0) = 2;
plotMat((snsGTEx + snsCTC) == 2) = 3;
filter2 = (sum(plotMat)  + sum(plotMat')) >0;
h = figure
heatmap(plotMat(filter2, filter2), 'colormap', mymap)


plotMat = zeros(size(snsGTEx));
plotMat(snsMed >0) = 1;
plotMat(snsSC > 0) = 2;
plotMat((snsMed + snsSC) == 2) = 3;
%filter2 = (sum(plotMat)  + sum(plotMat')) >0;
h = figure
heatmap(plotMat(filter2, filter2), 'colormap', mymap)

plotMat = zeros(size(snsMed));
plotMat(snsMed >0) = 1;
plotMat(snsSC > 0) = 2;
plotMat((snsMed + snsSC) == 2) = 3;
h = figure
heatmap(plotMat, 'colormap', mymap)


mymap = [1 1 1; .9 .9 0; 0 .3 .8; 0 .7 0]

plotMat = zeros(size(snGTEx));
plotMat(snGTEx >0) = 1;
plotMat(snMed > 0) = 2;
plotMat((snGTEx + snMed) == 2) = 3;
filter2 = (sum(plotMat)  + sum(plotMat')) >0;
h = figure
heatmap(plotMat(filter2, filter2), 'colormap', mymap)

% select the genes in GTEx brain expressed 

% give me the brain terms lost from bulk to CTC, examine their expression
% level in the sc - are they markers?

% from the terms that we lose from bulk to SC, do we also lose them
% in median? and is it more than expected?

% find the gene count, link count and significance of of each of
% the terms in each of the networks 

% get the enrichment for new median and new sc function

expMat = filDataSet.expMat(:, clusterMeta.inCells);
% normalizing for the million count
sumExp = sum(expMat);
milFac = sumExp ./ 1000000;
normExp = zeros(size(expMat));
sampleCount = size(expMat, 2);
for i = 1:sampleCount
    normExp(:, i) = expMat(:, i)./milFac(i);
end

[a, b] = ismember(funcRes.uniqueTermIDs, GTExGO.GOID);
fMat = GTExGO.matP(:, b(a));
fMatBrain = fMat(:, brainFilteredIDs);
sib = sum(fMatBrain);
potLink = sib .* (sib - 1) ./2;

fMatBrain = fMat(:, nonSpecificBrain);
sib = sum(fMatBrain);
potLink = sib .* (sib - 1) ./2;

[a_med, b_med] = sort(brainLinkCounts_med, 'descend');
myTerms = funcRes.uniqueTermNames(brainFilteredIDs);
myTermsSorted = myTerms(b_med);

% getting the ratio, the genes and the expression. 
for i = 1:length(myTerms)
    % get the GTEx genes from fMat
    myGenes = fMatBrain(:, b_med(i));
    myGeneSymbols = GTExGO.geneSymbols(logical(myGenes));
    [a, b] = ismember(myGeneSymbols, filDataSet.geneSyms);
    scGenesInds = b(a);
    myGenesCount = length(scGenesInds);
    potLinkCount = myGenesCount * (myGenesCount - 1)/2;
    r = a_med(i) / potLinkCount
    %    sampleCount = sum(clusterMeta.sortedObjCount(1:31))
    myExpMat = zeros(myGenesCount, 31*50);
    c = 1;
    for j = 1:inClusterCount
        
        j
        netNames{j} = clusterMeta.sortedClusterNames(j);
        % getting the exp  mat
        [a, b] = ismember(clusterMeta.clusters, ...
                          clusterMeta.sortedClusterNames(j));
        thisExpMat = log2(normExp(:, a)+1);
        
        randomSamples = datasample(1:size(thisExpMat, 2), 50);
        myExpMat(:, c:(c +49)) = ...
            thisExpMat(scGenesInds, randomSamples);
        c = c + 50;
    end
end

myMat = medMat(scGenesInds, :);
kado = std(myMat);
heatmap(myMat(:, [7,9,14:16, 19, 20, 23]))
heatmap(corr(medMat(scGenesInds, :)))
% get the quantiles for the median network and things 
% I am looking at 31 clusters > those with > 100 samples
inClusterCount = 31;
gCount = size(expMat, 1);
netNames = cell(1, 31);
medMat = zeros(gCount, inClusterCount);
qs = zeros(100, 31);
for i = 1:inClusterCount
    i
    netNames{i} = clusterMeta.sortedClusterNames(i);
    % getting the exp  mat
    [a, b] = ismember(clusterMeta.clusters, ...
                      clusterMeta.sortedClusterNames(i));
    thisExpMat = log2(normExp(:, a)+1);
    
    medMat(:, i) = median(thisExpMat');
    
    % sib = corr(thisExpMat');
    % upCorr = sib(logical(triu(ones(size(sib)), 1)));
    % qs(:, i) = quantile(upCorr, [0:.01:.99]);
end

heatmap(netNames, myGenes, medMat)
sib = corr(smallMedMat(:, 1:31)') - eye(6,6);
sib = corr(smallExpMat(:, 1:30)') - eye(6,6);

heatmap(smallMedMat)
h = figure
heatmap(sib)

% 10. Getting the modules which do not have overlap 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('~/data/general/GOdata_GTExV6_nonComp.mat')
GTExGO = GOdata;
clear GOdata;
load('~/data/GTEx/Brain_Cortex_expGenes.mat')

[a, b] = ismember(dataSet.genes, GTExGO.geneSymbols);
wholeFMat = GTExGO.matP(b, :);

inTerms = ((sib >=20) + (sib<=200)) == 2;
fMat = wholeFMat(:, inTerms);
sib = sum(fMat);

inTermNames = GTExGO.GOTerms(inTerms);
inTermsGOID = GTExGO.GOID(inTerms);

overlap = fMat' * fMat;
max(overlap(:))
myD = diag(overlap);
normMat = zeros(size(overlap));
for i = 1:length(normMat)
    i
    for j = i + 1:length(normMat)
        n = min(myD([i, j]));
        normMat(i, j) = overlap(i, j) / n;
        normMat(j, i) = overlap(i, j) / n;
    end
end
h = figure
hist(normMat(:))

lowOverlap = (normMat <= .05);
h = figure
heatmap(normMat(1:10, 1:10))

% get the meanR for each term 
load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_AllenSCBasedMarker.mat'])
fCount = length(normMat);
fTermMeanR = zeros(1, fCount);
for i = 1:fCount
    fTermMeanR(i) = mean(result.r2(logical(fMat(:, i))));
end
h = figure
hist(fTermMeanR)
lowMeanR = fTermMeanR < mean(fTermMeanR);
highMeanR = fTermMeanR > mean(fTermMeanR);

% get the link overlap for each term in GTEx brain
load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
 ...
      'correctedBinNets_logCorrected_jusResiduals_scBasedMarkers_7PCA.mat'])
ctcNet = ctc.net005 + ctc.net005';

load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat')
gNet = GTExFiveNets.nets(2).net005 + GTExFiveNets.nets(2).net005';

lcCountsG = zeros(fCount, fCount);
lcPortionG = zeros(fCount, fCount);
lcCountsCTC = zeros(fCount, fCount);
lcPortionCTC = zeros(fCount, fCount);
for i = 1:fCount
    i
    tic
    iInds = logical(fMat(:, i));
    iCount = sum(iInds);
    
    lcCountsG(i, i) = sum(sum(gNet(iInds, iInds)));
    lcCountsCTC(i, i) = sum(sum(ctcNet(iInds, iInds)));
    
    lcPortionG(i, i) = lcCountsG(i, i) / (iCount * (iCount-1));
    lcPortionCTC(i, i) = lcCountsCTC(i, i) / (iCount * (iCount-1));
    for j = i + 1:fCount
        jInds = logical(fMat(:, j));
        jCount = sum(jInds);
        
        lcCountsG(j, j) = sum(sum(gNet(jInds, jInds)));
        lcCountsCTC(j, j) = sum(sum(ctcNet(jInds, jInds)));
        
        lcCountsG(i, j) = sum(sum(gNet(iInds, jInds)));
        lcCountsCTC(i, j) = sum(sum(ctcNet(iInds, jInds)));
        
        lcPortionCTC(j, j) = lcCountsCTC(j, j) / (jCount * (jCount-1));
        lcPortionCTC(i, j) = lcCountsCTC(i, j) / (iCount * (jCount));
        
        lcPortionG(j, j) = lcCountsG(j, j) / (jCount * (jCount-1));
        lcPortionG(i, j) = lcCountsG(i, j) / (iCount * (jCount));
    end
    toc
end

sum(diag(lcCountsCTC(lowMeanR, lowMeanR))) /sum(diag(lcCountsG(lowMeanR, lowMeanR)))

% the intra cluster links between lowMeanR increased by 9.5%
% the inter cluster links between lowMeanR increased by 

lCTC = lcCountsCTC(lowMeanR, lowMeanR);
lowOverlapLowMeanR = lowOverlap(lowMeanR, lowMeanR);
lG = lcCountsG(lowMeanR, lowMeanR); 
sum(sum(triu(lCTC, 1) .* (lowOverlapLowMeanR+0))) / sum(sum(triu(lG, ...
                                                  1) .*(lowOverlapLowMeanR+0)))
sum(sum(triu(lCTC, 1))) / sum(sum(triu(lG, 1)))

% for those that were increased, how much were their neighbours
% lc

lcPortionGMod = lcPortionG;
lcPortionGMod(lcPortionGMod == 0) = -1;
porR = lcPortionCTC ./ lcPortionGMod;
porRD = diag(porR);

lowMeanRInds = find(lowMeanR);
highMeanRInds = find(highMeanR);

% getting the index of the tese terms with high increase
selected = (lowMeanR + (porRD>=2)') == 2;
selectedIDs = find(selected);

% now for these terms, I want the increased portion of overlap with
% terms which do not overlap with it
smallMat = porR(selected, :);
kado = zeros(size(smallMat));

tempInds = lowOverlap(selected, :);
kado(~tempInds) = nan;
kado(lowOverlap(selected, :)) = ...
    smallMat(lowOverlap(selected, :));

% it is enough to say that the link increase is not random, and is
% selective towards some processes. 


% now, I have for the terms with lowR2 (< mean) , with high
% Increase ratio (>= 2), the overlap ratio of themselves and
% neighbors (kado) and their index (highIncLCRatioLRInds);

myMat = kado;
myInds = highIncLCRatioLRInds;


porRDH = porRD(highMeanR);
h = figure
boxplot(porRD(lowMeanR))
ylim([0 15])

h = figure
boxplot(porRD(highMeanR))
ylim([0 15])


% get the link overlap for each term in GTEx ctc 
% do we see more intra function links than extra function links? in
% short, are the links pulled, or they are distributed all over?
% There must be a total count of link increase. Is it distributed
% inter cluster/intra cluster similarly, differnetly? for the
% clusters that are this or that. 



