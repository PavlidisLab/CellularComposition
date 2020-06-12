% This file has codes for obtaining specific results related to
% the manuscript. It has anything from late stage plots to
% modifications based on the review comments. 
% 
% 1. example for genes associated with brain specific functional
% terms (went into Figure 01)
% 2. the correlation of the CT profiles. This is to get rid of the
% surrogate CT profiles
% 3. writing the marker genes in the file
% 4. for genes, their cluster and their R2 from the model 
% 5. for clusters, the ID and average R2
% 6. for clusters, their genes, their functional terms
% 7. writing the high rep links into file
% 8. writing the list of enriched functions for high rep clusters
% into file
% 9. distribution of quasi markers
% 10. likelihood ratio of the overlap of links in the networks
% 11. Table for count of genes in the clusters
% 12. example of links for last result part
% 13. overlap between clusters: snuc network versus GTEx and CTC
% 14. overlapping links  
% 15. contamination in snuc-RNA seq data
% 17. removing the mito from the heatmap
% 18. just plot the R2
% 19. writhe the functional enrichment
% 20. snuc networks and their link counts
% 21. Review: correlation of genes in snuc one dataSet one cluster
% 22. Review: correlation of PCs AND marker expression
% 23. Review: the marker counts - comment 31
% 24. Review: count of genes in marker enriched cluster group
% 25. Review: if there are marker genes in the highRep clusters
% 26. Review: 5 versus 7 pcs and dists
% 27. that R2 thing

% 1. example for genes associated with brain specific functional terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load bulk data, functional terms association, the snc data
load('~/data/brainSingleCell/filDataSet_exon_V4.mat')
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exon_V4.mat'])
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exon_V4_clusterLabels.mat'])

load('~/data/GTEx/Brain_Cortex_expGenes.mat')

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_AllenSCBasedMarker_redo.mat'])

load('~/data/general/GOdata_GTExV6_nonComp.mat')
GTExGO = GOdata;
clear GOdata;

% fetch the genes
%1 
[a, b] = ismember(16575, GTExGO.GOID)
GTExGO.GOTerms(b)
[a, b] = ismember(31935, GTExGO.GOID)
GTExGO.GOTerms(b)
[a, b] = ismember(387, GTExGO.GOID)
GTExGO.GOTerms(b)
[a, b] = ismember(16577, GTExGO.GOID)
GTExGO.GOTerms(b)
geneList= GTExGO.geneSymbols(logical(GTExGO.matP(:, b)));
whos geneList
[a, b] = ismember(dataSet.genes, geneList);
mean(result.r2(a))

% 2
[a, b] = ismember(48167, GTExGO.GOID)
GTExGO.GOTerms(b)
[a, b] = ismember(8366, GTExGO.GOID)
GTExGO.GOTerms(b)
geneList= GTExGO.geneSymbols(logical(GTExGO.matP(:, b)));
whos geneList
[a, b] = ismember(dataSet.genes, geneList);
mean(result.r2(a))

% fetch the variance in the snc data for these genes
[a, b] = ismember(scExpGeneList, geneList);

book = thisExpMat(a, :);
inGenes = sum(book'>0) >= 100;

filteredBook = book(inGenes, :);

finalDS = zeros(sum(inGenes), 100);
for i = 1:sum(inGenes)
    mySelect = datasample(find(filteredBook(i, :) > 3),100);
    finalDS(i, :) = filteredBook(i, mySelect);
end

h = figure
heatmap(finalDS)
h = figure
hist(var(finalDS'))
vs1 = var(finalDS')
vs2 = var(finalDS')

% fetch the variance in the bulk tissue for these genes

bulkMat = log2(dataSet.mat);
[a, b] = ismember(dataSet.genes, geneList);

h = figure
hist(var(bulkMat(a, :)'))
vb1 = var(bulkMat(a, :)')
vb2 = var(bulkMat(a, :)')

% can we explain the variance in the bulk with high R2? is the
% variance in the snc lower?
% doing a boxplot

bulkVarGroup = [ones(1, length(vb1)), 2*ones(1, length(vb2))]
bulkVarData = [vb1, vb2];

sinVarGroup =  [ones(1, length(vs1)), 2*ones(1, length(vs2))]
sinVarData =  [vs1, vs2];

h = figure
boxplot(sinVarData, sinVarGroup, 'Labels', {'histone demethylation', ...
                    'Regulation of synaptic plasticity'})
title(['Observed variance in the single cell data - excitatory ' ...
       'cells'])
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%svarianceExample_SC', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

h = figure
boxplot(bulkVarData, bulkVarGroup, 'Labels', {'histone demethylation', ...
                    'Regulation of synaptic plasticity'})
title('Observed variance in the bulk tissue')
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%svarianceExample_bulk', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 2. the correlation of the CT profiles. This is to get rid of the
% surrogate CT profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thePattern = zscore(geneListExpMat')';

load('~/resultsAndFigures/secondProject/gtexClusters.mat')

myIDs = [231, 244, 57, 45, 107, 35, 53, 92, 161,242 239, 243, ...
         87, 136, 232, 174, 129]

clusterPatternSig = zeros(length(myIDs), 75);
for i = 1:length(myIDs)
    thisGenes = gtexCluster.syms(gtexCluster.cs1547 == myIDs(i));
    [a, b] = ismember(thisGenes, scExpGeneList);
    smallMat = thePattern(b(a), :);
    %    smallMat = geneListExpMat(b(a), :);
    clusterPatternSig(i, :) = mean(smallMat);
end

IDLabels = cell(1,17);
for i = 1:17
    IDLabels{i} = sprintf('ID%d', myIDs(i));
end

[s, p] = corr(clusterPatternSig');
h = figure
heatmap(IDLabels, IDLabels, s)
colormap(gray)
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sclusterCTProfiles_CorrelationHeatmap', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 3. writing the marker genes in the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% markers from Mancarci et al
load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])

fileName = ['~/resultsAndFigures/secondProject/suppFiles/' ...
            'Mancarci_filtered_Markers.txt'];
fid = fopen(fileName, 'w')
fprintf(fid, ['Astrocyte\tEndothelial\tMicroglia\tOligo\tPyramidal\' ...
              'n'])
markers = [1 2 5 8 10]
for i = 1:92
    m1 = markers.genes(markers.finalMarkerGenes{1}(i));
    m1 = m1{1}
    
    if i < 25
        m2 = markers.genes(markers.finalMarkerGenes{2}(i));
        m2 = m2{1}
    else
        m2 = ''
    end
    
    if i < 50
        m3 = markers.genes(markers.finalMarkerGenes{5}(i));
        m3 = m3{1}
    else 
        m3 = ''
    end
    
    if i < 72
        m4 = markers.genes(markers.finalMarkerGenes{8}(i));
        m4 = m4{1}
    else 
        m4 = ''
    end
    
    if i < 25
        m5 = markers.genes(markers.finalMarkerGenes{10}(i));
        m5 = m5{1}
    else 
        m5 = ''
    end
    
    fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', m1, m2, m3, m4, m5)
end

% Markers from SNC
load('~/resultsAndFigures/secondProject/scBasedMarkers.mat')

fileName = ['~/resultsAndFigures/secondProject/suppFiles/' ...
            'snuc_Markers.txt'];
fid = fopen(fileName, 'w')
fprintf(fid, ['Astrocyte\tEndothelial\tMicroglia\tOligo\tPyramidal\' ...
              'n'])
1 2 3 5 6
indLists = zeros(180, 5);
for i = [1 2 3 5 6]
    i
    kado = find(scBasedMarkers.markerMat(:, i));
    indLists(1:length(kado), i) = kado;
end
indLists(:, 4) = [];

for i = 1:size(indLists, 1)
    
    if i < 180
        m1 = scBasedMarkers.genes(indLists(i, 1));
        m1 = m1{1};
    else 
        m1 = ''
    end
    
    if i < 133
        m2 = scBasedMarkers.genes(indLists(i, 2));
        m2 = m2{1}
    else
        m2 = ''
    end
    
    if i < 127
        m3 = scBasedMarkers.genes(indLists(i, 3));
        m3 = m3{1}
    else 
        m3 = ''
    end
    
    if i < 142
        m4 = scBasedMarkers.genes(indLists(i, 4));
        m4 = m4{1}
    else 
        m4 = ''
    end
    
    m5 = scBasedMarkers.genes(indLists(i, 5));
    m5 = m5{1}
    
    fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', m1, m2, m3, m4, m5)
end

% 4. for genes, their cluster and their R2 from the model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear 
load('~/resultsAndFigures/secondProject/gtexClusters.mat')


load('~/data/GTEx/Brain_Cortex_expGenes.mat')

[a,b] = ismember(gtexCluster.syms, dataSet.genes);
myGeneR2s = result.r2(b);

fileName = ['~/resultsAndFigures/secondProject/suppFiles/' ...
            'GTExBulk_geneClusterAndR2_redo.txt'];
fid = fopen(fileName, 'w')
fprintf(fid, ['geneSymbol\tclusterID\tR2\n'])

for i = 1:length(gtexCluster.syms)
    fprintf(fid, '%s\t%d\t%.3f\n', gtexCluster.syms{i}, gtexCluster.cs1547(i), ...
            myGeneR2s(i));
end


% 5. for clusters, the ID and average R2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
load('~/data/GTEx/Brain_Cortex_expGenes.mat')
load('~/resultsAndFigures/secondProject/gtexClusters.mat')

clusterR2s = zeros(1, 253);
for i = 1:253
    inGenes = gtexCluster.syms(gtexCluster.cs1547 == i);
    [a, b] = ismember(dataSet.genes, inGenes);
    clusterR2s(i) = mean(result.r2(a));
end

fileName = ['~/resultsAndFigures/secondProject/suppFiles/redo/' ...
            'clusterInfo_redo.txt'];

fid = fopen(fileName, 'w')
fprintf(fid, ['clusterID\tmemberCount\taverageR2\n'])

myInds = find(gtexCluster.memberCounts >=20);

for i = 1:length(myInds)
    fprintf(fid, '%d\t%d\t%.2f\n', myInds(i), ...
            gtexCluster.memberCounts(myInds(i)), clusterR2s(myInds(i)))
end

% >>> 5.1 Genome Research comment 39 - marker gene information
% added to the file

% adding count of marker genes for each cluster
% also adding type of cluster
load('~/resultsAndFigures/secondProject/scBasedMarkers.mat')
markerMat = logical(scBasedMarkers.markerMat);

fileName = ['~/resultsAndFigures/secondProject/suppFiles/redo/' ...
            'clusterInfo_redo_GRcomment39.txt'];

fid = fopen(fileName, 'w')
fprintf(fid, ['clusterID\tmemberCount\taverageR2\tmarkerCount_Astrocyte\tmarkerCount_Endothelial\tmarkerCount_Microglia\tmarkerCount_Oligodendrocyte\tmarkerCount_Pyramidal\n'])

myInds = find(gtexCluster.memberCounts >=20);

ps = ones(length(myInds), 5);
for i = 1:length(myInds)
    
    thisClusterMembers = gtexCluster.syms(gtexCluster.cs1547 == myInds(i));

    % >>>> getting marker counts
    [a, b] = ismember(scBasedMarkers.genes(markerMat(:, 1)), thisClusterMembers);
    acount = sum(a);
    ps(i, 1) = 1 - hygecdf(length(a), 12000, sum(markerMat(:, 1)), );


    [a, b] = ismember(scBasedMarkers.genes(markerMat(:, 2)), thisClusterMembers);
    ecount = sum(a);

    [a, b] = ismember(scBasedMarkers.genes(markerMat(:, 3)), thisClusterMembers);
    mcount = sum(a);

    [a, b] = ismember(scBasedMarkers.genes(markerMat(:, 5)), thisClusterMembers);
    ocount = sum(a);

    [a, b] = ismember(scBasedMarkers.genes(markerMat(:, 6)), thisClusterMembers);
    pcount = sum(a);
    
    fprintf(fid, '%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\n', myInds(i), ...
            gtexCluster.memberCounts(myInds(i)), clusterR2s(myInds(i)), ...
            acount, ecount, mcount, ocount, pcount)
end

% 6. for clusters, their genes, their functional terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'functionalEnrichment.mat'])

myMat = clusterFRes.ps * length(clusterFRes.fTerms);
passFDR = ((myMat > 0) + (myMat < .1)) == 2;
myInds = find(sum(passFDR') >=3)

% write files with the clusters 
for i = 1:length(myInds)
    thisCind = myInds(i);
    myFile = sprintf(['~/resultsAndFigures/secondProject/suppFiles/clusterFTerms/' ...
                      'clusterID%d.txt'], clusterFRes.clusterIDs(thisCind));
    fid = fopen(myFile, 'w')
    fprintf(fid, ['GOID\tGOTerm\tpvalue\n'])
    myFInds = find(passFDR(thisCind, :));
    
    % writing the file
    for j = 1:length(myFInds)
        thisFind = myFInds(j);
        fprintf(fid, '%d\t%s\t%e\n', clusterFRes.fIDs(thisFind), ...
                clusterFRes.fTerms{thisFind}, clusterFRes.ps(thisCind, ...
                                                          thisFind));
    end

end

% >>>>>>>>>> Genome biology: merging the files
% write files with the clusters 

myFile = sprintf(['~/resultsAndFigures/secondProject/suppFiles/redo/' ...
                  'clusterFunctionTerm.txt']);
fid = fopen(myFile, 'w')

for i = 1:length(myInds)
    i
    thisCind = myInds(i);
    fprintf(fid, 'cluster ID: %d\n', clusterFRes.clusterIDs(thisCind))
    fprintf(fid, ['GOID\tGOTerm\tpvalue\n'])
    myFInds = find(passFDR(thisCind, :));
    
    % writing the file
    for j = 1:length(myFInds)
        thisFind = myFInds(j);
        fprintf(fid, '%d\t%s\t%e\n', clusterFRes.fIDs(thisFind), ...
                clusterFRes.fTerms{thisFind}, clusterFRes.ps(thisCind, ...
                                                          thisFind));
    end
    
    fprintf(fid, '\n')
end


% >>>>>>>>>> Genome biology: merging the files - second format
% write files with the clusters 

myFile = sprintf(['~/resultsAndFigures/secondProject/suppFiles/redo/' ...
                  'clusterFunctionTerm_secondFormat_try.txt']);
fid = fopen(myFile, 'w')
fprintf(fid, ['clusterID\tGOID\tGOTerm\tpvalue\n'])

originalIDs = find(gtexCluster.memberCounts > 20)

for i = 1:length(myInds)
    i
    thisCind = myInds(i);
    myFInds = find(passFDR(thisCind, :));
    
    % writing the file
    kado = originalIDs(thisCind)
    for j = 1:length(myFInds)
        thisFind = myFInds(j);
        fprintf(fid, '%d\t%d\t%s\t%e\n', kado, clusterFRes.fIDs(thisFind), ...
                clusterFRes.fTerms{thisFind}, clusterFRes.ps(thisCind, ...
                                                          thisFind));
    end
end



% 7. printing the high rep links
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each link, I want to have the R2 of the genes, it's
% coexpression in the GTEx and in the TAN-brain, also the count of
% links in exc and inh
% the file is there, I can also get examples
clear

% load bulk data, functional terms association, the snc data
load('~/data/brainSingleCell/filDataSet_exon_V4.mat')
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exon_V4.mat'])
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exon_V4_clusterLabels.mat'])



% sumNet comes from networksOverlapLR.m
load(['~/resultsAndFigures/secondProject/' ...
      'overlappingLinksBetweenSCNets/sumNets.mat'])
sumNet = sumNets.exc + sumNets.inh;
sumNet = triu(sumNet);
scGenes = filDataSet.geneSyms;

[a, b, c] = find(sumNet);
inLinks = c >=  10;
inLinks = c >=  5; % the file is printed for 5
inC = c(inLinks);
inA = a(inLinks);
inB = b(inLinks);

sum(sum(sumNet >= 5))

% get the GTEx net

% get the TAN
load( ['~/networks/tissues/brain/' ...
       'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])

% get the exc and inh counts
% >> just get it 
load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_AllenSCBasedMarker_redo.mat'])

% get the GTEx corr ranks 
load('~/data/GTEx/Brain_Cortex_expGenes.mat') 

load('~/data/general/GPL570GemmaMapNEW.mat')
affyGeneSyms = gpl570.uniqueSymbols;
clear gpl570

load('~/resultsAndFigures/secondProject/gtexClusters.mat')


[a, b]= ismember(dataSet.genes, gtexCluster.syms);
sib = corr(log2(dataSet.mat'+1));
upCorr = sib(logical(triu(ones(size(sib)), 1)));
qs = quantile(upCorr, 999);
sib = sib(a, a);
gtexGenes = dataSet.genes(a);
r2s = result.r2(a);

% getting the quantiles for the coexpression values
myCRank = zeros(1, length(inC));
tanP = zeros(1, length(inC));
R1 = zeros(1, length(inC));
R2 = zeros(1, length(inC));
for i = 1:length(inC)
    i
    [a, b1] = ismember(scGenes(inA(i)), gtexGenes);
    [a, b2] = ismember(scGenes(inB(i)), gtexGenes);

    if (b1 > 0)&&(b2 >0)
        myCorr = sib(b1, b2);
        if myCorr < qs(999)
            myCRank(i) = min(find(qs > myCorr));
        else
            myCRank(i) = qs(999);
        end
        R1(i) = result.r2(b1);
        R2(i) = result.r2(b2);
    end
    [a, b1] = ismember(scGenes(inA(i)), affyGeneSyms);
    [a, b2] = ismember(scGenes(inB(i)), affyGeneSyms);
    if (b1 > 0)&&(b2 >0)
        tanP(i) = binNet(b1, b2);
    end
end

% getting the regression values for the genes

% write the file
fileName = ['~/resultsAndFigures/secondProject/suppFiles/' ...
            'sc_robustLinks.txt'];

fid = fopen(fileName, 'w')
fprintf(fid, ['gene_1\tgene_2\ttotalRep\tinhRep\texcRep\tGTExCorrRank\tTAN_Presence\tgene1_r2\tgene2_r2\n'])

excC = zeros(1, length(inC));
inhC = zeros(1, length(inC));
for i = 1:length(myCRank)
    i
    gene1 = scGenes(inA(i));
    gene2 = scGenes(inB(i));
    excC(i) = full(sumNets.inh(inA(i), inB(i)));
    inhC(i) = full(sumNets.exc(inA(i), inB(i)));
    fprintf(fid, '%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\n', ...
            gene1{1}, gene2{1}, inC(i), ...
            excC(i), inhC(i),...
            myCRank(i), tanP(i), R1(i), R2(i));
end

h = figure
scatter(inC, myCRank)

maxRs = max([R1; R2]);

book = (((inC' >= 15) + (myCRank <= 800) + (maxRs >=.6)) == 3);
sum(book)

halva = [inC(book)'; myCRank(book); maxRs(book)];
finalR1 = R1(book);
finalR2 = R2(book);
finalIndsA = inA(book);
finalIndsB = inB(book);

r1s = R1(book);
r2s = R2(book);
g1s = scGenes(inA(book));
g2s = scGenes(inB(book));
incs = inC(book);
inhcs = inhC(book);
exccs = excC(book);

i = 10
r1s(i)
r2s(i)
g1s(i)
g2s(i)
incs(i)
inhcs(i)
exccs(i)

geneList = [g1s(i), g2s(i)]

maxRs = max([R1; R2]);

h = figure
scatter(R1, R2, 'filled')
alpha(.05)

h = figure
scatter(r1s, r2s, 'filled')
alpha(.05)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8. also writing list of genes
% 8. writing the list of enriched functions for high rep clusters
% into file

load(['~/resultsAndFigures/secondProject/' ...
      'overlappingLinksBetweenSCNets/highRepLinks.mat'])

fileName = ['~/resultsAndFigures/secondProject/suppFiles/' ...
            'clusterInfo_SCRobust.txt'];

fid = fopen(fileName, 'w')
fprintf(fid, ['gene\tcluster\n'])

[a, b] = ismember(overlapRes.clusters, [1 9 10]);
myGenes = overlapRes.geneSyms(a);
myClusters = overlapRes.clusters(a);
aind = find(a);

colors = {'pink', 'grey', 'black'}
for i = 1:length(aind)
    thisCluster = myClusters(i);
    thisColor = '';
    switch thisCluster
      case 9
        thisColor = colors(2)
      case 10
        thisColor = colors(3)
      case 1
        thisColor = colors(1)
    end
    thisGene = myGenes{i};
    fprintf(fid, '%s\t%s\n', thisGene, ...
            thisColor{1})
end

%%% writing the functions
myFile = ['~/resultsAndFigures/secondProject/suppFiles/clusterFTerms/' ...
         'snuc_RNAseqCluster_black.txt'];
fid = fopen(myFile, 'w')
fprintf(fid, ['GOID\tGOTerm\tFDR\n'])

myGenes = overlapRes.gtexGeneSyms(overlapRes.clusters == 10);
[pfdrIDs, pfdrTerms, inGenes, enrichMat, pfdr] = fenrich_function(myGenes, ...
                                                  .1, 16000);
% writing the file
for i = 1:length(pfdrIDs)
    fprintf(fid, '%d\t%s\t%e\n', pfdrIDs(i), pfdrTerms{i}, full(pfdr(i)));
end

% 9. distribution of quasi markers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for list of astrocyte and pyramidal clusters, do the distribution
% of the expression level of the cluster genes in them

% get the clusters
astroIDs = [231, 57, 45, 107, 35, 53, 92]

pyraIDs = [161,242 239, 243, ...
         87, 136, 232, 174, 129]

oligoIDs = [252]

endoIDs = [98]

microIDs = [65]

% get the genes 
[a, b] = ismember(gtexCluster.cs1547, astroIDs);
astroGenes = gtexCluster.syms(a);

[a, b] = ismember(gtexCluster.cs1547, pyraIDs);
pyraGenes = gtexCluster.syms(a);

[a, b] = ismember(gtexCluster.cs1547, endoIDs);
endoGenes = gtexCluster.syms(a);

[a, b] = ismember(gtexCluster.cs1547, oligoIDs);
oligoGenes = gtexCluster.syms(a);

[a, b] = ismember(gtexCluster.cs1547, microIDs);
microGenes = gtexCluster.syms(a);

[a, b] = ismember(astroGenes, filDataSet.geneSyms);
astroSCInds = b(a);

[a, b] = ismember(pyraGenes, filDataSet.geneSyms);
pyraSCInds = b(a);

[a, b] = ismember(microGenes, filDataSet.geneSyms);
microSCInds = b(a);

[a, b] = ismember(endoGenes, filDataSet.geneSyms);
endoSCInds = b(a);

[a, b] = ismember(oligoGenes, filDataSet.geneSyms);
oligoSCInds = b(a);

total = [astroSCInds, pyraSCInds];
groups = [ones(1, length(astroSCInds)), ones(1, length(pyraSCInds)).*2];

% do the distribution of the gene expression in the two clusters

% give me the ID of the cell population, I will give you the plot
thislabs = clusterMeta.sortedClusterNames(1:75);
[a, b] = sort(thislabs(1:75));
i = 1

means = geneListExpMat(total, b(i));
%h = figure
boxplot(means, groups)
book = sprintf('%s', a{i});
title(book)
i = i + 1

% plot all of them together

% all and everything - not using it 
h = figure
boxplot(geneListExpMat(astroSCInds, b), 'PlotStyle', 'compact')

h = figure
boxplot(geneListExpMat(pyraSCInds, b), 'PlotStyle', 'compact')

% the exc and astro populations >>>>>>>
thislabs = clusterMeta.sortedClusterNames;
[a, b] = sort(thislabs(1:75));
cellPops = b([1, 2, 4:27]);
xlabels = thislabs(b([1, 2, 4:27]));

cellPops = b;
xlabels = thislabs(b);

% EXC
x = 1:length(cellPops);
y =  mean(geneListExpMat(pyraSCInds, cellPops));
ypos = std(geneListExpMat(pyraSCInds, cellPops))/2;
yneg = zeros(1, length(x));
h = figure('units', 'centimeters', 'position', [0,0, 30, 12])
errorbar(x, y, yneg, ypos, 'o')

x = x + .2;
y = mean(geneListExpMat(logical(scBasedMarkers.markerMat(:,6)), ...
                        cellPops));
ypos = std(geneListExpMat(logical(scBasedMarkers.markerMat(:,6)), ...
                        cellPops))/2;
yneg = zeros(1, length(x));
hold on 
errorbar(x, y, yneg, ypos, 'ko');

% >>> V2 colorfil and no horbar
x = 1:length(cellPops);
y = mean(geneListExpMat(pyraSCInds, cellPops));
yneg = zeros(1, length(x));
ypos = std(geneListExpMat(pyraSCInds, cellPops))./2;
h = figure('units', 'centimeters', 'position', [0,0, 30, 12])
scatter(x, y, 40, [70 173 250]/256, 'filled')

for i  = 1:length(cellPops)
    line([x(i) x(i)], [y(i), y(i)+ypos(i)], 'color', [70 173 250]/256)
end

x = x + .2;
y = mean(geneListExpMat(logical(scBasedMarkers.markerMat(:,6)), ...
                        cellPops));
ypos = std(geneListExpMat(logical(scBasedMarkers.markerMat(:,6)), ...
                        cellPops))/2;
yneg = zeros(1, length(x));
hold on 
scatter(x, y, 30, 'k')

for i  = 1:length(cellPops)
    line([x(i) x(i)], [y(i), y(i)+ypos(i)], 'color', 'k')
end


xticks([1:length(cellPops)])
xticklabels(xlabels)
xtickangle(90)
xlim([0, length(x)+1])
title(['average expression level of genes in the GTExBulk Excitatory ' ...
       'cells clusters'])
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sexpressionLevelOfGTExBulkPyramidal_SCpops_V02', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% ASTRO
% >>>> original
x = 1:length(cellPops);
y = mean(geneListExpMat(astroSCInds, cellPops));
yneg = zeros(1, length(x));
ypos = std(geneListExpMat(astroSCInds, cellPops))./2;
h = figure('units', 'centimeters', 'position', [0,0, 30, 12])
errorbar(x, y, yneg, ypos, 'o')

x = x + .2;
y = mean(geneListExpMat(logical(scBasedMarkers.markerMat(:,1)), ...
                        cellPops));
ypos = std(geneListExpMat(logical(scBasedMarkers.markerMat(:,1)), ...
                        cellPops))/2;
yneg = zeros(1, length(x));
hold on 
errorbar(x, y, yneg, ypos, 'ko');

% >>> V2 colorfil and no horbar
x = 1:length(cellPops);
y = mean(geneListExpMat(astroSCInds, cellPops));
yneg = zeros(1, length(x));
ypos = std(geneListExpMat(astroSCInds, cellPops))./2;
h = figure('units', 'centimeters', 'position', [0,0, 30, 12])
scatter(x, y, 40, [70 173 250]/256, 'filled')

for i  = 1:length(cellPops)
    line([x(i) x(i)], [y(i), y(i)+ypos(i)], 'color', [70 173 250]/256)
end

x = x + .2;
y = mean(geneListExpMat(logical(scBasedMarkers.markerMat(:,1)), ...
                        cellPops));
ypos = std(geneListExpMat(logical(scBasedMarkers.markerMat(:,1)), ...
                        cellPops))/2;
yneg = zeros(1, length(x));
hold on 
scatter(x, y, 30, 'k')

for i  = 1:length(cellPops)
    line([x(i) x(i)], [y(i), y(i)+ypos(i)], 'color', 'k')
end


% extra and saving - same for two versions

hold on 
xticks([1:length(cellPops)])
xticklabels(xlabels)
xtickangle(90)
xlim([0, length(cellPops)+1])
title(['average expression level of genes in the GTExBulk Astrocyte ' ...
       'cells clusters'])
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sexpressionLevelOfGTExBulkAstrocyte_SCpops_V02', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% oligo
x = 1:length(cellPops);
y =  mean(geneListExpMat(oligoSCInds, cellPops));
ypos = std(geneListExpMat(oligoSCInds, cellPops))/2;
yneg = zeros(1, length(x));
h = figure('units', 'centimeters', 'position', [0,0, 30, 12])
errorbar(x, y, yneg, ypos, 'o')

x = x + .2;
y = mean(geneListExpMat(logical(scBasedMarkers.markerMat(:,5)), ...
                        cellPops));
ypos = std(geneListExpMat(logical(scBasedMarkers.markerMat(:,5)), ...
                        cellPops))/2;
yneg = zeros(1, length(x));
hold on 
errorbar(x, y, yneg, ypos, 'ko');


% >>> V2 colorfil and no horbar
x = 1:length(cellPops);
y = mean(geneListExpMat(oligoSCInds, cellPops));
yneg = zeros(1, length(x));
ypos = std(geneListExpMat(oligoSCInds, cellPops))./2;
h = figure('units', 'centimeters', 'position', [0,0, 30, 12])
scatter(x, y, 40, [70 173 250]/256, 'filled')

for i  = 1:length(cellPops)
    line([x(i) x(i)], [y(i), y(i)+ypos(i)], 'color', [70 173 250]/256)
end

x = x + .2;
y = mean(geneListExpMat(logical(scBasedMarkers.markerMat(:,5)), ...
                        cellPops));
ypos = std(geneListExpMat(logical(scBasedMarkers.markerMat(:,5)), ...
                        cellPops))/2;
yneg = zeros(1, length(x));
hold on 
scatter(x, y, 30, 'k')

for i  = 1:length(cellPops)
    line([x(i) x(i)], [y(i), y(i)+ypos(i)], 'color', 'k')
end


xticks([1:length(cellPops)])
xticklabels(xlabels)
xtickangle(90)
xlim([0, length(x)+1])
ylim([0, 8.3])
title(['average expression level of genes in the GTExBulk Oligodendrocyte ' ...
       'cells clusters'])
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sexpressionLevelOfGTExBulkOligo_SCpops_V02', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% micro
x = 1:length(cellPops);
y =  mean(geneListExpMat(microSCInds, cellPops));
ypos = std(geneListExpMat(pyraSCInds, cellPops))/2;
yneg = zeros(1, length(x));

h = figure('units', 'centimeters', 'position', [0,0, 30, 12])
errorbar(x, y, yneg, ypos, 'o')

x = x + .2;
y = mean(geneListExpMat(logical(scBasedMarkers.markerMat(:,3)), ...
                        cellPops));
ypos = std(geneListExpMat(logical(scBasedMarkers.markerMat(:,3)), ...
                        cellPops))/2;
yneg = zeros(1, length(x));
hold on 
errorbar(x, y, yneg, ypos, 'ko');


% >>> V2 colorfil and no horbar
x = 1:length(cellPops);
y = mean(geneListExpMat(microSCInds, cellPops));
yneg = zeros(1, length(x));
ypos = std(geneListExpMat(microSCInds, cellPops))./2;
h = figure('units', 'centimeters', 'position', [0,0, 30, 12])
scatter(x, y, 40, [70 173 250]/256, 'filled')

for i  = 1:length(cellPops)
    line([x(i) x(i)], [y(i), y(i)+ypos(i)], 'color', [70 173 250]/256)
end

x = x + .2;
y = mean(geneListExpMat(logical(scBasedMarkers.markerMat(:,3)), ...
                        cellPops));
ypos = std(geneListExpMat(logical(scBasedMarkers.markerMat(:,3)), ...
                        cellPops))/2;
yneg = zeros(1, length(x));
hold on 
scatter(x, y, 30, 'k')

for i  = 1:length(cellPops)
    line([x(i) x(i)], [y(i), y(i)+ypos(i)], 'color', 'k')
end



xticks([1:length(cellPops)])
xticklabels(xlabels)
xtickangle(90)
xlim([0, length(x)+1])
ylim([-.4 6.4])
title(['average expression level of genes in the GTExBulk Microglia ' ...
       'cells clusters'])
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sexpressionLevelOfGTExBulkMicro_SCpops_V02', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% endo
x = 1:length(cellPops);
y =  mean(geneListExpMat(endoSCInds, cellPops));
ypos = std(geneListExpMat(endoSCInds, cellPops))/2;
yneg = zeros(1, length(x));
h = figure('units', 'centimeters', 'position', [0,0, 30, 12])
errorbar(x, y, yneg, ypos, 'o')

x = x + .2;
y = mean(geneListExpMat(logical(scBasedMarkers.markerMat(:,2)), ...
                        cellPops));
ypos = std(geneListExpMat(logical(scBasedMarkers.markerMat(:,2)), ...
                        cellPops))/2;
yneg = zeros(1, length(x));
hold on 
errorbar(x, y, yneg, ypos, 'ko');


% >>> V2 colorfil and no horbar
x = 1:length(cellPops);
y = mean(geneListExpMat(endoSCInds, cellPops));
yneg = zeros(1, length(x));
ypos = std(geneListExpMat(endoSCInds, cellPops))./2;
h = figure('units', 'centimeters', 'position', [0,0, 30, 12])
scatter(x, y, 40, [70 173 250]/256, 'filled')

for i  = 1:length(cellPops)
    line([x(i) x(i)], [y(i), y(i)+ypos(i)], 'color', [70 173 250]/256)
end

x = x + .2;
y = mean(geneListExpMat(logical(scBasedMarkers.markerMat(:,2)), ...
                        cellPops));
ypos = std(geneListExpMat(logical(scBasedMarkers.markerMat(:,2)), ...
                        cellPops))/2;
yneg = zeros(1, length(x));
hold on 
scatter(x, y, 30, 'k')

for i  = 1:length(cellPops)
    line([x(i) x(i)], [y(i), y(i)+ypos(i)], 'color', 'k')
end

xticks([1:length(cellPops)])
xticklabels(xlabels)
xtickangle(90)
xlim([0, length(x)+1])
ylim([-.4 6.4])
title(['average expression level of genes in the GTExBulk Endo ' ...
       'cells clusters'])
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sexpressionLevelOfGTExBulkEndo_SCpops', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

load('~/resultsAndFigures/secondProject/scBasedMarkers.mat')

% 10. likelihood ratio of the overlap of links in the networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('~/resultsAndFigures/secondProject/allNetworkOverlaps_noVarAdded_redo.mat')

% list of cluster IDs vs list of networks
h = figure('units', 'centimeters', 'position', [0,0, 12, 25])
theNames = clusterMeta.sortedClusterNames(1:69);
[a, b] = sort(theNames);
b = 1:69
heatmap({'GTExBulk', 'GTEx_CTC'}, [clusterMeta.sortedClusterNames(b)], ...
        (noResult.netLikeRatioOverlap(b, [84, 86])+1))
colormap(1 - bone)
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%snetworkOverlap_snucRNAseq_GTEx_redo', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 11. Table for count of genes in the clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% give me the cluster IDs 

myIDs = [239, 242, 243, 129, 252, 98, 65, 244, 231 53 35];
myTable = zeros(length(myIDs)+1, 8); % one extra for total

[a, b] = ismember(myIDs, GTExClusterRep_itself.clusterInds)

myTable(1:11, 1) = GTExClusterRep_itself.gc(b);
myTable(1:11, 2) = GTExClusterRep_tan.gc(b);
myTable(1:11, 3) = GTExClusterRep_tsn.gc(b);
myTable(1:11, 4) = GTExClusterRep_ctc.gc(b);
myTable(1:11, 5) = GTExClusterRep_Gliver.gc(b);
myTable(1:11, 6) = GTExClusterRep_tanLiver.gc(b);
myTable(1:11, 7) = GTExClusterRep_Gblood.gc(b);
myTable(1:11, 8) = GTExClusterRep_tanBlood.gc(b);

book = sum(myTable);
finalT = [myTable; book]

figureIDs = gtexCluster.clusterLabels(myIDs);

netLabels = {'GTExBulk', 'TAN_brain', 'TSN_brain', 'GTEx_CTC', ...
             'GTEx_liver', 'TAN_liver', 'GTEx_blood', 'TAN_blood'}

h = figure
heatmap(netLabels, [myIDs, 10], finalT)
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sgeneCountsInMarkerClusters_heatmap', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 12. examples for the last part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scExpGeneList = {'ANKRD36', 'CHD5'}

% apart from these examples, we have halva, finalIndsA and
% finalIndsB for all the likewise potential links

% geneListExpMat % this has the two genes
% sortLabs

% boxplot for expression in SC
groups = [1 1, 2, ones(1, 24)*3, ones(1, 45)*4, 5 6 7];

h = figure('units', 'centimeters', 'position', [0,0, 20, 12])
i = 8
g1 = filDataSet.geneSyms(finalIndsA(i))
g2 = filDataSet.geneSyms(finalIndsB(i))

subplot(1, 3, 1)
boxplot(geneListExpMat(finalIndsA(i), sortLabsInd), groups)
title(sprintf('%s - R2: %.2f', g1{1}, finalR1(i)))

xticks([1:7])
xticklabels({'Astrocyte', 'Endothelial', 'Excitatory', 'Inhibitory', ...
             'Microglia', 'OPC', 'Oligo'})
xtickangle(60)

subplot(1, 3, 2)
boxplot(geneListExpMat(finalIndsB(i), sortLabsInd), groups)
title(sprintf('%s - R2: %.2f', g2{1}, finalR2(i)))

xticks([1:7])
xticklabels({'Astrocyte', 'Endothelial', 'Excitatory', 'Inhibitory', ...
             'Microglia', 'OPC', 'Oligo'})
xtickangle(60)

% get teh box plot of expraession for CALM3 and NRGN
[a, b] = ismember(g1, dataSet.genes)
e1 = dataSet.mat(b,:);

[a, b] = ismember(g2, dataSet.genes)
e2 = dataSet.mat(b,:);

subplot(1, 3, 3)
boxplot(log2([e1' e2']))
title(sprintf('Expression in GTEx \n correlation Rank: %d/1000', halva(2, ...
                                                  i)))
xticks([1:2])
xticklabels({g1{1}, g2{1}})
xtickangle(60)

ylim([0, 18])
halva(:, i)

figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%sExamples_RoustSNUClinks_Plus15_%s_%s_%d', figFolder, ...
               g1{1}, g2{1}, i)
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


% 13. overlap between clusters: snuc network versus GTEx and CTC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 

% load gtexCluster
load('~/resultsAndFigures/secondProject/gtexClusters.mat')

% load GTExnet
load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat') 
gNet = GTExFiveNets.nets(2).net005;
gSyms = GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(2).expGenes);

% load CTCnet 
load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'correctedBinNets_logCorrected_jusResiduals_scBasedMarkers_redo' ...
      '.mat'])
cNet = ctc.net005;

inClusters = find(gtexCluster.memberCounts >= 20);
gOv = zeros(69, length(inClusters));
cOv = zeros(69, length(inClusters));
gO = zeros(69, length(inClusters));
cO = zeros(69, length(inClusters));

lcs = zeros(69, length(inClusters));
glcs = zeros(69, length(inClusters));
clcs = zeros(69, length(inClusters));
scGeneSyms = filDataSet.geneSyms;

for n = 1:10
    n
    load(sprintf(['~/networks/allenBrainSC/rpmAllFiveNets/' ...
                  'rpm_binNets_Exon_V4_allFive_net%d.mat'], n))
    snucNet = net.net005;
    snucSyms = scGeneSyms(net.expGenes);

    for i = 1:length(inClusters)
        thisCluster = inClusters(i);
        thisCGenes = gtexCluster.syms(gtexCluster.cs1547 == thisCluster);
        sNet = zeros(length(thisCGenes), length(thisCGenes));
        
        [a, b] = ismember(snucSyms, thisCGenes);
        sum(a)
        if (sum(a) > 10)
            sNet(b(a), b(a)) = snucNet(a, a);
            
            [ag, bg] = ismember(thisCGenes, gSyms);

            gOv(i, n) = sum(sum((sNet + gNet(bg, bg)) == 2))/ sum(sum(gNet(bg, ...
                                                              bg)));
            cOv(i, n) = sum(sum((sNet + cNet(bg, bg)) == 2))/ sum(sum(cNet(bg, ...
                                                              bg)));

            gO(i, n) = sum(sum((sNet + gNet(bg, bg)) == 2));
                                                            
            cO(i, n) = sum(sum((sNet + cNet(bg, bg)) == 2));

            lcs(i, n) = sum(sum(sNet))/2;
            glcs(i, n) = sum(sum(gNet(bg, bg)));
            clcs(i, n) = sum(sum(cNet(bg, bg)));
        end
    end
end

h = figure
heatmap(cOv'./gOv')

finalFil = ((lcs(:, 1) >=10) + (((gOv(:,1))+(cOv(:,1))) > 0)) == 2;
h = figure
heatmap(clusterMeta.sortedClusterNames(1:10), gtexCluster.clusterLabels(inClusters(finalFil)), ...
        cOv(finalFil, 1:10)./gOv(finalFil,1:10))

heatmap(1, gtexCluster.clusterLabels(inClusters(finalFil)), ...
        cOv(finalFil, 1)./gOv(finalFil,1))
colormap(gray)

figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = [figFolder 'resNetResult_snc01_ratioOfGTExBulkAndResidual_testNew']
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

h  = figure
heatmap({'gtexBulk', 'residual', clusterMeta.sortedClusterNames{1}}, gtexCluster.clusterLabels(inClusters(finalFil)), ...
        [glcs(finalFil,1), clcs(finalFil, 1), lcs(finalFil, 1)])
colormap(gray)
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = [figFolder 'resNetResult_actualLinkCountsInGTExAndCTCAndsnc01_redo']
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

h  = figure
heatmap({'gtexBulk', 'residual'}, ...
        gtexCluster.clusterLabels(inClusters(finalFil)), [gOv(finalFil, ...
                                                  1), cOv(finalFil, ...
                                                  1)])
colormap(gray)
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = [figFolder 'resNetResult_retrievedRatioInBulkAndResidual_snc01_redo']
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

h  = figure
heatmap({'gtexBulk', 'residual', clusterMeta.sortedClusterNames{1}}, gtexCluster.clusterLabels(inClusters(finalFil)), ...
        [gO(finalFil,1), cO(finalFil, 1), lcs(finalFil, 1)])
colormap(gray)
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = [figFolder 'resNetResult_actuaLinkCountsRetrieved_snc01_redo']
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% report this. 

%14. overlapping links  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this file is made in gettingTheTSSDistForBiPrimer.m
load(['~/resultsAndFigures/secondProject/' ...
      'overlappingLinksBetweenSCNets/highRepLinks.mat'])
overlapRes

[a, b] = ismember(overlapRes.gtexGeneSyms, gSyms);

fullCt = cNet(b(a), b(a)) + cNet(b(a), b(a))';
fullGt = overlapRes.gtexBinNet005(a, a);
inGenesC = overlapRes.clusters(a);
c1ind = inGenesC ==1;
c9ind = inGenesC == 9;
c10ind = inGenesC == 10;
c910ind = inGenesC > 8;

% I am getting the clusters of mitochondria, the unclusters genes
% and 9 and 10. >> 2 and three are gene families. PCHD@

[plotGenes, b] = ismember(inGenesC, [1 2 5 7 8 9 10]);
[plotGenes, b] = ismember(inGenesC, [9 10]);
[g9 , b] = ismember(inGenesC, [9]);
[g10 , b] = ismember(inGenesC, [10]);

h = figure
bar([sum(inGenesC == 9) sum(inGenesC == 10)])

figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
%file = [figFolder 'gtexNet_overlappingLinks_clusters910']
file = [figFolder 'gtexNet_overlappingLinks_clustersAll_noMito_colorBar']
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

h = figure
% heatmap(fullCt(c910ind, c910ind), 'GridVisible', 'off')
% colormap(1-gray/2)

heatmap(fullGt(plotGenes, plotGenes), 'GridVisible', 'off')
colormap(1-gray)

figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
%file = [figFolder 'gtexNet_overlappingLinks_clusters910']
file = [figFolder 'gtexNet_overlappingLinks_clustersAll_noMito']
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


h = figure
% heatmap(fullGt(c910ind, c910ind), 'GridVisible', 'off')
% colormap(1-gray/2)

heatmap(fullCt(plotGenes, plotGenes), 'GridVisible', 'off')
heatmap(fullCt, 'GridVisible', 'off')
colormap(1-gray)

figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
%file = [figFolder 'gtexNet_overlappingLinks_clusters910']
file = [figFolder 'gtexCTC_overlappingLinks_clustersAll_noMito_new_redo']
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% getting clusters 9 and 10 only
halva = sib > 8;
halva9 = sib == 9;
halva10 = sib == 10;

ct9 = fullCt(c9ind, c9ind);
ct10 = fullCt(c10ind, c10ind);
ct910 = fullCt(c9ind, c10ind);

ca1 = sum(sum(ct9))/(sum(c9ind)*(sum(c9ind)-1))
ca2 = sum(sum(ct10))/(sum(c10ind)*(sum(c10ind)-1))
ca3 = sum(sum(ct910))/(sum(c9ind)*(sum(c10ind)-1)/2)

ca1/ca3
ca2/ca3

gt9 = fullGt(c9ind, c9ind);
gt10 = fullGt(c10ind, c10ind);
gt910 = fullGt(c9ind, c10ind);

ga1 = sum(sum(gt9))/(sum(c9ind)*(sum(c9ind)-1))
ga2 = sum(sum(gt10))/(sum(c10ind)*(sum(c10ind)-1))
ga3 = sum(sum(gt910))/(sum(c9ind)*(sum(c10ind)-1)/2)

ga1/ga3
ga2/ga3


% >>> getting the overlap of the 9 ant 10
load(['~/resultsAndFigures/secondProject/' ...
      'overlappingLinksBetweenSCNets/highRepLinks.mat'])
overlapRes

ind9 = overlapRes.clusters == 9;
ind10 = overlapRes.clusters ==10;

d9 = sum(sum((overlapRes.sumNet(ind9, ind9) >= 10)));
d10 = sum(sum((overlapRes.sumNet(ind10, ind10) >= 10)));
d910 = sum(sum((overlapRes.sumNet(ind10, ind9) >=10 )));


ra1 = sum(sum(d9))/(sum(ind9)*(sum(ind9)-1))
ra2 = sum(sum(d10))/(sum(ind10)*(sum(ind10)-1))
ra3 = sum(sum(d910))/(sum(ind9)*(sum(ind10)-1)/2)

ra1/ra3
ra2/ra3

% 15. contamination in snuc-RNA seq data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% get the expression level for the first cell pop
myMat = thisExpMat; % from expressoinInSC.m

% >>>> get Ogan marker genes. 
load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])
markers

markerInds = [];
for i = [1 2 5 6 7 8 9]
    markerInds = [markerInds, markers.finalMarkerGenes{i}'];
end

markerGenes = markers.genes(markerInds);
scGeneSyms = filDataSet.geneSyms;
[a, b] = ismember(markerGenes, scGeneSyms);
markerExp = thisExpMat(b(a), :);

% >>> get it with scbased markers
load('~/resultsAndFigures/secondProject/scBasedMarkers.mat')
markerGenes = logical(scBasedMarkers.markerMat(:, 6));
markerExp = myMat(markerGenes, :);

% do the model and get the R2s 
[w, score, latent, tsquared, explained, mu] = pca(markerExp', ...
                                                  'NumComponents', 10);
T = markerExp' * w;
% correct it for all the genes
featureMat = score(:, 1:7);
regOut = zeros(size(thisExpMat));
coeffs = zeros(size(thisExpMat, 1), 8);
r2 = zeros(size(thisExpMat, 1), 1);
for i = 1:size(thisExpMat, 1)
    i
    y = thisExpMat(i, :);
    % mean(thisGeneExp)
    % var(thisGeneExp)

    lm = LinearModel.fit(featureMat, y, 'Intercept', true);
    coeff = lm.Coefficients.Estimate;
    coeffs(i, :) = coeff;
    
    r2(i) = lm.Rsquared.Adjusted;
    % pr = coeff(1) + coeff(2:end)' * featureMat';
    % %corrected = y - pr;
    % sib = [pr' y'];

    % e1 = sum((y - mean(y)).^2);
    % e2 = sum((y - pr).^2);

    % get the residual
    regOut(i, :) = (coeff(2:end)' * featureMat');
end

% give me the first cell pop network 
scGeneSyms = filDataSet.geneSyms;
load(sprintf(['~/networks/allenBrainSC/rpmAllFiveNets/' ...
              'rpm_binNets_Exon_V4_allFive_net%d.mat'], i))
% scNet = net.net005;
% thisNetSyms = scGeneSyms(net.expGenes);
firstNet = triu(net.net005, 1);
firstSyms = scGeneSyms(net.expGenes);

% get the small network
% >>> getting the overlap of the 9 ant 10
load(['~/resultsAndFigures/secondProject/' ...
      'overlappingLinksBetweenSCNets/highRepLinks.mat'])
overlapRes

[a, b] = ismember(overlapRes.geneSyms, firstSyms);
h = figure
fullPlot = firstNet(b(a), b(a)) + firstNet(b(a), b(a))';
heatmap(fullPlot,'GridVisible','off')

% build the network using the residual

% getting the distribution of the average R2 for marker enriched
% clusters and their neighbors 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('~/resultsAndFigures/secondProject/gtexClusters.mat')

myIDs = [231, 244, 57, 45, 107, 35, 53, 92, 161,242 239, 243, ...
         87, 136, 232, 174, 129, 98, 65, 252,250, 217, 233, 164]
myIDs = [242, 243, 239, 129, 252, 98,65, 231, 244, 53 35]

myIDs = [253, 189, 197, 163, 159, 113, 40, 88, 234, 248, 130, 83, ...
         78, 167, 241, 190, 147, 233, 164, 105, 52, 170, 168, 124, ...
         169, 136, 250, 107]
inCs = gtexCluster.memberCounts >= 20;

kado = find(inCs);

[a, b] = ismember(kado, myIDs);

mbds = [myIDs, kado(~a)]
groupings = [ones(1,28), ones(1, 41)*2]

boxplot(gtexCluster.meanRs(mbds), groupings)

quantile(gtexCluster.meanRs(mbds), 3)
median(gtexCluster.meanRs(mbds))

sum(gtexCluster.meanRs(myIDs) < 0.2753)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 1:69
    i
    load(sprintf(['~/networks/allenBrainSC/rpmAllFiveNets/' ...
                  'rpm_binNets_Exon_V4_allFive_net%d.mat'], n))
    % scNet = net.net005;
    % thisNetSyms = scGeneSyms(net.expGenes);
    sumInNetGenes(i) = sum(sum(net.net005)>0);
end

% 17. removing the mito from the heatmap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['~/resultsAndFigures/secondProject/' ...
      'overlappingLinksBetweenSCNets/highRepLinks.mat'])

ins = overlapRes.sumNet(overlapRes.clusters > 8, overlapRes.clusters>8);
h = figure
% heatmap(fullGt(c910ind, c910ind), 'GridVisible', 'off')
% colormap(1-gray/2)

heatmap(log2(ins), 'GridVisible', 'off')
colormap(1-gray)

figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
%file = [figFolder 'gtexNet_overlappingLinks_clusters910']
file = [figFolder 'sumNetLog2_overlappingLinks_noMito']
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 18. getting rank of CALM3 and NRGN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the GTEx corr ranks 
load('~/data/GTEx/Brain_Cortex_expGenes.mat') 
sib = corr(log2(dataSet.mat+1)');
upCorr = sib(logical(triu(ones(size(sib)), 1)));
qs = quantile(upCorr, 1000);

[a, b] = ismember({'NRGN', 'CALM3'}, dataSet.genes)

sib(b,b)

% 19. enrichment of functional terms in GTExBulk vs CTC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['~/resultsAndFigures/secondProject/moduleRepresentation/' ...
      'moduleRep_GTExBrainCortex_CTC_logCorrected.mat'])
gtexBinNetCTC

load(['~/resultsAndFigures/secondProject/moduleRepresentation/' ...
      'moduleRep_GTExBrainCortex.mat'])
gtexBinNet

load(['~/resultsAndFigures/secondProject/moduleRepresentation/' ...
      'funcRes.mat'])
funcRes % tan tsn GTExBin GTExCTC GTExBlood GTExLiver

% how many of the enriched terms have 2 or more links - all of
% them,

gtexInTerms = (funcRes.uniqueTermIDs(logical(funcRes.funPresence(3, ...
                                                  :))));

[a, b] = ismember(gtexBinNet.inTermsGOID, gtexInTerms);

gtexIDs = gtexBinNet.inTermsGOID(a);
gtexPs = gtexBinNet.ps(a);

ctcInTerms = (funcRes.uniqueTermIDs(logical(funcRes.funPresence(4, ...
                                                  :))));

[a, b] = ismember(gtexBinNetCTC.inTermsGOppID, ctcInTerms);

nctcIDs = gtexBinNetCTC.inTermsGOID(a);
ctcPs = gtexBinNetCTC.ps(a);

uniqueIDs = unique([ctcIDs, gtexIDs]);

gfps = ones(1, length(uniqueIDs))*-1;
cfps = ones(1, length(uniqueIDs))*-1;

[a, b] = ismember(ctcIDs, uniqueIDs);
cfps(b) = ctcPs;

[a, b] = ismember(gtexIDs, uniqueIDs);
gfps(b) = gtexPs;

[a, b] = ismember(uniqueIDs, funcRes.uniqueTermIDs);
uniqueTerms = funcRes.uniqueTermNames(b);

%%% writing the file
myFile = ['~/resultsAndFigures/secondProject/suppFiles/' ...
         'GTExBulk_CTC_functionsEnriched.txt'];
fid = fopen(myFile, 'w')
fprintf(fid, ['GOID\tGOTerm\tGTExBulk_p\tGTEx_residual_p\n'])

% writing the file
for i = 1:length(uniqueTerms)
    fprintf(fid, '%s\t%d\t%e\t%e\n', uniqueTerms{i}, uniqueIDs(i), gfps(i), ...
            cfps(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'GTExitselfClusterRep.mat'])

load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'GTExCTCClusterRep.mat'])
GTExClusterRep_ctc

ctc
load('~/resultsAndFigures/secondProject/gtexClusters.mat')
gtexCluster

book = zeros(69, 3);
for i = 1:69
    clusterID = GTExClusterRep_ctc.clusterInds(i);
    cgenes = gtexCluster.syms(gtexCluster.cs1547 == clusterID);
    [a, b] = ismember(cgenes, ctc.geneSyms);

    smallMat = ctc.net005(b, b);
    book(i, 1) = sum(sum(smallMat));
    book(i, 2) = GTExClusterRep_ctc.lc(i);
    book(i, 3) = GTExClusterRep_itself.lc(i);
end

kado = [GTExClusterRep_ctc.clusterInds', book];

corr(book)

% 18. just plot the R2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('~/data/GTEx/Brain_Cortex_expGenes.mat') 

load('~/resultsAndFigures/secondProject/gtexClusters.mat')

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_AllenSCBasedMarker.mat'])

[a, b] = ismember(dataSet.genes, gtexCluster.syms);
sum(a)

sib = result.r2(a);

% 19. writhe the functional enrichment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% got them from moduleRepresentation.m
load('~/resultsAndFigures/secondProject/moduleRepresentation/moduleRep_GTExBrainCortex.mat')

load(['~/resultsAndFigures/secondProject/' ...
      'moduleRep_GTExBrainCortex_CTC_redo.mat'])

% get the whole fterms 

% file will have six columns, one for ctc for for gtex bulk, ids
% terms , ps and pass

bothps = [gtexBinNet.ps' gtexBinNetCTC.ps'];

bothpsAd = bothps * length(bothps);
passps = bothpsAd <= .1;
halva = (sum(passps') > 0);

terms = gtexBinNet.inTermNames(halva);
ids = gtexBinNet.inTermsGOID(halva);
finalps = bothpsAd(halva, :);
 
myFile = ['~/resultsAndFigures/secondProject/suppFiles/redo/' ...
         'GTExBulk_CTC_functionsEnriched_redo.txt'];
fid = fopen(myFile, 'w')
fprintf(fid, ['GOID\tGOTerm\tGTExBulk_p\tGTEx_residual_p\n'])

% writing the file
for i = 1:length(terms)
    fprintf(fid, '%s\t%d\t%e\t%e\n', terms{i}, ids(i), finalps(i,1), ...
            finalps(i,2));
end

% 20. snuc networks and their link counts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gcs = zeros(69, 1);
lcs = zeros(69, 1);
netGenes = zeros(69, 1);
for n = 1:69
    n
    load(sprintf(['~/networks/allenBrainSC/rpmAllFiveNets/' ...
                  'rpm_binNets_Exon_V4_allFive_net%d.mat'], n))
    gcs(n) = sum(net.expGenes);
    lcs(n) = sum(sum(net.net005));
    netGenes(n) = sum(sum(net.net005) > 0);
end

h = figure('units', 'centimeters', 'position', [0,0, 10, 25])
heatmap(1, clusterMeta.sortedClusterNames(1:69), netGenes , 'GridVisible', 'off')
colormap(1-gray)

figFolder = ['~/resultsAndFigures/secondProject/']
file = [figFolder 'scNets_geneCounts_actual']
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

h = figure('units', 'centimeters', 'position', [0,0, 10, 25])
heatmap(1, clusterMeta.sortedClusterNames(1:69), lcs , 'GridVisible', 'off')
colormap(1-gray)

figFolder = ['~/resultsAndFigures/secondProject/']
%file = [figFolder 'gtexNet_overlappingLinks_clusters910']
file = [figFolder 'scNets_linkCounts']
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 21. Review: correlation of genes in snuc one dataSet one cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load GTExCluster 

load('~/data/GTEx/Brain_Cortex_expGenes.mat') 

load('~/resultsAndFigures/secondProject/gtexClusters.mat')

scExpGeneList;
thisExpMat;

% get the genes in cluster 242 
temp = gtexCluster.syms(gtexCluster.cs1547 == 242);
[as, bs] = ismember(temp, scExpGeneList);
myGenes = temp(as);

% get their correlaion and rank with each other 
[a, b] = ismember(myGenes ,dataSet.genes);
sib = corr(dataSet.mat(b(a), :)');

upC = sib(logical(triu(ones(size(sib)), 1)));

% in single cell, get the correlation and ranks 
[a, b] = ismember(myGenes, scExpGeneList);
smallMat = thisExpMat(b, :);

% % get the corr with zero values 
smallMatNAN = smallMat;
smallMatNAN(smallMatNAN == 0) = nan;
tic
[sib1, p1] = corr(smallMatNAN', 'rows', ...
                  'pairwise', 'Tail', 'right');
toc

% Bonferroni correction - looking at n >= 20
kado = (smallMat > 0) + 0;
book = kado * kado';
fdrFac = sum(sum(triu(book, 1) >= 20));
fdrFacRatio(i) = fdrFac / (expGeneCounts(i) * (expGeneCounts(i) - 1)/2);

BonBinNet = triu(p1 < (.1/fdrFac), 1);

% BH correction
throwOff = book < 20;
p1(throwOff) = 1;

upP1 = p1(logical(triu(ones(size(p1)), 1)));
upSib1 = sib1(logical(triu(ones(size(sib1)), 1)));

[sa, sb] = sort(upC, 'descend');

corr(sa(1:10000), upSib1(sb(1:10000)))
scps = upP1(sb(1:10000));
sccs = upSib1(sb(1:10000));

h = figure
plot(sa(1:10000), sccs, '.')
hist(scps, 30)

inds = datasample(1:length(upP1), 10000);

[a, b] = sort(upP1, 'ascend');
qs = quantile(a(1:fdrFac), 1000);

% 22. Review: correlation of PCs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('~/data/GTEx/Brain_Cortex_expGenes.mat')
dataSet.mat = log2(dataSet.mat+1);

% load markers 
load('~/resultsAndFigures/secondProject/scBasedMarkers.mat')

% for individual markers
for i = [1 2 3 5 6]
    i
    myMarkers = ...
    scBasedMarkers.genes(logical(scBasedMarkers.markerMat(:,i)));
    [a, b] = ismember(myMarkers, dataSet.genes);
    markerExp = dataSet.mat(b(a), :);
    [w, score, ~, ~, exp, ~] = pca(markerExp', 'NumComponents', 10);
    scores{i} = score(:, 1:3);
    exps{i} = exp; % variation explained
end

kado1 = corr(myFirstPCs);
kado1 = kado1 - eye(size(kado1));
kado(kado < 0) = 0;

kado2 = corr(mySecondPCs);
kado2 = kado2 - eye(size(kado2));

% now for all the markers together
geneList = scBasedMarkers.sortedGenes;

[a, b] = ismember(dataSet.genes, geneList);
sum(a)
markerExp = dataSet.mat(a, :);
[w, score, latent, tsquared, explained, mu] = pca(markerExp', ...
                                                  'NumComponents', 10);

% putting together the matrix
wholeMat = [scores{1} scores{2} scores{3} scores{5} scores{6} ...
            score];
sib = corr(wholeMat);
sib = sib - eye(25);
h = figure
heatmap(sib, 'GridVisible', 'off')
colormap(1- [1-gray' copper']')

figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sPCcorrs_draft_GR_V3', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

expMat = [exps{1}(1:10), exps{2}(1:10) exps{3}(1:10) exps{5}(1:10) ...
          exps{6}(1:10), explained(1:10)];
h = figure('units', 'centimeters', 'position', [0,0, 20, 7])
heatmap(expMat)

figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sExplained_cellTypesAndAll', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

%%%%% pattern of expression with PCs

% for individual markers
meanMarkerExps = zeros(5, 114);
order = [1 2 3 5 6];
for i = 1:5
    myMarkers = ...
    scBasedMarkers.genes(logical(scBasedMarkers.markerMat(:,order(i))));
    [a, b] = ismember(myMarkers, dataSet.genes);
    markerExp = zscore(dataSet.mat(b(a), :)');
    meanMarkerExps(i, :) = mean(markerExp');
end

% now for all the markers together get PCs
geneList = scBasedMarkers.sortedGenes;

[a, b] = ismember(dataSet.genes, geneList);
sum(a)
markerExp = dataSet.mat(a, :);
[w, score, latent, tsquared, explained, mu] = pca(markerExp', ...
                                                  'NumComponents', 10);

% putting together the matrix
myMat = [meanMarkerExps', score(:, 1:7)];
sib = corr(myMat);
sib = sib - eye(size(sib));
sib = sib + (eye(size(sib)).* min(sib(:)));
h = figure
heatmap(sib, 'GridVisible', 'off')
colormap(1- [1-gray' copper']')

figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sPCcorrs_markerMeanExps_7pcs_diagone', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

expMat = [exps{1}(1:10), exps{2}(1:10) exps{3}(1:10) exps{5}(1:10) ...
          exps{6}(1:10), explained(1:10)];
h = figure('units', 'centimeters', 'position', [0,0, 20, 7])
heatmap(expMat)

figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sExplained_cellTypesAndAll', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 23. Review: the marker counts - comment 31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% markers from Mancarci et al
load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])
IDs = [1 2 5 8 10];
markers.finalMarkerGenes
counts = [92, 24, 49, 71, 24]
h = figure
bar(counts)
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%scountOfMarkerGenes_GR_comment31', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 24. Review: count of genes in marker enriched cluster group
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('~/resultsAndFigures/secondProject/gtexClusters.mat')

astroIDs = [231, 244, 57, 45, 107, 35, 53, 84, 161]
[aas, b] = ismember(gtexCluster.cs1547, astroIDs);
sum(aas)

pyraIDs = [242 239, 243, ...
         87, 136, 232, 174, 129]
[apy, b] = ismember(gtexCluster.cs1547, pyraIDs);
sum(apy)

oligoIDs = [252]
[aol, b] = ismember(gtexCluster.cs1547, oligoIDs);
sum(aol)

endoIDs = [98]
[aen, b] = ismember(gtexCluster.cs1547, endoIDs);
sum(aen)

microIDs = [65]
[ami, b] = ismember(gtexCluster.cs1547, microIDs);
sum(ami)

sum(ami + aen + aol + apy + aas)

% 25. Review: if there are marker genes in the highRep clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get Ogan marker genes. 
load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])

load(['~/resultsAndFigures/secondProject/' ...
      'overlappingLinksBetweenSCNets/highRepLinks.mat'])

load('~/resultsAndFigures/secondProject/scBasedMarkers.mat')

inMs = sum(scBasedMarkers.markerMat(:, [1:3 5 6])');

inMGenes = scBasedMarkers.genes(inMs > 0);

inMGenes = scBasedMarkers.genes(logical(scBasedMarkers.markerMat(:, 6)));

[a, b]= ismember(inMGenes, overlapRes.geneSyms);
myInds = b(a);
myInds(11) = []; % to keep them within the two clusters 9 and 10
smyInds = sort(myInds)

smallNet = full(overlapRes.SCSumBinaryNet(smyInds, smyInds));
[a, b, c] = find(smallNet);

% distribution of R2 for genes in those clusters
load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_AllenSCBasedMarker_redo.mat'])

load('~/data/GTEx/Brain_Cortex_expGenes.mat')
dataSet.mat = log2(dataSet.mat+1);

[a, b] = ismember(overlapRes.geneSyms(overlapRes.clusters == 9), ...
                  dataSet.genes);

srs = result.r2(b(a));
mean(srs)

[a, b] = ismember(gtexCluster.syms(gtexCluster.cs1547 == 242), dataSet.genes);
mrs = result.r2(b(a));
mean(mrs)

% 26. Review: 5 versus 7 pcs and dists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distribution of R2 for genes in those clusters

R7 = load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_AllenSCBasedMarker_redo.mat'])


% R5 = load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
%       'regressionResult_5PCA_AllenSCBasedMarker.mat'])

load('~/resultsAndFigures/secondProject/gtexClusters.mat')

[a, b] = ismember(dataSet.genes, gtexCluster.syms);

myvals = R7.result.r2(a);

[f, xi] = ksdensity(myvals);
h = figure;
plot(xi, f)
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sbrainR2dist_ksdensity', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

h = figure;
histogram(myvals, 100, 'DisplayStyle', 'stairs')
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sbrainR2dist_histStairs', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 27. get the correlations of R2, bulkVar, CTVar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% >>>>> get the CT
% loading data - exon
load('~/data/brainSingleCell/filDataSet_exon_V4.mat')
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exon_V4.mat'])
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exon_V4_clusterLabels.mat'])

expMat = filDataSet.expMat(:, clusterMeta.inCells);
% normalizing for the million count
sumExp = sum(expMat);
milFac = sumExp ./ 1000000;
normExp = zeros(size(expMat));
sampleCount = size(expMat, 2);
for i = 1:sampleCount
    normExp(:, i) = expMat(:, i)./milFac(i);
end

scExpGeneList = filDataSet.geneSyms;

% from expressionInSC.m
expMat = normExp;
geneListExpMat = zeros(length(scExpGeneList), 69); 
for i = 1:75 % for each cluster - here I get the count of sc types
             % as 69
    i
    net(i).clusterName = clusterMeta.sortedClusterNames(i);
    [a, b] = ismember(clusterMeta.clusters, ...
                      clusterMeta.sortedClusterNames(i));
    %    thisExpMat = log2(expMat(:, a)+1);
     thisExpMat = log2(expMat(:, a)+1);    

    geneListExpMat(:, i) = mean(thisExpMat');
end

myMeans = mean(geneListExpMat');
scsVar = var(geneListExpMat');

% >>>> get the bulk data Var
load('~/data/GTEx/Brain_Cortex_expGenes.mat')
dataSet.mat = log2(dataSet.mat+1);

load('~/resultsAndFigures/secondProject/gtexClusters.mat')

myGenest = gtexCluster.syms;

% >>>> get the R2
R7 = load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_AllenSCBasedMarker_redo.mat'])

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_AllenSCBasedMarker.mat'])

[a, b] = ismember(filDataSet.geneSyms, dataSet.genes);
sib = var(geneListExpMat');
scVar = sib(a);
bulkVar = var(dataSet.mat(b(a), :)');

r2 = result.r2(b(a)); % with 5pca 
corr(r2, scVar') % 
corr(bulkVar', scVar')

r2 = R7.result.r2(b(a)); % with 7pca
corr(r2, scVar') % 
corr(bulkVar', scVar')


