% here I want to compare the dense clusters in the GTEx bulk, GTEx,
% TAN and TSN
% corrected and single cells
% 0. test: SC clustering 
% 1. building gtex net and comparing the clusters (the marker ones)
% 1.5. testing random stuff
% 2. testing presence of markers
% 3. comparing the GTEx cluster with TSN and TAN, report the R2 for
% each of the GTEx clusters (load this)
% 3.5. comparing the GTEx cluster with Simulated networks 
% 4.RESULT: The base clustering plot: I have the average r2 for each GTEx
% cluster and the reproducibility of the clusters in SC data  
% 5. Enrichment of functional terms in each cluster (GTEx)
% 5.5. HKG in clusters 
% 6. plotting the distribution of r2 in different clusters
% 7. 
% 8. Identification of the CT clusters based on the SC markers and
% Ogan Markers - bar plots of LC go here
% 9. Density of the clusters in GTEx
% 10. The big fat plot of link proportion 
% 11. The big fat bar plot 
% 12. comparing the count of links in GTEx and ctc
% 13. building the parts of the big plot
% $. random testing ... 
% 14. getting the joint function of the genes in some clusters

% 0. test: SC clustering 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% let's try and cluster the single cells with topological overlap  
load('~/networks/allenBrainSC/binNets_Exon_V4_allFive_net7.mat') 

scNet7 = net.net01;
scSyms7 = filDataSet.geneSyms(net.expGenes);

scNet = net.net01;
scSyms = filDataSet.geneSyms(net.expGenes);

sum(sum(net.net01))
myNet = net.net01;
mySyms = filDataSet.geneSyms(net.expGenes);
nds = sum(myNet);
sum(nds >= 2)  % I think I can do that, hierarchical clustering
               % using that 

% get the topological overlap - 
smallNet = myNet(nds >=2, nds >=2) + 0;
smallSyms = mySyms(nds >=2);
smallnds = nds(nds >=2);
tic
overlap = smallNet * smallNet;
toc

book = diag(overlap);
removed = find(book == 0);
overlap(removed, :) = [];
overlap(:, removed) = [];
smallSyms(removed) = [];
smallNetRM = smallNet;
smallNetRM(:, removed) = [];
smallNetRM(removed, :) = [];
smallNDs = smallnds;
smallNDs(removed) = [];

% get the normalization factor (since I know that the minimum
% degree is 2, I don't have to plus one it)
normMat = zeros(size(overlap));
for i = 1:size(overlap, 1);
    i
    for j = 1:size(overlap, 1);
        normMat(i, j) = min([overlap(i,i), overlap(j, j)]);
        normMat(j, i) = normMat(i, j);
    end
end

top = (overlap + smallNet) ./ (normMat - smallNet + 1);

myDist = 1 - top;
myVDist = squareform(tril(myDist, -1));

z = linkage(myVDist, 'average');
dendrogram(z, 2000)

t = cluster(z, 'cutoff', .3, 'criterion', 'distance');
t = cluster(z, 'cutoff', 1.1547);
max(t)
book = hist(t, max(t)); % count of members in each cluster
sum(book > 40)
sum(book >= 20)
sum(book(book >= 5)) % 3636 genes belong to clusters with greater
                    % than five members and 446 genes are not "that
                    % close" to anything (.3 is the threshold, they
                    % are single clusters: but what does a .3
                    % distance mean in this network? what is the
                    % density of the clusters I am getting, and how
                    % far away are they from other clusters? 
% 
scCluster.cs = t;
scCluster.sysm = smallSyms;
save('~/resultsAndFigures/secondProject/scClusters_net1_01.mat', ...
     'scCluster')

plot(sort(book(book >= 5)))

myClusterIDs = find(book >= 5);
length(myClusterIDs)
lr = zeros(1, length(myClusterIDs));

for i = 1:length(myClusterIDs)
    myClusterGenes = t == myClusterIDs(i);
    sum(myClusterGenes)
    myClusterMat = smallNetRM(myClusterGenes, myClusterGenes);
    lr(i) = (sum(sum(myClusterMat)) / (sum(myClusterGenes)*(sum(myClusterGenes)-1)/2)) ...
        /.005;
end

hist(lr, 20)
% a strage observation is that some clusters have very low count of
% links between them, while ~60% of them have more than 10 times
% than expected by chance

% where are the marker genes, where are the ribosome genes and etc
% etc

% marker genes do cluster in GTEx, how about in this one.

% what are the clusters enriched in (finding the relevant
% clusters...)

% do the whole f enrichment... 

% and then the conclusion is that which clusters, which are
% identified in bulk, are also identified in SC.  

% 1. building gtex net and comparing the clusters (the marker ones)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat')

gtexGeneSyms = GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(2).expGenes);
gtexNet = GTExFiveNets.nets(2).net005;
gtexND = sum(gtexNet + gtexNet');

gtexSNet = gtexNet(gtexND > 0, gtexND > 0);
gtexsfNet = gtexSNet + gtexSNet';
gtexFSyms = gtexGeneSyms(gtexND >0);

tic
gtexOverlap = gtexsfNet * gtexsfNet;
toc

gtexNormMatUp = zeros(size(gtexOverlap));
maxOs = diag(gtexOverlap);
tic
for i = 1:size(gtexOverlap, 1);
    i
    for j = i:size(gtexOverlap, 1);
        a = maxOs(i);
        b = maxOs(j);
        tempMin = min(a, b);
        gtexNormMatUp(i, j) = tempMin;
        % gtexNormMat(j, i) = tempMin;
    end
end
toc
gtexNormMat = triu(gtexNormMatUp, 1) + gtexNormMatUp' + 1 - gtexsfNet;

kado = diag(gtexOverlap);
gtexRemove = find(kado == 0);
gtexTop = gtexOverlap ./ gtexNormMat;

gtexDist = 1 - gtexTop;
gtexDist(:, gtexRemove)= [];
gtexDist(gtexRemove, :)= [];
gtexFSyms(gtexRemove) = [];
gtexVDist = squareform(tril(gtexDist, -1));
 
z = linkage(gtexVDist, 'average');
dendrogram(z, 2000)

t = cluster(z, 'cutoff', 1.1546); % go for 1.1547 to get the 300 clusters
%t = cluster(z, 'cutoff', .3, 'criterion', 'distance');
max(t)
h = figure
hist(t, max(t))
book = hist(t, max(t)); % count of members in each cluster

myGeneSyms = gtexFSyms(t == 115);

[a, b] = ismember(myGeneSyms, gtexGeneSyms);
myrs=result.r2(b);
h = figure
hist(myrs)

[a, b] = ismember(veryCorrGenesSyms, gtexCluster.syms);

[a, b] = ismember(veryCorrGenesSyms, gtexGeneSyms);

sum(result.r2(b) > .5)

% I think these two big clusters in GTEx are the ones that I had
% found before: dictated by the ctp

gtexCluster.cs1547 = t;
gtexCluster.cs1546 = t;
gtexCluster.syms = gtexFSyms;
clusterLabels = cell(1, max(gtexCluster.cs1547));
clusterMemberCount = zeros(1, max(gtexCluster.cs1547));
for i = 1:max(gtexCluster.cs1547)
    clusterMemberCount(i) = sum(gtexCluster.cs1547 == i);
    clusterLabels{i} = sprintf('ID%d - %d', i, clusterMemberCount(i));
end
gtexCluster.clusterLabels = clusterLabels;
gtexCluster.memberCounts = clusterMemberCount;
save('~/resultsAndFigures/secondProject/gtexClusters.mat', ...
     'gtexCluster')

gtexCluster.cs1546 = t;
gtexCluster.syms = gtexFSyms;
save('~/resultsAndFigures/secondProject/gtexClusters_14102Genes.mat', ...
     'gtexCluster')

load('~/resultsAndFigures/secondProject/gtexClusters_14102Genes.mat')

load('~/resultsAndFigures/secondProject/gtexClusters.mat')
% 1.5. testing random stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('~/resultsAndFigures/secondProject/gtexClusters.mat')
load('~/resultsAndFigures/secondProject/scClusters.mat')

% getting the count of genes in each cluster
gtexcCounts = hist(gtexCluster.cs1547, max(gtexCluster.cs1547));
sccCounts = hist(scCluster.cs, max(scCluster.cs));

% getting the index of clusters with more than 5 genes in sc

sccInds = find(sccCounts >= 20);

%gtexCInds = find(gtexcCounts >= 20);

kado = zeros(1, length(sccInds));
memberCounts = zeros(1, length(sccInds));
cellType = zeros(1, length(sccInds));
belonging = zeros(length(sccInds), max(gtexCluster.cs));
threeMaxValues = zeros(length(sccInds), 3);
threeMaxIDs = zeros(length(sccInds), 3);
for i = 1:length(sccInds)
    thisInds = scCluster.cs == sccInds(i);
    sum((thisInds))
    thisSyms = scCluster.sysm(thisInds);

    [a, b] = ismember(thisSyms, gtexCluster.syms);
    gClusters = (gtexCluster.cs(b(a)));
    
    % d = length(unique(gtexCluster.cs(b(a))));
    % kado(i) = d/sum(a);
    
    memberCounts(i) = sum(a);
    
    sib = hist(gClusters, unique(gClusters));
    values = sib ./ sum(a);
    belonging(i, unique(gClusters)) = values;    
    
    thisGTExIDs = unique(gClusters);
    
    [sa, sb] = sort(values, 'descend')
    if length(sa) > 2
        threeMaxValues(i, :) = sa(1:3);
        threeMaxIDs(i, :) = thisGTExIDs(sb(1:3));
    else
        threeMaxValues(i, 1:length(sa)) = sa;
        threeMaxIDs(i, 1:length(sa)) = thisGTExIDs(sb);
    end
end

temp = sum(belonging);
belongingRed = belonging(:, temp >0);
IDs = [1:314];
IDsRed = IDs(temp>0);

heatmap(threeMaxValues')

% 2. testing presence of markers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do marker genes fall in the same cluster? 
% loading the GTEx cluster
load('~/data/GTEx/Brain_Cortex_expGenes.mat')
load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])

i = 1
[a, b] = ismember(markers.genes(markers.finalMarkerGenes{i}), ...
                  scCluster.sysm);
sum(a)/length(a)
sum(a)
h = figure
hist(scCluster.cs(b(a)), unique(scCluster.cs))

i = i + 1
% >> are marker genein same GTEx cluster?
i = 10
[a, b] = ismember(markers.genes(markers.finalMarkerGenes{i}), ...
                  gtexCluster.syms);
sum(a)/length(a)
sum(a)
%h = figure
book = hist(gtexCluster.cs1546(b(a)), ...
            unique(gtexCluster.cs1546(b(a))))
% hist(gtexCluster.cs1547(b(a)), ...
%             unique(gtexCluster.cs1547(b(a))))

kado = unique(gtexCluster.cs1547(b(a)));
[book;kado']

length(unique(gtexCluster.cs1547(b(a))))
i = i + 1

% Checking the presence of the marker genes
myInds = find(highCount); % clusters with count > 20
majorCID = zeros(1, 10);
for i = 1:10
    myGenes = dataSet.genes(markers.finalMarkerGenes{i});
    [a, b] = ismember(myGenes, gtexCluster.syms);
    sum(a)
    sib = unique(cs(b(a)));
    book = hist(cs(b(a)), sib);
    [a, b] = max(book);
    majorCID(i) = sib(b);
end

% I assume that most of the clusters in GTEx will "contain"
% clusters in the SC, so let's find the belongingness of the SC to
% the GTEx

% for each cluster in SC, find the gene symbols, 

% for each cluster in GTEx, give me the percentage of SC which is
% in it: I suspect that some clusters in SC will belong to the same
% GTEx cluster ... ALSO, marker genes in GTEx should belong to the
% same cluster ...

% So, that's that. Also, what happens to the extra dense clusters
% in GTEx? was it for the more dense network? cause in that network
% I think I had nds between only some enes. 

% 3. comparing the GTEx cluster with TSN and TAN, report the R2 for
% each of the GTEx clusters:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('~/data/GTEx/Brain_Cortex_expGenes.mat')

load('~/resultsAndFigures/secondProject/gtexClusters.mat')

% load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
%       'regressionResult_5PCA_10CellTypes_MarkersV2_cthr7.mat'])

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_AllenSCBasedMarker.mat'])

% for each cluster with N > 10, give me R2 
cs = gtexCluster.cs1547;
cgCount = hist(cs, 1:max(cs));
meanR2 = zeros(1, max(cs));
memc = zeros(1, max(cs));
for i = 1:max(cs)
    if(length(find(cs == i)) > 1000)
        i
    end
    
    if (length(find(cs == i)) >= 10)
        mySyms = gtexCluster.syms(cs == i);
        [a, b] = ismember(mySyms, dataSet.genes);
        meanR2(i) = mean(result.r2(b(a)));
        memc(i) = sum(a);
    end
end
 
[a, b] = sort(meanR2(meanR2>0));
h = figure
plot(a, '.')

book = memc(meanR2>0);
h = figure
plot((book(b)), '.')

% getting the density of networks in affy and SC 
load('~/data/general/GPL570GemmaMapNEW.mat')

load('~/data/general/tissueExpGenes/brainExpGenes0.8.mat')
affyBrainExpG = expGenesInd;

load('~/data/general/tissueExpGenes/bloodExpGenes0.8.mat')
affyBloodExpG = expGenesInd;

load('~/data/general/tissueExpGenes/liverExpGenes0.8.mat')
affyLiverExpG = expGenesInd;

load(['~/networks/tissues/blood/binaryNet_FDR5e-' ...
      '5_0.8Expr_Ind0.10.mat'])

load(['~/networks/tissues/brain/binaryNet_FDR5e-' ...
      '5_0.8Expr_Ind0.10.mat'])
tan = binNet;
clear binNet;

load(['~/networks/tissues/blood/binaryNet_FDR5e-' ...
      '5_0.8Expr_Ind0.10.mat'])
bloodTAN = binNet;

load(['~/networks/tissues/liver/binaryNet_FDR5e-' ...
      '5_0.8Expr_Ind0.10.mat'])
liverTAN = binNet;

FDR = '0012'
FC = '3'
load(sprintf(['~/resultsAndFigures/firstProject/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))
tsn = finalTable(2).wholeNet;

% clusters with R2 > .4
highR2 = meanR2 > .4;
highCount = memc > 20;
both = (highR2 + highCount) == 2;

% examine the rep of clusters, enrichment of functional terms and
% presence of marker genes
load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat')

gtexGeneSyms = GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(2).expGenes);
gtexNet = GTExFiveNets.nets(2).net005;
gtexND = sum(gtexNet + gtexNet');

[a, b] = ismember(gtexCluster.syms, gtexGeneSyms);
gtexSNet = gtexNet(b, b);

gtexBGS = ...  % blood gene symbol
    GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(1).expGenes);
gtexBN =  GTExFiveNets.nets(1).net005;
gtexBND = sum(gtexBN' + gtexBN');

gtexLGS = ... % liver gene symbol
    GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(3).expGenes);
gtexLN =  GTExFiveNets.nets(3).net005;
gtexLND = sum(gtexLN' + gtexLN');

%gtexSNet = gtexNet(gtexND > 2, gtexND > 2);

% get the ctcNet 
load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
 ...
      'correctedBinNets_logCorrected_jusResiduals_scBasedMarkers_7PCA.mat'])

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'correctedBinNets_logCorrected_jusResiduals_scBasedMarkers_redo.mat'])

[a, b] = ismember(gtexCluster.syms, ctc.geneSyms);
ctcNet = ctc.net005(b, b);

% the rest

myInds = find(highCount);
myMeanR2s = meanR2(highCount);
gCounts = zeros(1, length(myInds));
ctcDensity = zeros(1, length(myInds));
ctclc = zeros(1, length(myInds));
ctcgc = zeros(1, length(myInds));
gtexDensity = zeros(1, length(myInds));
gtexlc = zeros(1, length(myInds));
TANDensity = zeros(1, length(myInds));
TANlc = zeros(1, length(myInds));
TANgc = zeros(1, length(myInds));
TSNDensity = zeros(1, length(myInds));
TSNlc = zeros(1, length(myInds));
% TSNgc = zeros(1, length(myInds)); % this is same as TAN GC
bloodTANDensity = zeros(1, length(myInds));
bloodTANlc = zeros(1, length(myInds));
bloodTANgc = zeros(1, length(myInds));
liverTANDensity = zeros(1, length(myInds));
liverTANlc = zeros(1, length(myInds));
liverTANgc = zeros(1, length(myInds));
GBD = zeros(1, length(myInds));
GBlc = zeros(1, length(myInds));
GBgc = zeros(1, length(myInds));
GLD = zeros(1, length(myInds));
GLlc = zeros(1, length(myInds));
GLgc = zeros(1, length(myInds));

cs = gtexCluster.cs1547;
for i = 1:length(myInds)
    tempIDs = find(cs == myInds(i));
    length(tempIDs)
    
    myGenes = gtexCluster.syms(tempIDs);
    gCount = length(myGenes)
    gCounts(i) = gCount;
    
    % ctc density
    ctcDensity(i) = (sum(sum(ctcNet(tempIDs, tempIDs)))*2) / (gCount * (gCount - ...
                                                      1));
    ctcgc(i) = gCount;
    ctclc(i) = sum(sum(ctcNet(tempIDs, tempIDs)));

    % Density in GTEx:
    gtexDensity(i) = (sum(sum(gtexSNet(tempIDs, tempIDs)))*2) / (gCount * (gCount - ...
                                                      1));
    gtexlc(i) = sum(sum(gtexSNet(tempIDs, tempIDs)));

    [a, b] = ismember(myGenes, gpl570.uniqueSymbols(affyBloodExpG));
    
    bloodTANgc(i) = sum(a);
    subBloodTAN = bloodTAN(affyBloodExpG, affyBloodExpG);
    affyBloodS = subBloodTAN(b(a), b(a));
    bloodTANlc(i) = sum(affyBloodS(:));
    bloodTANDensity(i) = bloodTANlc(i) / (sum(a)*(sum(a)-1)/ ...
                                               2);
    
    
    [a, b] = ismember(myGenes, gpl570.uniqueSymbols(affyLiverExpG));
    
    liverTANgc(i) = sum(a);
    subLiverTAN = liverTAN(affyLiverExpG, affyLiverExpG);
    affyLiverS = subLiverTAN(b(a), b(a));
    liverTANlc(i) = sum(affyLiverS(:));
    liverTANDensity(i) = liverTANlc(i) / (sum(a)*(sum(a)-1)/ ...
                                               2);

    
    clear a b
    [a, b] = ismember(myGenes, gpl570.uniqueSymbols(affyBrainExpG));
    TANgc(i) = sum(a);
    subTAN = tan(affyBrainExpG, affyBrainExpG);
    TANs = subTAN(b(a), b(a));
    TANlc(i) = sum(TANs(:));
    TANDensity(i) = TANlc(i) / (sum(a)*(sum(a)-1)/2);
    
    subTSN = tsn(affyBrainExpG, affyBrainExpG);
    TSNs = subTSN(b(a), b(a));
    TSNlc(i) = sum(TSNs(:));
    TSNDensity(i) = TSNlc(i) / (sum(a)*(sum(a)-1)/2);
    
    clear a b
    [a, b] = ismember(myGenes, gtexBGS);
    GBgc(i) = sum(a);
    tempSmallNet = gtexBN(b(a), b(a));
    GBlc(i) = sum(tempSmallNet(:));
    GBD(i) = GBlc(i) / (sum(a)*(sum(a)-1)/2);

    clear a b tempSmallNet
    [a, b] = ismember(myGenes, gtexLGS);
    GLgc(i) = sum(a);
    tempSmallNet = gtexLN(b(a), b(a));
    GLlc(i) = sum(tempSmallNet(:));
    GLD(i) = GLlc(i) / (sum(a)*(sum(a)-1)/2);
end

GTExClusterRep_itself.clusterInds = myInds;
GTExClusterRep_itself.clusterGC = gCounts;
GTExClusterRep_itself.lc = gtexlc;
GTExClusterRep_itself.gc = gCounts;
GTExClusterRep_itself.density = gtexDensity;
save('~/resultsAndFigures/secondProject/gtexClusterRep/GTExitselfClusterRep_V2.mat', ...
     'GTExClusterRep_itself')

GTExClusterRep_ctc.clusterInds = myInds;
GTExClusterRep_ctc.clusterGC = gCounts;
GTExClusterRep_ctc.lc = ctclc;
GTExClusterRep_ctc.gc = ctcgc;
GTExClusterRep_ctc.density = ctcDensity;
save('~/resultsAndFigures/secondProject/gtexClusterRep/GTExCTCClusterRep_SC7PCA_redo.mat', ...
     'GTExClusterRep_ctc')

% saving reproducibility results
GTExClusterRep_Gblood.clusterInds = myInds;
GTExClusterRep_Gblood.clusterGC = gCounts(1, :);
GTExClusterRep_Gblood.lc = GBlc;
GTExClusterRep_Gblood.gc = GBgc;
GTExClusterRep_Gblood.density = GBD;
save('~/resultsAndFigures/secondProject/gtexClusterRep/GTExBloodClusterRep_V2.mat', ...
     'GTExClusterRep_Gblood')

% saving reproducibility results
GTExClusterRep_Gliver.clusterInds = myInds;
GTExClusterRep_Gliver.clusterGC = gCounts(1, :);
GTExClusterRep_Gliver.lc = GLlc;
GTExClusterRep_Gliver.gc = GLgc;
GTExClusterRep_Gliver.density = GLD;
save('~/resultsAndFigures/secondProject/gtexClusterRep/GTExLiverClusterRep_V2.mat', ...
     'GTExClusterRep_Gliver')

% saving reproducibility results
GTExClusterRep_tan.clusterInds = myInds;
GTExClusterRep_tan.clusterGC = gCounts(1, :);
GTExClusterRep_tan.lc = TANlc;
GTExClusterRep_tan.gc = TANgc;
GTExClusterRep_tan.density = TANDensity;
save('~/resultsAndFigures/secondProject/gtexClusterRep/TANClusterRep_V2.mat', ...
     'GTExClusterRep_tan')

% saving reproducibility results
GTExClusterRep_tsn.clusterInds = myInds;
GTExClusterRep_tsn.clusterGC = gCounts(1, :);
GTExClusterRep_tsn.lc = TSNlc;
GTExClusterRep_tsn.gc = TANgc;
GTExClusterRep_tsn.density = TSNDensity;
save('~/resultsAndFigures/secondProject/gtexClusterRep/TSNClusterRep_V2.mat', ...
     'GTExClusterRep_tsn')

% saving reproducibility results
GTExClusterRep_tanBlood.clusterInds = myInds;
GTExClusterRep_tanBlood.clusterGC = gCounts(1, :);
GTExClusterRep_tanBlood.lc = bloodTANlc;
GTExClusterRep_tanBlood.gc = bloodTANgc;
GTExClusterRep_tanBlood.density = bloodTANDensity;
save('~/resultsAndFigures/secondProject/gtexClusterRep/affyBloodClusterRep_V2.mat', ...
     'GTExClusterRep_tanBlood')

% saving reproducibility results
GTExClusterRep_tanLiver.clusterInds = myInds;
GTExClusterRep_tanLiver.clusterGC = gCounts(1, :);
GTExClusterRep_tanLiver.lc = liverTANlc;
GTExClusterRep_tanLiver.gc = liverTANgc;
GTExClusterRep_tanLiver.density = liverTANDensity;
save('~/resultsAndFigures/secondProject/gtexClusterRep/affyLiverClusterRep_V2.mat', ...
     'GTExClusterRep_tanLiver')

scDensity =  zeros(76, length(myInds));
sclc =  zeros(76, length(myInds));
scGC =  zeros(76, length(myInds));
cellTypeNames = cell(1,69);
expGeneCount = zeros(1, 69);
for i = 1:76
    i
    load(sprintf(['~/networks/allenBrainSC/allFiveNets/' ...
                  'binNets_Exon_V4_allFive_net%d.mat'], i)) 
    expGeneCount(i) = sum(net.expGenes)
    cellTypeNames(i) = net.clusterName;
    scNet = net.net01;
    scSyms = filDataSet.geneSyms(net.expGenes);
    clear net
    
    for j = 1:length(myInds)
        j
        tempIDs = find(cs == myInds(j));
        length(tempIDs)
        myGenes = gtexCluster.syms(tempIDs);
        gCount = length(myGenes)
        gCounts(i, j) = gCount;
        
        [a, b] = ismember(myGenes, scSyms);
        tgc = sum(a);
        scGC(i, j) = tgc;
        smallSC = scNet(b(a), b(a));
        scDensity(i, j) = sum(smallSC(:)) / (tgc * (tgc -1)/2);
        sclc(i, j) = sum(smallSC(:));
    end
end

% saving SC results
scGTExClusterRep.expGeneCount = expGeneCount;
scGTExClusterRep.clusterInds = myInds;
scGTExClusterRep.clusterGC = gCounts(1, :);
scGTExClusterRep.lc = sclc;
scGTExClusterRep.gc = scGC;
scGTExClusterRep.density = scDensity;
scGTExClusterRep.ctn = cellTypeNames;

save('~/resultsAndFigures/secondProject/gtexClusterRep/scClusterRep.mat', ...
     'scGTExClusterRep')

load('~/resultsAndFigures/secondProject/gtexClusterRep/scClusterRep.mat')

% 3.5. comparing the GTEx cluster with Simulated networks 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('~/data/brainSingleCell/scGeneIDfier.mat')
% geneSyms = sc.updateGeneSyms;

load('~/resultsAndFigures/secondProject/gtexClusterRep/GTExBloodClusterRep_V2.mat')
myInds = GTExClusterRep_Gblood.clusterInds;


cs = gtexCluster.cs1547;

files = dir(['~/resultsAndFigures/secondProject/' ...
             'SimBulkNetworksFromSC/withCounts/bulk*.mat'])
files = dir(['~/resultsAndFigures/secondProject/' ...
             'SimBulkNetworksFromSC/*180.mat'])

simDensity =  zeros(length(files), length(myInds));
simlc = zeros(length(files), length(myInds));
netFileNames = cell(1, length(files));

gCounts = zeros(1, 69);
scGeneCount = zeros(1, 69);
for j = 1:length(files)
    j
    load([files(j).folder '/' files(j).name])
    myBinNet = bulkFromSC.binNet005;
    %    geneSyms = bulkFromSC.geneSyms;
    for i = 1:length(myInds)
        i
        tempIDs = find(cs == myInds(i));
        length(tempIDs)
        
        myGenes = gtexCluster.syms(tempIDs);
        gCount = length(myGenes)
        gCounts(i) = gCount;
        
        [a, b] = ismember(myGenes, geneSyms);
        sum(a)
        scGeneCount(i) = sum(a);
        % simDensity
        simDensity(j, i) = (sum(sum(myBinNet(b(a), b(a))))) / (scGeneCount(i) * (scGeneCount(i) - ...
                                                          1));
        simlc(j, i) = (sum(sum(myBinNet(b(a), b(a)))));
    end
end

lowVarLC = simlc;
highVarLC15 = simlc;
highVarLC12 = simlc;

j = 1
load([files(j).folder '/' files(j).name])

(clusterMeta.clusters(bulkFromSC.selectedSamples{1}))


%%%%%%%%%%%%
% >> I need to fill in the gaps: 
% where are the cell type markers,
% what are the functions
% correlation of the markers (not in the heatmpa, kind of side
% product)
% division of the big clusters. 

% get the distriubtion of r2 for the cluster genes
[a, b] = ismember(gtexCluster.syms, dataSet.genes);
clustersR2 = result.r2(b);
h = figure
hist(clustersR2, 100)
h = figure
hist(result.r2, 100)

bigR2s = find(myMeanR2s > .4);
smallR2s = find(myMeanR2s < .2);

h = figure
heatmap(myMeanR2s)

h = figure
inds = myMeanR2s > .4;
heatmap([gtexDensity(inds).*2; TANDensity(inds); TSNDensity(inds); ...
         bloodTANDensity(inds)])

% 4. The base clustering plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% scGTExClusterRep
load('~/resultsAndFigures/secondProject/gtexClusterRep/scClusterRep.mat')
cellTypeNames = scGTExClusterRep.ctn;

% 4.1 adding the count of genes in each cluster to the  myInds 
genesInCluster = zeros(1, length(myInds));
clusterLabels = cell(1, 69);
for i = 1:length(myInds)
    genesInCluster(i) = length(find(cs == myInds(i)));
    clusterLabels{i} = sprintf('ID%d - %d', myInds(i), genesInCluster(i));
end

% 4.2 adding the network names for each network 
netLabels = cell(1, 76);
netLabels{1} = 'Rsqrd';
netLabels{2} = 'GTEx_cortex';
netLabels{3} = 'TAN';
netLabels{4} = 'TSN';
netLabels{5} = 'TAN_blood';
netLabels{6} = 'GTEx_blood';
netLabels{7} = 'GTEx_liver';
for i = 1:length(cellTypeNames)
    netLabels{(i+7)} = sprintf('%s -%d', cellTypeNames{i}, i);
end

netLabels(8:end) = cellTypeNames;

ctInds = [1:7, ([1:6, 8,10:12, 21:24,26:28, 7 9, 13:15,19,20 25, 29, 30, ...
          16:18] + 7)];
% 4.3 plotting 
plotMat = [myMeanR2s; GBD; GLD; bloodTANDensity; TANDensity; TSNDensity; ...
           scDensity(1:69, :)];
h = figure
heatmap(clusterLabels, [1:75], plotMat)
colormap(jet)

plotMat = [myMeanR2s;gtexDensity.*2; TANDensity; TSNDensity; ...
         bloodTANDensity; GBD; GLD; scDensity(1:69, :)];

densPlot.plotMat = plotMat;
densPlot.clusterLabels = clusterLabels;
densPlot.netLabels = netLabels;

save(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'clusterDensityPlot'], 'densPlot')

h = figure

heatmap(clusterLabels, netLabels(ctInds), plotMat(ctInds, :))
heatmap(clusterLabels, netLabels, plotMat)

colormap(jet)
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sclusterRep_30SC', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

plotMat = [myMeanR2s;gtexDensity.*2; TANDensity; TSNDensity; ...
         bloodTANDensity; scDensity];
sib = corr(plotMat', 'type', 'Spearman', 'rows','pairwise');
h = figure
heatmap(sib)

h = figure
inds = myMeanR2s < .4;
heatmap([gtexDensity(inds).*2; TANDensity(inds); TSNDensity(inds); ...
         bloodTANDensity(inds)])

% 5. Enrichment of functional terms in each sc (GTEx) cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('~/data/general/GOdata_GTExV6_nonComp.mat')
GTExGO = GOdata;
clear GOdata;

[a, b] = ismember(gtexCluster.syms,GTExGO.geneSymbols);
wholeFMat = GTExGO.matP(b, :);
% terms with <=200 and >= 5 genes
sib = sum(wholeFMat);
inTerms = ((sib >=3) + (sib<=200)) == 2;

fMat = wholeFMat(:, inTerms);
inTermNames = GTExGO.GOTerms(inTerms);
inTermsGOID = GTExGO.GOID(inTerms);
rawGeneCount = sib(inTerms);
expGeneCounts = sum(fMat);

cs = gtexCluster.cs1547;
csGenes = gtexCluster.syms;
clusterCount = max(cs);
fCount = length(inTermNames);

% members of each cluster 
cmc = hist(cs, [1:max(cs)]); % clusterMemberCount

inClusters = find(cmc >= 20);
inClusterCount = length(inClusters);

% getting the count of genes in each of the clusters
cfMat = zeros(inClusterCount, fCount);
for i = 1:length(inClusters)
    myGenes = (cs == inClusters(i));
    smallFMat = fMat(myGenes, :);
    cfMat(i, :) = sum(smallFMat);
end

gCounts= sum(fMat);
totalGCount = length(dataSet.genes)
geneCountInCluster = cmc(inClusters);
ps = ones(length(inClusters), fCount) * -1;
for i = 1:length(inClusters)
    % count of genes in each cluster with f 
    inFs = cfMat(i, :) >=3;
    ps(i, inFs) = 1 - hygecdf(cfMat(i, inFs) - 1, totalGCount, gCounts(inFs), ...
                           geneCountInCluster(i));
end

fCounts = sum(inTerms)
% printing list of terms enriched in each cluster: 
selectedPs = ((ps >=0) + ((ps* fCounts) <= .1)) == 2;
kado = sum(selectedPs');
h = figure
plot(sort(kado))

% saving the results:
clusterFRes.ps = ps
clusterFRes.totalGCounts = gCounts;
clusterFRes.geneCountEachCluster = cfMat;
clusterFRes.clusterIDs = inClusters;
clusterFRes.fTerms = inTermNames;
clusterFRes.fIDs = inTermsGOID;
save(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'functionalEnrichment.mat'], 'clusterFRes')

load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'functionalEnrichment.mat'])

myMat = clusterFRes.ps * length(clusterFRes.fTerms);

ins = ((myMat < .1) + (myMat >= 0)) == 2;
book = sum(ins');

thisTerms = clusterFRes.fTerms(ins(19, :));
thisTerms(:)
thisIDs = clusterFRes.fIDs(ins(68, :));

% 5.5 hkg in clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('~/resultsAndFigures/secondProject/gtexClusterRep/GTExBloodClusterRep_V2.mat')
myInds = GTExClusterRep_Gblood.clusterInds;

load('~/data/general/hkgInd.mat')
load('~/data/general/GPL570GemmaMapNEW.mat')

hkgSyms = gpl570.uniqueSymbols(hkgInd);

% count of HKG in clusters and what functional terms do they have
hkgPr = zeros(length(hkgInd), 69);
for i = 1:length(myInds) 
   [a, b] = ismember(hkgSyms, gtexCluster.syms(gtexCluster.cs1547 == ...
                                               myInds(i)));
   hkgPr(:, i) = a;
end

book = sum(hkgPr');
sum(book > 0)
h = figure
plot(sum(hkgPr), 'o')
set(gca, 'XTick', [1:69],'XTickLabel', gtexCluster.clusterLabels(myInds))
xtickangle(60)

% for clusters with >= 10 HKG, what are the functions enriched in
% those terms
for i = 1:69
    i
    if sum(hkgPr( :, i)) >= 10;
        [passedFDRIDs{i}, passedFDRTerms{i}, inGenes{i}, enrichMat{i}] = ...
            fenrich_function(hkgSyms(logical(hkgPr(:, i))), .1);
    end
end

passedFDRL = zeros(1, 69);
for i = 1:69
    passedFDRL(i) = length(passedFDRIDs{i});
end

meanR2(myInds(passedFDRL>0))

i = 9
finGenes = inGenes{i}(sum(enrichMat{i}') > 0)

h = figure
heatmap(passedFDRTerms{i}, inGenes{i}, enrichMat{i}, 'GridVisible', ...
        'off')

% 6. Studying the clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clusters are already sorted based on the R2 values. I will look
% at 4 other things in each cluster: 1. reproducibility, 2. markers
% 3. functional terms 4. protein complexes 

% 0. put the gene count in the name of clusters and name the
% networks 
% >>> This is done in section 4. 

% >>>>>>>>>>>>>>>>>> 0- fetch the rep results and find it out for each cluster for
% each network 
% >>> Single cell networks

load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'scClusterRep.mat'])

myPs = zeros(69, 69);
allTestValues = zeros(69, 69);
for i = 1:69  % i is for the SC cluster/network
             % loading linkCounts
    load(sprintf(['~/networks/allenBrainSC/randomNetLC_perCluster/' ...
                  'gtexCluster_1547/rpm_clusterLinkCount_cellType%d_10k.mat'], ...
                 i ))
    for j = 1:69 % j is for the gtexCluster (here we only need that
                 % count of genes in that cluster) 
        d = scGTExClusterRep.gc(i, j);
        if d >= 5
            myInd = find(linkCounts.geneCounts == d);
            
            myPs(i, j) = sum(linkCounts.thr01(myInd, :) >= scGTExClusterRep.lc(i,j))/length(linkCounts.thr01);
        else 
            myPs(i,j)= nan;
        end
    end
    clear linkCounts
end

h = figure
heatmap(myPs)

% >>>>>>>>> GTEx networks

% GTEx blood
load('~/resultsAndFigures/secondProject/gtexClusterRep/GTExBloodClusterRep.mat')
load(['~/networks/allenBrainSC/randomNetLC_perCluster/' ...
      'gtexCluster_1547/gtexBlood_clusterLinkCount.mat'])

GBPs = zeros(1, 69);
for i = 1:69
    if linkCounts.geneCounts(i) >=5
        GBPs(i) = sum(linkCounts.thr005(i, :) > ...
                        GTExClusterRep_Gblood.lc(i))/length(linkCounts.thr005);
    else
        GBPs(i) = nan;
    end
end

GBFDR = GBPs * 69;

% GTEx liver
load('~/resultsAndFigures/secondProject/gtexClusterRep/GTExLiverClusterRep.mat')
load(['~/networks/allenBrainSC/randomNetLC_perCluster/' ...
      'gtexCluster_1547/gtexLiver_clusterLinkCount.mat'])

GLPs = zeros(1, 69);
for i = 1:69
    if linkCounts.geneCounts(i) >=5
        GLPs(i) = sum(linkCounts.thr005(i, :) > ...
                     GTExClusterRep_Gliver.lc(i))/length(linkCounts.thr005);
    else
        GLPs(i) = nan;
    end
end
GLFDR = GLPs * 69;

% GTEx CTC
clear linkCounts
load('~/resultsAndFigures/secondProject/gtexClusterRep/GTExCTCClusterRep_V2.mat')
load(['~/networks/allenBrainSC/randomNetLC_perCluster/' ...
      'gtexCluster_1547/gtexCTC_clusterLinkCount.mat'])

GCPs = zeros(1, 69);
for i = 1:69
    if linkCounts.geneCounts(i) >=5
        GCPs(i) = sum(linkCounts.thr005(i, :) > ...
                        GTExClusterRep_ctc.lc(i))/length(linkCounts.thr005);
    else
        GCPs(i) = nan;
    end
end
GCFDR = GCPs * 69;

h = figure 
heatmap(GCPs)
% >>>>>>>>> TAN networks
% get the count of genes for each cluster in that network 
clear linkCounts
load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'affyLiverClusterRep.mat'])

load(['~/networks/allenBrainSC/randomNetLC_perCluster/' ...
      'gtexCluster_1547/liverTAN_clusterLinkCount.mat'])

tanLiverPs = zeros(1, 69);
for i = 1:69
    if GTExClusterRep_tanLiver.gc(i) >=5
        tanLiverPs(i) = sum(linkCounts(i, :) > ...
                        GTExClusterRep_tanLiver.lc(i))/length(linkCounts(i, ...
                                                          :));
    else
        tanLiverPs(i) = nan;
    end
end
tanLiverFDR = tanLiverPs * 69;

% TAN blood
clear linkCounts
load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'affyBloodClusterRep.mat'])

load(['~/networks/allenBrainSC/randomNetLC_perCluster/' ...
      'gtexCluster_1547/bloodTAN_clusterLinkCount.mat'])

tanBloodPs = zeros(1, 69);
for i = 1:69
    if GTExClusterRep_tanBlood.gc(i) >=5
        tanBloodPs(i) = sum(linkCounts(i, :) > ...
                        GTExClusterRep_tanBlood.lc(i))/length(linkCounts(i, ...
                                                          :));
    else
        tanBloodPs(i) = nan;
    end
end
tanBloodFDR = tanBloodPs * 69;
sum(tanBloodFDR < .1)

% TAN brain
clear linkCounts
load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'TANClusterRep.mat'])

load(['~/networks/allenBrainSC/randomNetLC_perCluster/' ...
      'gtexCluster_1547/brainTAN_clusterLinkCount.mat'])

tanBrainPs = zeros(1, 69);
for i = 1:69
    if GTExClusterRep_tan.gc(i) >=5
        tanBrainPs(i) = sum(linkCounts(i, :) > ...
                        GTExClusterRep_tan.lc(i))/length(linkCounts(i, ...
                                                          :));
    else
        tanBrainPs(i) = nan;
    end
end
tanBrainFDR = tanBrainPs * 69;
sum(tanBrainFDR < .1)

% >>>>>> TSN brain
clear linkCounts
load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'TSNClusterRep.mat'])

load(['~/networks/allenBrainSC/randomNetLC_perCluster/' ...
      'gtexCluster_1547/brainTSN_clusterLinkCount.mat'])

tsnBrainPs = zeros(1, 69);
for i = 1:69
    if GTExClusterRep_tsn.gc(i) > 5
        tsnBrainPs(i) = sum(linkCounts(i, :) > ...
                        GTExClusterRep_tsn.lc(i))/length(linkCounts(i, ...
                                                          :));
    else
        tsnBrainPs(i) = nan;
    end
end
tsnBrainFDR = tsnBrainPs * 69;
sum(tsnBrainFDR < .1)

book = [GBPs; GLPs; tanBloodPs; tanLiverPs; tanBrainPs; tsnBrainPs; GCPs; ...
         myPs];

allPs.mat = book;
allPs.clusterLabels = densPlot.clusterLabels;
allPs.netLabels = {'GTEx_blood', 'GTEx_liver', 'TAN_blood', ...
                   'TAN_liver', 'TAN_brain', 'TSN_brain', 'GTEx_CTC'}
allPs.netLabels = [allPs.netLabels, densPlot.netLabels(7:75)];

save(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'repPValues_76Nets_structure.mat'], 'allPs');
book = allPs * 60;
h = heatmap(allPs)

sum(allPs' == 0)

% 6.2 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
clear
% I want for each cluster: 
% 0- setup the heatmap plot for the reproducibility
load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'GTExBloodClusterRep.mat']) % GTExClusterRep_Gblood.clusterInds forthe
                                  % cluster index

load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'clusterDensityPlot'])

load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'repPValues_75Nets.mat']);
h = figure('units', 'centimeters', 'position', [0,0, 33, 20])
netLabels = ['GTExBlood', 'GTExLiver', 'tanBlood', 'tanLiver', 'tan', ...
             'tsn', densPlot.netLabels(7:end)];
heatmap(densPlot.clusterLabels, netLabels,allPs, 'FontSize', 7)
title('pvalues for the reproducibility of the clusters')
figFolder = ['~/resultsAndFigures/secondProject/clusterProfiling/']
file = sprintf('%sGTEx1547ClusterReproducibility', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

h = figure
heatmap(allPs(:, faceThic

% plotting the density plot
load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'clusterDensityPlot'])

h = figure('units', 'centimeters', 'position', [0,0, 35, 18])
heatmap(densPlot.clusterLabels, [1:76], ...
        densPlot.plotMat)
title('density of the networks')
colormap(jet)

% plotting the geneCount for all the SC clusters 
load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'scClusterRep.mat']) % scGTExClusterRep
load('~/resultsAndFigures/secondProject/gtexClusters.mat') % 'gtexCluster'
tempClusterGeneCounts = hist(gtexCluster.cs1547, ...
                             unique(gtexCluster.cs1547));
clusterGeneCounts = tempClusterGeneCounts(tempClusterGeneCounts >= ...
                                          20);
expGenePercent = zeros(75, 69);
for i = 1:69
    expGenePercent(7:75, i) = scGTExClusterRep.gc(1:69, i) ./ clusterGeneCounts(i);
end
expGenePercent(1, :) = GTExClusterRep_Gblood.gc ./ ...
    clusterGeneCounts;
expGenePercent(2, :) = GTExClusterRep_Gliver.gc ./ ...
    clusterGeneCounts;
expGenePercent(3, :) = GTExClusterRep_tanBlood.gc ./ ...
    clusterGeneCounts;
expGenePercent(4, :) = GTExClusterRep_tanLiver.gc ./ ...
    clusterGeneCounts;
expGenePercent(5, :) = GTExClusterRep_tan.gc ./ ...
    clusterGeneCounts;
expGenePercent(6, :) = GTExClusterRep_tsn.gc ./ ...
    clusterGeneCounts;

h = figure('units', 'centimeters', 'position', [0,0, 35, 18])
heatmap(densPlot.clusterLabels, netLabels, expGenePercent);
title('percent of genes expressed')

h = figure('units', 'centimeters', 'position', [0,0, 33, 20])
netLabels = ['GTExBlood', 'GTExLiver', 'tanBlood', 'tanLiver', 'tan', ...
             'tsn', densPlot.netLabels(7:end)];
heatmap(densPlot.clusterLabels, netLabels, expGenePercent, 'FontSize', ...
        7);
title('Percent of genes expressed')
figFolder = ['~/resultsAndFigures/secondProject/clusterProfiling/']
file = sprintf('%sGTEx1547ClusterGeneExp', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

myInds = GTExClusterRep_Gblood.clusterInds;
% get the clusterID
% ----------------
thic = 202% this cluster ID
faceThic = find(myInds == thic)
% ----------------

% 1- the gene count, 2- the gene list, 3- the sub clusters (for the big ones)
% load the clustering files (both of them)
load('~/resultsAndFigures/secondProject/gtexClusters.mat') % 'gtexCluster'

geneCount = sum(gtexCluster.cs1547 == thic)
geneList = gtexCluster.syms(gtexCluster.cs1547 == thic);
length(geneList)
subClusters = gtexCluster.cs1546(gtexCluster.cs1547 == thic);
subgCounts = hist(subClusters, unique(subClusters));
% getting a specific gene list from sub clusters
find(subgCounts > 500)
book = unique(subClusters);
book(find(subgCounts > 30))
subGeneList = geneList(subClusters == 112);
geneList = subGeneList
geneList = subGeneList;
h = figure
hist(subClusters, unique(subClusters))

% this is the dense part of cluster 242 in TSN
save(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'tsnDenseGeneListInCluster242_cs1547.mat'], ...
     'tsnDenseGeneList')

% do the sub clusters present different meanR2, and if they do,
% does it affect their rep in the other tissues. 
% so basically, why do I see this cluster reproduced in other
% networks and what does it mean

% 1.25 - subcluster investigation
% TODO get the ID for clusters with > 30 genes
% TODO get the similar info for the sub clusters (I just need to
% save the cluster in the GTEx cluster format and run the code for it)

% 1.5 - count of genes in each network
load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'scClusterRep.mat']) % scGTExClusterRep
h = figure
heatmap(1, netLabels(7:end),(scGTExClusterRep.gc(1:69, faceThic))./scGTExClusterRep.clusterGC(faceThic)');
heatmap(1, netLabels(7:end),(scGTExClusterRep.lc(1:69, faceThic)));
heatmap(1, netLabels(7:end),(scGTExClusterRep.gc(1:69, faceThic)));

% 3.5- the r2 % for this one I need the whole dataSte for the gene
% Symbols 
% load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
%       'regressionResult_5PCA_10CellTypes_MarkersV2_cthr7.mat']) %
%       not good, I switched to a more stringent filter

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_MarkersV1_cthr6.mat'])

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_AllenSCBasedMarker.mat'])

load('~/data/GTEx/Brain_Cortex_expGenes.mat') % 'dataSet'

% get genes with the top 20% r2

[a, b] = ismember(geneList, dataSet.genes);
myR2s = result.r2(b);
qs = quantile(myR2s, 100);
r2Selected = find(myR2s >= qs(80));
highR2Genes = geneList(r2Selected);

coeffMat = result.coeff(b, :);
h = figure
heatmap([1:5], geneList, (coeffMat(:, 2:end)))

[a, b] = ismember(geneList, dataSet.genes);
[a, b] = ismember(presentMarkers, dataSet.genes);
myRs = result.r2(b);
h = figure
hist(myRs)

h = figure
heatmap(1, geneList, myRs')
[sa, sb] = sort(myRs, 'descend');
h = figure
heatmap(1, geneList(sb), myRs(sb)')

r2SortedGL = geneList(sb);

% 3.75 marker genes

load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])
                                               
% does the cluster have any markers and what are they. 
markerCounts = zeros(1, 10);
presentMarkers = zeros(1, length(geneList));
for i = 1:10
    thisMarkerGenes = markers.genes(markers.finalMarkerGenes{i});
    [a, b] = ismember(thisMarkerGenes, geneList);
    markerCounts(i) = sum(a);
    presentMarkers(b(a)) = i;
end

h = figure
plot(markerCounts)

% identify the marker genes present in the geneList for a give cell
% type:
thisMarkerGenes = markers.genes(markers.finalMarkerGenes{1});
[a, b] = ismember(thisMarkerGenes, geneList);
presentMarkers = thisMarkerGenes(a);

[a, b] = ismember(presentMarkers, tsnGeneListNDBasedSorted);

% 4- the functions enriched
load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'functionalEnrichment.mat']) % 'clusterFRes'
fCounts = length(clusterFRes.fTerms)
selectedPs = ((clusterFRes.ps >=0) + ((clusterFRes.ps* fCounts) <= .1)) == 2;
posF = sum(selectedPs') >= 3;

h = figure
plot(sort(kado))

myFunctionList = clusterFRes.fTerms(selectedPs(faceThic, :)); 
myFunctionInds = find(selectedPs(faceThic,:));

% second way: using the fenrich function
[passedFDRIDs, passedFDRTerms, inGenes, enrichMat] = ...
    fenrich_function(geneList, .1);
whos inGenes
% cluster 239 two three nine
ttn.IDs = passedFDRIDs;
ttn.terms = passedFDRTerms;
ttn.genes = inGenes;
ttn.mat = enrichMat;

% getting clusters with >=3 functions
fSelectedCs = kado >=3;
kado = sum(allPs <= (.1/69));

repSelectedCs = kado >= 5;

h = figure
heatmap(densPlot.clusterLabels, [1:3],[densPlot.plotMat(1,:); fSelectedCs; repSelectedCs]+ 0)
colormap(jet)

% for a given function, how many genes have it and how many of them
% are in my cluster
ppf = clusterFRes.geneCountEachCluster(faceThic, myFunctionInds) ./ ...
      clusterFRes.totalGCounts(myFunctionInds);

find(ppf >.3)
selectedFIDs = clusterFRes.fIDs(myFunctionInds(ppf >.3));
selectedFNames = myFunctionList(ppf>.3);

% 5- the functions represented
% GTEx bin, tan, tsn, GTEx blood, GTEx liver functions represented
load(['~/resultsAndFigures/secondProject/moduleRepresentation/' ...
      'funcRes.mat'])

% getting the brain specific and non specific terms (in moduleRepresentation.m)
load(['~/resultsAndFigures/secondProject/' ...
      'moduleRep_GTExBlood.mat']) % gtexBinNetBlood

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

% SC function representation
load(['~/resultsAndFigures/secondProject/' ...
      'moduleRep_netSCAll_exon3_V4.mat']) % ''

% 6- the complexes represented
load(['~/resultsAndFigures/secondProject/networkComplexes/' ...
      'complexPresentation_fiveNets_Exons_V4.mat']) % complexRes

% get name of all complexes
allCompNames = cell(1, 2003);
for i = 1:2003
    allCompNames{i} = complexRes.complexes(i).Name;
    book = strfind(allCompNames{i}, 'ATP');
    if book 
        i
    end
end
kado = unique(allCompNames);

genesInComplexes = zeros(1, 2003);
for i = 1:2003 % count of complexes
    [a, b] = ismember('ATP6V1F', ...
                      complexRes.complexes(i).atomGeneSyms);
    if(a)
        i
    end
end

% to find out which complexes are present in my cluster and to what
% extent. 
genesInComplexes = zeros(1, 2003);
for i = 1:2003 % count of complexes
    [a, b] = ismember(geneList, ...
                      complexRes.complexes(i).atomGeneSyms);
    genesInComplexes(i) = sum(a);
end

% getting the present complexes
presentComps = find(genesInComplexes > 3);
complexRes.complexes(presentComps).Name

% how many of the genes from each complex
[genesInComplexes(presentComps) ./ ...
    [complexRes.complexes(presentComps).pCount]; complexRes.complexes(presentComps).pCount]

% TODO how many of the links from each complex in each network (what do I have in
% my cluster?)
[complexRes.linkCounts([2, 7, 19, 23:40], presentComps)] % just testing

% for a given complex: (presentComps >= 3)
[presentGenes, linksInNets, linkPortions] = ...
    preoteinComplexes_clusterPresence(presentComps, geneList)

h = figure
heatmap(linkPortions)

% TODO: are they subcomplexes of each other? (give me the parents
% of each complex) 

% Are the complexes enriched in these networks?
tempMat =  complexRes.fdrs([2, 6, 19, 23:end], presentComps);
tempMat(tempMat > .1) = nan;
plotMat = tempMat;
h = figure
heatmap(plotMat)

% the complexes are not present in the TSN, and not so much in the
% SC data (which they should be): are there any substantial
% differences between links present in TSN and TAN and GTEx blood?
% are the links specific? 

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Go to the file >> clustersInNetworks.m to investigate clusters furhter
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% 7- Distribution of R2 for individual clusters
% >>>>>>>>>>>>>>>>>>>>>>>>>>>

gtexCluster
dataSet
result

[a, b] = ismember(gtexCluster.cs1547, myInds);
clusterIDsForPreGenes = gtexCluster.cs1547(a);
[a, b] = ismember(gtexCluster.syms(a), dataSet.genes);
presentGeneInds = b(a); % this is for the result.r2

myR2s = result.r2(presentGeneInds);

h = figure
subplot(2, 1, 1)
boxplot(myR2s, clusterIDsForPreGenes, 'plotstyle', 'compact');
ylabel('R2')

myMeanExps = mean(log2(dataSet.mat(presentGeneInds, :) + 1)');
subplot(2, 1, 2)
boxplot(myMeanExps, clusterIDsForPreGenes, 'plotstyle', 'compact');
ylabel('expression')

% find clusters with high R2
length(myInds(find(densPlot.plotMat(1, :) >= .55)))

% 7.5- Getting results from this cluster 

% given a cluster ID, give me the genes expressed, link counts and
% density in all the networks. 

load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'scClusterRep.mat']) % scGTExClusterRep
scGTExClusterRep
load('~/resultsAndFigures/secondProject/gtexClusterRep/TANClusterRep_V2.mat')
GTExClusterRep_tan
load('~/resultsAndFigures/secondProject/gtexClusterRep/TSNClusterRep_V2.mat')
GTExClusterRep_tsn
load('~/resultsAndFigures/secondProject/gtexClusterRep/GTExBloodClusterRep_V2.mat')
GTExClusterRep_Gblood
load('~/resultsAndFigures/secondProject/gtexClusterRep/GTExLiverClusterRep_V2.mat')
GTExClusterRep_Gliver
load('~/resultsAndFigures/secondProject/gtexClusterRep/affyLiverClusterRep_V2.mat')
GTExClusterRep_tanLiver
load('~/resultsAndFigures/secondProject/gtexClusterRep/affyBloodClusterRep_V2.mat')
GTExClusterRep_tanBlood

load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'GTExitselfClusterRep.mat'])
GTExClusterRep_itself

clusterInds = GTExClusterRep_Gblood.clusterInds; % can get the
                                                   % clusterInds
                                                   % from any of
                                                   % the clusterRep structures.
% myInds
myInds1 = find(densPlot.plotMat(1, :) >= .5)

% markers 
load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])

markerCountMat = zeros(10, 69);
for j = 1:69
    % does the cluster have any markers and what are they. 
    geneList = gtexCluster.syms(gtexCluster.cs1547 == (clusterInds(j)));
    for i = 1:10
        thisMarkerGenes = markers.genes(markers.finalMarkerGenes{i});
        [a, b] = ismember(thisMarkerGenes, geneList);
        markerCountMat(i, j) = sum(a);
    end
end

myInds2 = find(sum(markerCountMat) > 3);

myInds3 = find(allPs(6, :) == 0);
myInds3Opp = find(allPs(6, :) > 0);

finalInds = unique([myInds1, myInds2, myInds3]);

h = figure
heatmap(densPlot.plotMat(1, finalInds), markers.cellTypes, markerCountMat(:, ...
                                                  finalInds))
heatmap(clusterInds(finalInds), markers.cellTypes, markerCountMat(:, ...
                                                  finalInds))

sib = find(max(markerCountMat) > 3)
densPlot.plotMat(1, sib)

lcs = zeros(7, length(clusterInds));
for i = 1:length(clusterInds)
    % h = figure
    % hist(scGTExClusterRep.lc(:, i))
    lcs(:, i) = [GTExClusterRep_itself.lc(i)
                 GTExClusterRep_tan.lc(i)
                 GTExClusterRep_tsn.lc(i)
                 GTExClusterRep_Gblood.lc(i)
                 GTExClusterRep_Gliver.lc(i)
                 GTExClusterRep_tanLiver.lc(i)
                 GTExClusterRep_tanBlood.lc(i)];
end
h = figure
heatmap(log10(lcs))
heatmap(clusterInds(myInds3), [1:6], (lcs(:, myInds3)))

% myInds has the IDs
myInds = GTExClusterRep_Gblood.clusterInds;
markerClusterIDs = [35 53 92 231 244, 252, 65, 129 239 242 243, 125];
[a, clusterInds] = ismember(markerClusterIDs, myInds)
markerClusterLCs = lcs(:, b);
sumLCs = sum(markerClusterLCs');

% now a table for SC 

% plotting the marker presence 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% let's stack them for each of the 7 networks for each of the cell
% types, same goes for genes and others

% >> fix the bottome two elements based on the marker set
plotInds = 12 % marker cluster
fromSC = [7 9 13 14 15]
cellTypeMarker = 'Inhibitory'
% <<
sib = markerClusterLCs(: , plotInds)
sib = [sib; scGTExClusterRep.lc(fromSC, clusterInds(plotInds))]
netLabels = {'GTEx-brain-cortex', 'TAN-brain', 'TSN-brain', 'GTEx-blood', ...
             'GTEx-liver', 'TAN-liver', 'TAN-blood'}
netLabels = [netLabels, clusterMeta.sortedClusterNames(fromSC)']
h = figure
[as, bs ] = sort(sib(1, :), 'descend')
bar(sib(: , bs), 'stack')
set(gca, 'XTickLabel', netLabels)
xtickangle(45)
ylabel('count of links')
legend(gtexCluster.clusterLabels(markerClusterIDs(plotInds(bs))))
title(sprintf(['Count of intra-cluster links for the marker-enriched \n clusters ' ...
       'in different networks - %s'], cellTypeMarker))
set(gca, 'FontSize', 12)
figFolder = ['~/resultsAndFigures/secondProject/clusterProfiling/']
file = sprintf('%slinkCountsInMarkerEnrichedClusters_%s', figFolder, cellTypeMarker)
%set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

gcs = zeros(7, length(clusterInds));
for i = 1:length(clusterInds)
    % h = figure
    % hist(scGTExClusterRep.gc(:, i))
    gcs(:, i) = [GTExClusterRep_itself.gc(i)
                 GTExClusterRep_tan.gc(i)
                 GTExClusterRep_tsn.gc(i)
                 GTExClusterRep_Gblood.gc(i)
                 GTExClusterRep_Gliver.gc(i)
                 GTExClusterRep_tanLiver.gc(i)
                 GTExClusterRep_tanBlood.gc(i)];
end
h = figure
heatmap(log10(gcs))

markerClusterGCs = gcs(:, b);
sumGCs = sum(markerClusterGCs');

gcs(:, 2)

% for the rest of the clusters ... 
[a, b] = ismember(clusterInds, clusterInds(finalInds));

oppInds = find(~a);

lcs = zeros(6, length(oppInds))
for i = 1:length(oppInds)
    % h = figure
    % hist(scGTExClusterRep.lc(:, oppInds(i)))
    lcs(:, i) = [GTExClusterRep_tan.lc(oppInds(i))
                 GTExClusterRep_tsn.lc(oppInds(i))
                 GTExClusterRep_Gblood.lc(oppInds(i))
                 GTExClusterRep_Gliver.lc(oppInds(i))
                 GTExClusterRep_tanLiver.lc(oppInds(i))
                 GTExClusterRep_tanBlood.lc(oppInds(i))]
end
h = figure
heatmap(clusterInds(oppInds), [1:6],lcs)

% get the count of links in GTEx cluster
netSyms = ...
    GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(2).expGenes);
mainLCs = zeros(1, 69);
for i = 1:69
    tempGeneList = gtexCluster.syms(gtexCluster.cs1547 == ...
                                    clusterInds(i));
    [a, b] = ismember(tempGeneList, netSyms);
    subNet = GTExFiveNets.nets(2).net005(b, b);
    mainLCs(i) = sum(subNet(:));
end
 
sum(mainLCs(finalInds))
sum(mainLCs(myInds3))
sum(mainLCs(oppInds))

% marker validation in SC 
scGeneRep = zeros(30, 69);
for i = 1:69
    scGeneRep(:, i) = scGTExClusterRep.gc(labModInds, i)./ scGTExClusterRep.clusterGC(i);
end

h = figure
h = figure('units', 'centimeters', 'position', [0,0, 22, 15])
heatmap(densPlot.clusterLabels(finalInds([2, 3, 5 ,6, 8 9 11:15, 17:end])), thislabs(labModInds), ...
        scGeneRep(:, finalInds([2, 3, 5 ,6, 8 9 11:15, 17:end])))
colormap(bone)
title('Percent of genes expressed in each cluster')
figFolder = ['~/resultsAndFigures/secondProject/clusterProfiling/']
file = sprintf('%sGTEx1547ClusterGeneExp_highR2OrMarkers', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

h = figure
h = figure('units', 'centimeters', 'position', [0,0, 35, 15])
heatmap(densPlot.clusterLabels(oppInds([1:7, 9:19, 21:31, 33, 36, 41, 42, 44:end])), thislabs(labModInds), ...
        scGeneRep(:, oppInds([1:7, 9:19, 21:31, 33, 36, 41, 42, 44:end])))
colormap(bone)
colormap(bone)
title('Percent of genes expressed in each cluster')
figFolder = ['~/resultsAndFigures/secondProject/clusterProfiling/']
file = sprintf('%sGTEx1547ClusterGeneExp_lowerR2', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

h = figure
heatmap(1:48, thislabs(labModInds), ...
        scGeneRep(:, oppInds))

book = sum(scGeneRep(:, oppInds));

% coexpression in SC and significance (based on expression)

% There are these count of rep TSN between these genes 

% also their expression 


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

myGeneList = load(['~/resultsAndFigures/secondProject/' ...
      'veryCorrGenesSyms_GTEx.mat'])
% just checking to see count of links
lc = zeros(1, 69);
gc = zeros(1, 69);
for n = 1:69
    n
    load(sprintf(['~/networks/allenBrainSC/allFiveNets/' ...
                  'binNets_Exon_V4_allFive_net%d.mat'], n))
    myNet = net.net01;
    thisNetSyms = filDataSet.geneSyms(net.expGenes);
    [a, b] = ismember(geneList, thisNetSyms);
    gc(n) = sum(a);
    smallNet = myNet(b(a), b(a));
    lc(n) = sum(smallNet(:));
end

[a, b] = ismember(geneList, dataSet.genes);
myRs = result.r2(b);
h = figure
heatmap(1, thislabs, lc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myInds; % has the ID of clusters 
find(myInds == 57)
find(myInds == 248)

theGenes = cs == 236; % 125 and 170 are the inh clusters
sum(theGenes)

myTerms = inTermNames(logical(selectedPs(66, :)))

for i = 1:length(cellTypeNames)
    netLabels{(i+5)} = sprintf('%s -%d', cellTypeNames{i}, i);
end
% 4.3 plotting 
plotMat = [myMeanR2s;gtexDensity.*2; TANDensity; TSNDensity; ...
         bloodTANDensity; scDensity(1:69, :)];
h = figure
h = figure('units', 'centimeters', 'position', [0,0, 35, 20])
heatmap(clusterLabels, netLabels(1:35), plotMat(1:35, :))
book = (zscore(scGC(1:30, :)));
book = scGC(1:30, :);
normSCGC = zeros(30, 69);
for  i = 1:69
    normSCGC(:, i) = scGC(1:30, i)./genesInCluster(i)';
end

h = figure
heatmap(clusterLabels, ['rsqrd', netLabels(6:35)], [myMeanR2s; ...
                    normSCGC])
colormap(jet)

% 6.2 give me the presence of markers for each cluster (this part
% might be extended to which markers are correlated)
markers; 
majorCID; % ID of the cluster with the major group of markers for
          % each cell type 

% 6.3 give me the functional terms for each cluster: the identity
% of the clusters 
selectedPs; % from section 5

% 6.4 give me the protein complexes

% 6.5 the expression level of genes in each cluster

% 6.6 plotting the distribution of R2 
clusterGenesSyms = gtexCluster.syms;
[a, b] = ismember(clusterGenesSyms, dataSet.genes);

clusterR2s = result.r2(b);

[f, x] = ksdensity(clusterR2s(nullInds));
h = figure
plot(x, f)
hold on
umajor = unique(majorCID);
for i = 1:length(umajor)
    ins = gtexCluster.cs1547 == majorCID(i);
           ins = gtexCluster.cs1547 == 252;
    sum(ins)
    [f, x] = ksdensity(clusterR2s(ins));
    plot(x, f)
    hold on
    i = i +1
end

ins = (gtexCluster.cs1547 == 53) + ... % astro
      (gtexCluster.cs1547 == 252) + ...%oligo
      (gtexCluster.cs1547 == 239) + ... % second main Pyramidal
      (gtexCluster.cs1547 == 242); % thie big 2018 one

nullInds = ~ins;

figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sr2Dists', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 7. Studying the clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 7.1 reproducing GTEx networks in blood 

% load the GTEx blood and liver networks 

% get the permutation for the density: for 1000 samples of the
% genes.  

% $. random testing ... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 3
markerPattern = find(cs == kado(i));
length(markerPattern)
mpGenes = gtexCluster.syms(markerPattern);

[a, b] = ismember(mpGenes, gpl570.uniqueSymbols);
sum(a)

smallNet = binNet(b(a), b(a));
sum(smallNet(:))
smallNet = tan(b(a), b(a));
sum(smallNet(:))
smallNet = tsn(b(a), b(a));
sum(smallNet(:))
i = i + 1

% clustering TSN
highRC = find(meanR2 > .4);

% I want to document the reproducibility of the clusters with high
% or low R2 - their functional relevance and the count of genes and
% links they have - and do we see the correlation of the markers 

% then remains the coexpression of marker genes and clusters 

ribInds = (9110:9220);

length(unique(cs(ribInds)))
hist(cs(ribInds), unique(cs(ribInds)))

% report the main functional theme of the clussters, 

% give me the terms for each cluster: 


% I can have the brain terms (present exclusively in brain) and see
% which cluster they end up being in.

% for bigger clusters, can you dissect them to their subparts, can
% you differentiate the functional terms in each subtype 

% The conclusion is that much of the observed 

% random, just finding a specific gene
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

key = 'GABA'
founds = zeros(1, length(gtexCluster.syms));
for i = 1:length(gtexCluster.syms)
    por = strfind(gtexCluster.syms{i}, key);
    if(length(por))
        founds(i) =1;
    end
end


% 8. Identification of the CT clusters based on the SC markers and
% Ogan Markers 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 

% get the SC markers 
load('~/resultsAndFigures/secondProject/scBasedMarkers_withInh.mat')

% get the Ogan mouse markers 
load('~/data/cellTypeVarianceFiles/finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat')

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_MarkersV2_cthr7.mat'])

% get the GTEx cluster and figure out which cluster is which 
load('~/resultsAndFigures/secondProject/gtexClusters.mat') % 'gtexCluster'

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_AllenSCBasedMarker.mat'])

load('~/data/GTEx/Brain_Cortex_expGenes.mat')
% for each cluster with N > 10, give me R2 
cs = gtexCluster.cs1547;
cgCount = hist(cs, 1:max(cs));
meanR2 = zeros(1, max(cs));
memc = zeros(1, max(cs));
for i = 1:max(cs)
    if(length(find(cs == i)) > 1000)
        i
    end
    
    if (length(find(cs == i)) >= 10)
        mySyms = gtexCluster.syms(cs == i);
        [a, b] = ismember(mySyms, dataSet.genes);
        meanR2(i) = mean(result.r2(b(a)));
        memc(i) = sum(a);
    end
end

% Do we have high R2 without marker? Do we have markers without
% high R2?

% for each set of markers, which clusters have them, and what is
% the R2 for each of the clusters, and how many of the markers does
% the cluster have

% >> for SC
sum(scBasedMarkers.markerMat)
scBasedMarkers.cellTypes
presentMarkerCounts = zeros(1, 7);
for i = 1:7
    cellType = scBasedMarkers.cellTypes(i);
    ctGeneSyms = scBasedMarkers.genes(logical(scBasedMarkers.markerMat(:, ...
                                                      i)));
    [a, b] = ismember(ctGeneSyms, gtexCluster.syms);
    presentMarkerCounts(i) = sum(a);
    presentClusters = gtexCluster.cs1547(b(a));
    presentClusterIDs = unique(presentClusters);
    % h = figure
    % hist(presentClusters, unique(presentClusterIDs));
    % title(sprintf('cellType: %s, totalMarkers: %d, presentMarkers: %d', ...
    %               scBasedMarkers.cellTypes{i}, length(ctGeneSyms), ...
    %               sum(a)))
    
    % I plot all clusters which have >= 5% of the markers in them 
    myBar = hist(presentClusters, unique(presentClusterIDs));
    plotClusterInds = find(myBar >= (.05 * presentMarkerCounts(i)));
    
    plotClusterIDs = presentClusterIDs(plotClusterInds); % Xtick
                                                         % label of
                                                         % the plot
    
    % text in plot the R2 of the cluster, 
    h = figure
    plotBar = myBar(plotClusterInds);
    bar([1:length(plotBar)], plotBar, .8/8*(length(plotBar)))
    set(gca, 'xticklabel', ...
             gtexCluster.clusterLabels(plotClusterIDs));
    xtickangle(45)
    
    set(gca, 'FontSize', 16)
    mm = max(plotBar);
    for j = 1:length(plotClusterIDs)
        text(j-((.8/8*(length(plotBar)))/2), plotBar(j)+(mm*.05), ...
             sprintf('%.2f', meanR2(plotClusterIDs(j))), 'FontSize', ...
             16)
    end
    
    xlabel('clusterID - memberCount')
    ylabel('Marker count in cluster')
    title(sprintf('%s \n total count of markers: %d ', cellType{1}, ...
                  presentMarkerCounts(i)))
    ylim([0, mm + (mm*.1)])
    figFolder = ['~/resultsAndFigures/secondProject/markerRepInGTExClusters/']
    file = sprintf('%sSCbasedMarkers_%s_largeFont', figFolder, cellType{1});
    set(h, 'PaperOrientation', 'landscape')
    print(h, '-deps', [file '.eps'])
    print(h, '-dpdf', [file '.pdf'])
    saveas(h, [file '.eps'], 'epsc')
    % save all the plot Information
end

% >> for OganMarkers
markers
presentMarkerCounts = zeros(1, length(myCellTypes));
indCounter = 1
for i = [1:5, 8:10]
    cellType = markers.cellTypes(i);
    ctGeneSyms = markers.genes(markers.finalMarkerGenes{i});
    [a, b] = ismember(ctGeneSyms, gtexCluster.syms);
    presentMarkerCounts(indCounter) = sum(a); % can't use i
    presentClusters = gtexCluster.cs1547(b(a));
    presentClusterIDs = unique(presentClusters);
    % h = figure
    % hist(presentClusters, unique(presentClusterIDs));
    % title(sprintf('cellType: %s, totalMarkers: %d, presentMarkers: %d', ...
    %               scBasedMarkers.cellTypes{i}, length(ctGeneSyms), ...
    %               sum(a)))
    
    % I plot all clusters which have >= 5% of the markers in them 
    if length(unique(presentClusterIDs)) > 1
        myBar = hist(presentClusters, unique(presentClusterIDs));
        plotClusterInds = find(myBar >= (.05 * presentMarkerCounts(indCounter)));
        
        plotClusterIDs = presentClusterIDs(plotClusterInds); % Xtick
                                                             % label of
                                                             % the plot
    else 
        myBar = hist(presentClusters, [unique(presentClusterIDs), 300]);
        plotClusterInds = find(myBar >= (.05 * presentMarkerCounts(indCounter)));
        
        plotClusterIDs = presentClusterIDs(plotClusterInds); % Xtick
                                                             % label of
    end

    % text in plot the R2 of the cluster, 
    h = figure
    plotBar = myBar(plotClusterInds);
        mm = max(plotBar)
    bar([1:length(plotBar)], plotBar, .8/8*(length(plotBar)))
    set(gca, 'xticklabel', ...
             gtexCluster.clusterLabels(plotClusterIDs));
    xtickangle(45)
    
    for j = 1:length(plotClusterIDs)
        text(j-((.8/8*(length(plotBar)))/2), plotBar(j)+(mm*.05), sprintf('%.2f', ...
                                             meanR2(plotClusterIDs(j))), ...
             'FontSize', 16)
    end
    
    xlabel('clusterID - memberCount')
    ylabel('Marker count in cluster')
    title(sprintf('%s \n total count of markers: %d ', cellType{1}, ...
                  presentMarkerCounts(indCounter)))
    set(gca, 'FontSize', 16)
    ylim([0, mm + mm*.3])
    figFolder = ['~/resultsAndFigures/secondProject/markerRepInGTExClusters/']
    file = sprintf('%smouseBasedMarkers_%s_largeFont', figFolder, cellType{1});
    set(h, 'PaperOrientation', 'landscape')
    print(h, '-deps', [file '.eps'])
    print(h, '-dpdf', [file '.pdf'])
    saveas(h, [file '.eps'], 'epsc')
    % save all the plot Information
    
    indCounter = indCounter + 1;
end

% TODO 
% density of links in the GTEx 

% density of links in the  SC 
% 9. Density of the clusters in GTEx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the GTEx cluster and figure out which cluster is which 
load('~/resultsAndFigures/secondProject/gtexClusters.mat') % 'gtexCluster'

load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat') 
gtexNet = GTExFiveNets.nets(2).net005;
gtexFullNet = gtexNet + gtexNet';
clear gtexNet
gtexNetSyms = ...
    GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(2).expGenes);

tempInds = 1:253;
myInds = tempInds(find(gtexCluster.memberCounts>20));
clear tempInds

overlapDensity = zeros(69, 69);
for i = 1:69
    i
    geneSet1 = gtexCluster.syms(gtexCluster.cs1547 == myInds(i)); 
    [a, b1] = ismember(geneSet1, gtexNetSyms);
    overlapDensity(i, i) = sum(sum(gtexFullNet(b1, b1)))/(length(a) *(length(a)-1));
    for j = i+1:69
        geneSet2 = gtexCluster.syms(gtexCluster.cs1547 == myInds(j)); 
        [a, b2] = ismember(geneSet2, gtexNetSyms);
        chunk = gtexFullNet(b1, b2);
        overlapDensity(i, j) = sum(sum(chunk))/(length(b1) *(length(b2)));
    end
end

plotMat = triu(overlapDensity,1) + overlapDensity';
h = figure
heatmap(gtexCluster.clusterLabels(myInds), ...
        gtexCluster.clusterLabels(myInds), plotMat, 'FontSize', 8)
colormap(jet)
title('Overlap of the GTEx cluster with each other')
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sDensityOverlapOfClusters', figFolder, cellType{1});
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 9.5 writing the file for Cytoskape : cluster file and the node
% file 

% >>>>> cluster file: 
% for each cluster, I should have 1: the cluster ID, 2.the count of
% genes. 3. the intra cluster links (density) 4. the R2 

fileName = 'clusterInfo.txt'
file = ['~/resultsAndFigures/secondProject/cytoFilesForClusters/' fileName]
fid = fopen(file, 'w')
fprintf(fid, ['clusterID\t geneCount\t density \t R2\n'])
for i = 1:length(myInds)
    ID = myInds(i);
    gCount = gtexCluster.memberCounts(ID);
    density = plotMat(i, i);
    R2 = gtexCluster.meanRs(ID);
    fprintf(fid, '%d\t%d\t%.2f\t%.2f\n', ID, gCount, density, R2);
end
fclose(fid)

% >>>>> clusternetwork file
% for each two clusters, I need to have their link density - that's
% all.
fileName = 'interClusterInfo.txt'
file = ['~/resultsAndFigures/secondProject/cytoFilesForClusters/' fileName]
fid = fopen(file, 'w')
[a, b, c] = find(triu(plotMat, 1));
fprintf(fid, ['clusterID_1\t clusterID_2\t density\n'])
for i = 1:length(c)
    i
    ID1 = myInds(a(i));
    ID2 = myInds(b(i));
    d = c(i);
    fprintf(fid, '%d\t%d\t%.2f\n', ID1, ID2, d);
end
fclose(fid)

fileName = 'clusterInfo.txt'
file = ['~/resultsAndFigures/secondProject/cytoFilesForClusters/' fileName]
fid = fopen(file, 'w')
fprintf(fid, ['clusterID\t geneCount\t density \t R2\n'])
for i = 1:length(myInds)
    ID = myInds(i);
    gCount = gtexCluster.memberCounts(ID);
    density = plotMat(i, i);
    R2 = gtexCluster.meanRs(ID);
    fprintf(fid, '%d\t%d\t%.2f\t%.2f\n', ID, gCount, density, R2);
end
fclose(fid)

% color of clusters 
colors = zeros(1, 69);
[a, b] = ismember([98], myInds)
colors(b) = 1;
[a, b] = ismember([252], myInds)
colors(b) = 2;
[a, b] = ismember([65], myInds)
colors(b) = 3;
[a, b] = ismember([35 53 231 244], myInds)
colors(b(a)) = 4;
[a, b] = ismember([129, 239, 242, 243], myInds)
colors(b(a)) = 5;
[a, b] = ismember([64, 163, 168, 169, 189, 248, 253], myInds)
colors(b(a)) = 6;
[a, b] = ismember([52, 88, 105 107, 124, 136, 147, 164, 170, 190, ...
                   197, 198, 233 250], myInds)
colors(b) = 7;
[a, b] = ismember([40, 78, 83 113, 130, 159, 167, 234, 241], myInds)
colors(b) = 8;


fileName = 'clusterInfo_colors.txt'
file = ['~/resultsAndFigures/secondProject/cytoFilesForClusters/' fileName]
fid = fopen(file, 'w')
fprintf(fid, ['clusterID\t color\n'])
for i = 1:length(myInds)
    ID = myInds(i);
    fprintf(fid, '%d\t%d\n', ID, colors(i));
end
fclose(fid)

% R2 for border thing
borderWidth = zeros(1, 69);
myR2s = gtexCluster.meanRs(myInds);
printR2s = zeros(1,69);
printR2s(myR2s < .3) = 1;
printR2s(myRs >= 3) = 2;
printR2s(myRs >= 4.5) = 3;

fileName = 'clusterInfo_r2s.txt'
file = ['~/resultsAndFigures/secondProject/cytoFilesForClusters/' fileName]
fid = fopen(file, 'w')
fprintf(fid, ['clusterID\t borderWidth\n'])
for i = 1:length(myInds)
    ID = myInds(i);
    fprintf(fid, '%d\t%d\n', ID, printRs(i));
end
fclose(fid)


% 9.8 getting the correlation of R2 and density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myRs = gtexCluster.meanRs(myInds);
myDs = diag(plotMat);
h = figure
scatter(myRs, myDs)

% density of links in the  SC 
% 10. study of the microglia cluster 
% TODO: just get me the correlaion of the node degree in other
% tissues and write up about it. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort the marker genes and non marker genes
mouseMarkerMicro = markers.genes(markers.finalMarkerGenes{5});
scMarkerMicro = scBasedMarkers.genes(logical(scBasedMarkers.markerMat(:, 3)));

t = 2
brainGTExNet = GTExFiveNets.nets(t).net005;
thisNetSyms = ...
    GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(t).expGenes);
[ins, b] = ismember(geneList, thisNetSyms);
smallBGN = brainGTExNet(b(ins), b(ins));
fullsbgn = smallBGN + smallBGN';
fullsbgnNDs = sum(fullsbgn);

[sa, sb] = sort(fullsbgnNDs, 'descend');
sortedGeneList = geneList(sb);

% modify the sortedGeneList to have markerGenes first
[ma, mb] = ismember(sortedGeneList, mouseMarkerMicro);
sortedGeneListMMB = [sortedGeneList(ma), ...
                    sortedGeneList(~ma)]; % mouse marker based
MMBIns = [find(ma), find(~ma)]

[ma, mb] = ismember(sortedGeneList, scMarkerMicro);
sortedGeneListSCB = [sortedGeneList(ma), ...
                    sortedGeneList(~ma)]; % mouse marker based
SCBIns = [find(ma), find(~ma)]

h = figure
heatmap(sortedGeneList, sortedGeneList, fullsbgn(sb, sb))

h = figure
heatmap(sortedGeneListMMB, sortedGeneListMMB, fullsbgn(sb(MMBIns), ...
                                                  sb(MMBIns)))
h = figure
heatmap(sortedGeneListSCB, sortedGeneListSCB, fullsbgn(sb(SCBIns), sb(SCBIns)))
title('gtex brain Microglia')

t = 1
bloodGTExNet = GTExFiveNets.nets(t).net005;
thisNetSyms = ...
    GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(t).expGenes);
[ins, b] = ismember(sortedGeneList, thisNetSyms);
smallLGN = bloodGTExNet(b(ins), b(ins));
fullslgn = smallLGN + smallLGN';
h = figure
heatmap(fullslgn(sb(MMBIns), sb(MMBIns)))

t = 3
liverGTExNet = GTExFiveNets.nets(t).net005;
thisNetSyms = ...
    GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(t).expGenes);
[ins, b] = ismember(sortedGeneList, thisNetSyms);
sum(ins)
smallVGN = liverGTExNet(b(ins), b(ins));
fullvlgn = smallVGN + smallVGN';
h = figure
heatmap(fullvlgn(sb(MMBIns), sb(MMBIns)))

% 11. The big fat bar plot 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'};
affyExpMat = zeros(18494, 5);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    affyExpMat(:, t) = expGenesInd;
end

t = 1
tissue = 'blood';
load( ['~/networks/tissues/' tissue '/' ...
       'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
myNet = binNet;


% for TSN
FDR = '0012'
FC = '3'
load(sprintf(['~/resultsAndFigures/firstProject/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))

% for GTEx 
load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat') 

% ctc
load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'correctedBinNets_logCorrected_jusResiduals_scBasedMarkers' ...
      '.mat'])

load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exonAndIntron_V4_clusterLabels.mat'])

load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'scClusterRep.mat']) % scGTExClusterRep
scGTExClusterRep
load('~/resultsAndFigures/secondProject/gtexClusterRep/TANClusterRep_V2.mat')
GTExClusterRep_tan
load('~/resultsAndFigures/secondProject/gtexClusterRep/TSNClusterRep_V2.mat')
GTExClusterRep_tsn
load('~/resultsAndFigures/secondProject/gtexClusterRep/GTExBloodClusterRep_V2.mat')
GTExClusterRep_Gblood
load('~/resultsAndFigures/secondProject/gtexClusterRep/GTExLiverClusterRep_V2.mat')
GTExClusterRep_Gliver
load('~/resultsAndFigures/secondProject/gtexClusterRep/affyLiverClusterRep_V2.mat')
GTExClusterRep_tanLiver
load('~/resultsAndFigures/secondProject/gtexClusterRep/affyBloodClusterRep_V2.mat')
GTExClusterRep_tanBlood
load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'GTExCTCClusterRep_V2.mat'])
GTExClusterRep_ctc
load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'GTExitselfClusterRep.mat'])
GTExClusterRep_itself

% get the GTEx cluster and figure out which cluster is which 
load('~/resultsAndFigures/secondProject/gtexClusters.mat') %
                                                           % 'gtexCluster'
cs = gtexCluster.cs1547;
load('~/data/brainSingleCell/scGeneIDfier.mat')
geneSyms = sc.geneSyms;

% I need the list of clusters, the inside and the outside number
% list of selected clusters and their IDs, and their relevant cell
% type

myInds = GTExClusterRep_tanBlood.clusterInds;
myClusterIDs = [242 239 129 243, 125, 244 53 231 35 92, 65, 252]
[a, localInds] = ismember(myClusterIDs, myInds);
myCellTypes = {'Pyramidal', 'Pyramidal', 'Pyramidal','Pyramidal', ...
              'Inhibitory',...
              'Astrocyte', 'Astrocyte', 'Astrocyte', 'Astrocyte', ...
               'Astrocyte',... 
              'Microglia',...
              'Oligodendrocyte'}

% For the 8 networks of bulk: gtex *3 , TAN * 3, TSN, CTC 
% recort the total count of links in the network, to report the LC
% 1- intra cluster link count 
% 2- expected intra cluster link count basd on the density 

% >>>>>> GTEx brain gbc
% get the lcs 
gbcTotal = sum(GTExFiveNets.nets(2).net005(:))
plotTotalLC(1) = gbcTotal;
gbcLCs = sum(GTExClusterRep_itself.lc(localInds)) % **
plotgc(1) = sum(GTExClusterRep_itself.gc(localInds)) % **
plotlc(1) = gbcLCs;

% get the expected count
gbcNullExpected = ((GTExClusterRep_itself.gc .* (GTExClusterRep_itself.gc ...
                                                - 1))./2) .* .005;
gbcNullExpectedSum = sum(gbcNullExpected(localInds)) % **
plotlcNull(1) = gbcNullExpectedSum;

% >>>>>> GTEx liver gl
% get the lcs 
glTotal = sum(GTExFiveNets.nets(3).net005(:))
plotTotalLC(2) = glTotal;
glLCs = sum(GTExClusterRep_Gliver.lc(localInds)) % **
plotgc(2) = sum(GTExClusterRep_Gliver.gc(localInds))
plotlc(2) = glLCs ;

% get the expected count
glNullExpected = ((GTExClusterRep_Gliver.gc .* (GTExClusterRep_Gliver.gc ...
                                                - 1))./2) .* .005;
glNullExpectedSum = sum(glNullExpected(localInds)) % **
plotlcNull(2) = glNullExpectedSum;

% >>>>>> GTEx blood gbl
% get the bl 
gblTotal = sum(GTExFiveNets.nets(1).net005(:))
plotTotalLC(3) = gblTotal;
gblLCs = sum(GTExClusterRep_Gblood.lc(localInds)) % **
plotgc(3) = sum(GTExClusterRep_Gblood.gc(localInds))
plotlc(3) = gblLCs;

% get the expected count
gblNullExpected = ((GTExClusterRep_Gblood.gc .* (GTExClusterRep_Gblood.gc ...
                                                - 1))./2) .* .005;
gblNullExpectedSum = sum(gblNullExpected(localInds)) % **
plotlcNull(3) = gblNullExpectedSum;

% >>>>>> GTEx ctc
% get the bl 
ctcTotal = sum(ctc.net005(:))
plotTotalLC(4) = ctcTotal;
ctcLCs = sum(GTExClusterRep_ctc.lc(localInds)) % **
plotgc(4) = sum(GTExClusterRep_ctc.gc(localInds))
plotlc(4) = ctcLCs;

% get the expected count
ctcNullExpected = ((GTExClusterRep_ctc.gc .* (GTExClusterRep_ctc.gc ...
                                                - 1))./2) .* .005;
ctcNullExpectedSum = sum(ctcNullExpected(localInds)) % **
plotlcNull(4) = ctcNullExpectedSum;

% >>>>>> TAN blood
tissue = 'blood';
load( ['~/networks/tissues/' tissue '/' ...
       'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
myNet = binNet;

% get the bl 
tanblTotal = sum(myNet(:))
plotTotalLC(5) = tanblTotal;
tanblDens = tanblTotal / (sum(affyExpMat(:,1)) * (sum(affyExpMat(:, ...
                                                  1))-1)/2)
clear myNet
tanblLCs = sum(GTExClusterRep_tanBlood.lc(localInds)) % **
plotgc(5) = sum(GTExClusterRep_tanBlood.gc(localInds))
plotlc(5) = tanblLCs;

% get the expected count
tanblNullExpected = ((GTExClusterRep_tanBlood.gc .* (GTExClusterRep_tanBlood.gc ...
                                                - 1))./2) .* tanblDens;
tanblNullExpectedSum = sum(tanblNullExpected(localInds)) % **
plotlcNull(5) = tanblNullExpectedSum;

% >>>>>> TAN liver
tissue = 'liver';
load( ['~/networks/tissues/' tissue '/' ...
       'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
myNet = binNet;

% get the lv
tanlvTotal = sum(myNet(:))
plotTotalLC(6) = tanlvTotal;
tanlvDens = tanblTotal / (sum(affyExpMat(:,3)) * (sum(affyExpMat(:, ...
                                                  3))-1)/2)
clear myNet
tanlvLCs = sum(GTExClusterRep_tanLiver.lc(localInds)) % **
plotgc(6) = sum(GTExClusterRep_tanLiver.gc(localInds))
plotlc(6) = tanlvLCs;

% get the expected count
tanlvNullExpected = ((GTExClusterRep_tanLiver.gc .* (GTExClusterRep_tanLiver.gc ...
                                                - 1))./2) .* tanlvDens;
tanlvNullExpectedSum = sum(tanlvNullExpected(localInds)) % **
plotlcNull(6) = tanlvNullExpectedSum;

% >>>>>> TAN brain
tissue = 'brain';
load( ['~/networks/tissues/' tissue '/' ...
       'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
myNet = binNet;

% get the br
tanbrTotal = sum(myNet(:))
plotTotalLC(7) = tanbrTotal;
tanbrDens = tanbrTotal / (sum(affyExpMat(:,2)) * (sum(affyExpMat(:, 2))-1)/2)
clear myNet
tanbrLCs = sum(GTExClusterRep_tan.lc(localInds)) % **
plotlc(7) = tanbrLCs;

% get the expected count
tanbrNullExpected = ((GTExClusterRep_tan.gc .* (GTExClusterRep_tan.gc ...
                                                - 1))./2) .* tanbrDens;
tanbrNullExpectedSum = sum(tanbrNullExpected(localInds)) % **
plotlcNull(7) = tanbrNullExpectedSum;

% >>>>>> TSN brain
myNet = finalTable(2).wholeNet;

% get the br
tsnTotal = sum(myNet(:))
plotTotalLC(8) = tsnTotal;
tsnDens = tsnTotal / (sum(affyExpMat(:,2)) * (sum(affyExpMat(:, ...
                                                  2))-1)/2)
clear myNet
tsnLCs = sum(GTExClusterRep_tsn.lc(localInds)) % **
plotlc(8) = tsnLCs;

% get the expected count
tsnNullExpected = ((GTExClusterRep_tsn.gc .* (GTExClusterRep_tsn.gc ...
                                                - 1))./2) .* tsnDens;
tsnNullExpectedSum = sum(tsnNullExpected(localInds)) % **
plotlcNull(8) = tsnNullExpectedSum;

% For the three sets of sim networks 
% 1- intra cluster LC
% 2- expected intra cluster link count

files = dir(['~/resultsAndFigures/secondProject/' ...
             'SimBulkNetworksFromSC/withCounts/*newComb3*.mat'])

simDensity =  zeros(length(files), length(myInds));
simlc = zeros(length(files), length(myInds));
netFileNames = cell(1, length(files));

gCounts = zeros(1, 69);
scGeneCount = zeros(1, 69);
for j = 1:length(files)
    j
    load([files(j).folder '/' files(j).name])
    myBinNet = bulkFromSC.binNet005;
    %    geneSyms = bulkFromSC.geneSyms;
    for i = 1:length(myInds)
        i
        tempIDs = find(cs == myInds(i));
        length(tempIDs)
        
        myGenes = gtexCluster.syms(tempIDs);
        gCount = length(myGenes)
        gCounts(i) = gCount;
        
        [a, b] = ismember(myGenes, geneSyms);
        sum(a)
        scGeneCount(i) = sum(a);
        % simDensity
        simDensity(j, i) = (sum(sum(myBinNet(b(a), b(a))))) / (scGeneCount(i) * (scGeneCount(i) - ...
                                                          1));
        simlc(j, i) = (sum(sum(myBinNet(b(a), b(a)))));
    end
end
totalLC = sum(myBinNet(:))

h = figure
heatmap(simlc)

lowVar.lc = simlc;
lowVar.density = simDensity;
lowVar.gc = scGeneCount;
lowVar.info = 'newComb3'

highVar.lc = simlc;
highVar.density = simDensity;
highVar.gc = scGeneCount;
highVar.info = 'newComb4'

noVar.lc = simlc;
noVar.density = simDensity;
noVar.gc = scGeneCount;
noVar.info = 'noVar'

nullExpected = (scGeneCount) .*  (scGeneCount-1) .* .005;

simNetLCs{1} = noVar;
simNetLCs{2} = lowVar;
simNetLCs{3} = highVar;

save(['~/resultsAndFigures/secondProject/SimBulkNetworksFromSC/' ...
      'GTEx_clusterInfo.mat'] ,'simNetLCs')

load(['~/resultsAndFigures/secondProject/SimBulkNetworksFromSC/' ...
      'GTEx_clusterInfo.mat'])

noVar = simNetLCs{1};
lowVar = simNetLCs{2};
highVar = simNetLCs{3};

noVarSums = sum(noVar.lc(:, localInds)')
lowVarSums = sum(lowVar.lc(:, localInds)')
highVarSums = sum(highVar.lc(:, localInds)')

nullExpectedSum = sum(nullExpected(localInds))

sib = [noVarSums(1:15); lowVarSums(1:15); highVarSums(1:15)]

h = figure
boxplot((sib')./totalLC)
hold on
line([3.75, 4.25], ([nullExpectedSum, nullExpectedSum]./totalLC))
xlim([0 5])
ylim([0 .55])

title('Density')
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%slinkContributionSims', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% >>>>>>>>>>>> Final plot
% 5 >> for the eight networks
plotlc
plotlcNull
plotTotalLC

indOrder = [1 7 8 4 2 5 3 6];
h = figure
plot([1:8], plotlc(indOrder)./plotTotalLC(indOrder), 'o')
hold on
plot([1.2:1:8.2], plotlcNull(indOrder)./plotTotalLC(indOrder), 'o')
xlim([0 9])
ylim([0 .55])

title('Density')
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%slinkContributionThe8', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


% for the singleCell networks: 
% NOTE Densities are .01 and lcs should be divided by two (it is the
% full net not the triu)

% getting the sort for the celltypes
[a, sortSC] = sort(clusterMeta.sortedClusterNames);

% getting the lcs for the clusters
scGTExClusterRep.lc()

% 12. comparing the count of links in GTEx and ctc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scGTExClusterRep
load('~/resultsAndFigures/secondProject/gtexClusterRep/TANClusterRep_V2.mat')
GTExClusterRep_tan
load('~/resultsAndFigures/secondProject/gtexClusterRep/TSNClusterRep_V2.mat')
GTExClusterRep_tsn
load('~/resultsAndFigures/secondProject/gtexClusterRep/GTExBloodClusterRep_V2.mat')
GTExClusterRep_Gblood
load('~/resultsAndFigures/secondProject/gtexClusterRep/GTExLiverClusterRep_V2.mat')
GTExClusterRep_Gliver
load('~/resultsAndFigures/secondProject/gtexClusterRep/affyLiverClusterRep_V2.mat')
GTExClusterRep_tanLiver
load('~/resultsAndFigures/secondProject/gtexClusterRep/affyBloodClusterRep_V2.mat')
GTExClusterRep_tanBlood
load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'GTExCTCClusterRep_V2.mat'])
GTExClusterRep_ctc
load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'GTExitselfClusterRep.mat'])
GTExClusterRep_itself

load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'GTExCTCClusterRep.mat'])
GTExClusterRep_ctc
load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'GTExitselfClusterRep.mat'])
GTExClusterRep_itself

myInds = GTExClusterRep_itself.clusterInds;

% passedFDRL comes from 5.5
h = figure
heatmap([1, 2], gtexCluster.clusterLabels(myInds(passedFDRL>0)), ...
        ([GTExClusterRep_itself.lc(passedFDRL>0); ...
          GTExClusterRep_ctc.lc(passedFDRL>0)]'))


h = figure
heatmap([1:4], gtexCluster.clusterLabels(myInds(passedFDRL>0)), ...
        [GTExClusterRep_itself.lc(passedFDRL>0); ...
          GTExClusterRep_ctc.lc(passedFDRL>0); GTExClusterRep_Gblood.lc(passedFDRL>0); ...
         GTExClusterRep_Gliver.lc(passedFDRL>0)]')

percent = (GTExClusterRep_ctc.lc(passedFDRL>0) - GTExClusterRep_itself.lc(passedFDRL>0)) ...
    ./ GTExClusterRep_itself.lc(passedFDRL>0);


percent = (GTExClusterRep_ctc.lc - GTExClusterRep_itself.lc) ...
    ./ GTExClusterRep_itself.lc;

h = figure
heatmap(gtexCluster.clusterLabels(myInds), 1, percent*100)

h = figure('units', 'centimeters', 'position', [0,0, 33, 15])
heatmap(gtexCluster.clusterLabels(myInds), 1, percent*100)
colormap(1-gray)
title('percent of loss or gain of linkCount from GTEx_brain to GTEx_ctc')
figFolder = ['~/resultsAndFigures/secondProject/clusterProfiling/']
file = sprintf('%slossOrGainLCforGTExBrainToCTC', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 13. building the parts of the big plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'functionalEnrichment.mat']) % 'clusterFRes'
fCounts = length(clusterFRes.fTerms)
selectedPs = ((clusterFRes.ps >=0) + ((clusterFRes.ps* fCounts) <= .1)) == 2;
posF = sum(selectedPs') >= 3;

posF; % comes from line 1135, functional enrichment of the
      % clusters, it is the list of clusters with functions

% fixing the order for the HK and cell type specific modules
clusterCat = zeros(1, 69) + posF;
[a, b] = ismember([64, 163, 168, 169, 189, 248, 253], myInds)  % hk
clusterCat(b(a)) = 2; % HK CAT

[a, b] = ismember([129, 239, 243 242], myInds)  % pyra
clusterCat(b(a)) = 3; % pyra CAT

[a, b] = ismember([35, 53 43, 231, 244], myInds)  % astro
clusterCat(b(a)) = 4; % astro CAT

[a, b] = ismember([65], myInds)  % micro
clusterCat(b(a)) = 5; % micro cat

[a, b] = ismember(252, myInds)  % oligo
clusterCat(b(a)) = 6; % oligo CAT

[a, b] = ismember(98, myInds)  % endo
clusterCat(b(a)) = 7; % endo CAT

[as, bs] = sort(clusterCat, 'descend')

h = figure
heatmap(allPs.clusterLabels(bs), 1, clusterCat(bs), 'GridVisible', 'off')
myCmat = [160, 160, 160;
          219, 146, 117;
          221, 103, 51;
          140 191 136;
          32, 102, 242;
          125, 46, 141;
          237, 230, 105;
          76, 189, 237
         ]./256;
colormap(myCmat)

title('clusterCategories')
figFolder = ['~/resultsAndFigures/secondProject/clusterProfiling/']
file = sprintf('%sclusterCategoriesOrdered', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% p1: replot the allps 

load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'repPValues_76Nets_structure.mat']);
h = figure('units', 'centimeters', 'position', [0,0, 25, 22])
heatmap(allPs.clusterLabels(bs), allPs.netLabels, allPs.mat(:, bs), 'FontSize', ...
        7)
colormap((1-gray*.5))
title(['p-values for the reproducibility of GTEx-brain-cortex clusters ' ...
       'in other networks - ordered'])
figFolder = ['~/resultsAndFigures/secondProject/clusterProfiling/']
file = sprintf('%sPsrefined_ordered', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% p2: plot the allps for non SC
h = figure('units', 'centimeters', 'position', [0,0, 25, 5])
heatmap(allPs.clusterLabels(bs), allPs.netLabels(1:7),((allPs.mat(1:7, bs).*69)<=.1)+0, 'FontSize', ...
        8)
colormap((1-gray*.5))
title(['p-values for the reproducibility of GTEx-brain-cortex clusters ' ...
       'in other networks_ordered'])
figFolder = ['~/resultsAndFigures/secondProject/clusterProfiling/']
file = sprintf('%sPsrefined_nonSC_FDR.1_ordered', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% p3: plot the count SC nets which have each cluster, with pref on
% each of the four cell type: bar for Inh and Exc, The two Asc, the endothelial
scPs = allPs.mat(8:end, :);
scLabels = allPs.netLabels(8:end);
[a, b] = sort(scLabels);

kado = ~isnan(scPs);
notNANCount = sum(kado');

scAstro = scPs(b(1:2), :);
scAstro(1, :) = scAstro(1,:) .* notNANCount(b(1));
scAstro(2, :) = scAstro(2,:) .* notNANCount(b(2));
scOligo = scPs((b(end)), :) .* notNANCount(b(end));
scOPC = scPs(b(68), :) .* notNANCount(b(68));

draftPM = [scAstro; scOligo];
plotMat = (draftPM < .1) +0;
plotLabs = a([1:2, 69]);
h = figure('units', 'centimeters', 'position', [0,0, 25, 4])
heatmap(allPs.clusterLabels(bs), plotLabs, plotMat(:, bs), 'FontSize', ...
        8)
colormap(1-gray*.5)
title(['p-values for the reproducibility of GTEx-brain-cortex clusters ' ...
       'in other networks - ordered'])
figFolder = ['~/resultsAndFigures/secondProject/clusterProfiling/']
file = sprintf('%sPsrefined_AstroOligo_FDR.1_ordered', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

inhInds = b(3:25);
scInhs = zeros(size(inhInds,2), 69);
for i = 1:size(inhInds, 2)
    scInhs(i, :) = scPs(inhInds(i), :) .* notNANCount(inhInds(i));
end
h = figure
heatmap(scInhs)
inhBar = nansum(scInhs <= .1);

excInds = b(26:68);
scExc = zeros(size(excInds,2), 69);
for i = 1:size(excInds, 2)
    scExc(i, :) = scPs(excInds(i), :) .* notNANCount(excInds(i));
end

h = figure
heatmap(scExc)
excBar = nansum(scExc <= .1);

h = figure('units', 'centimeters', 'position', [0,0, 25, 5])
plotMat= [inhBar; excBar];
bar(plotMat(:, bs)', 'stacked')
set(gca, 'XTick', [1:69])
set(gca, 'XTickLabels', allPs.clusterLabels(bs))
xtickangle(90)
legend({'inh', 'exc'})
title(['Reproducibility of the GTEx-brainCortex clusters in Inh and ' ...
       'exc cell types - ordered'])
colormap(gray)
figFolder = ['~/resultsAndFigures/secondProject/clusterProfiling/']
file = sprintf('%sexcInhRep_barplot_ordered', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% p4 the bar colored for clusters identidified for each cell type
% and common

% TODO: >>>>> this fat thing remains
h = figure
book = heatmap((1:69))
title(['just the Cluster info'])
figFolder = ['~/resultsAndFigures/secondProject/clusterProfiling/']
file = sprintf('%scolormapTempForClusterInfo', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% P4.5  the R2 plot
h = figure
heatmap(allPs.clusterLabels(bs) , 1,(gtexCluster.meanRs(myInds(bs))))
colormap((1-(bone*.7)))
title(['meanRs'])
figFolder = ['~/resultsAndFigures/secondProject/clusterProfiling/']
file = sprintf('%smeanRsHeatmapForEachCluster_ordered', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% p5 : on the CTC: 
% the correlaion of the pvalues between the networks. CTC is
% negatively correlated with TAN brain and TSN, while more or less
% correlated with other values
h = figure('units', 'centimeters', 'position', [0,0, 25, 22])
h = figure
sib = corr(allPs.mat(1:7, :)', 'rows', 'pairwise')
heatmap(allPs.netLabels(1:7), allPs.netLabels(1:7), sib, 'FontSize', ...
        10)
colormap((1-gray*.5))
title(['correlation of the cluster Ps in different networks'])
figFolder = ['~/resultsAndFigures/secondProject/clusterProfiling/']
file = sprintf('%sclusterPsCorrelation', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


% the stem plot for increas or decrease of the links in the clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('~/resultsAndFigures/secondProject/gtexClusterRep/GTExCTCClusterRep_SC7PCA_redo.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
posF; % comes from line 1135, functional enrichment of the
      % clusters, it is the list of clusters with functions

% fixing the order for the HK and cell type specific modules
clusterCat = zeros(1, 69) + posF;
[a, b] = ismember([64, 163, 168, 169, 189, 248, 253], myInds)  % hk
clusterCat(b(a)) = 2; % HK CAT

[a, b] = ismember([129, 239, 243 242], myInds)  % pyra
clusterCat(b(a)) = 3; % pyra CAT

[a, b] = ismember([35, 53 43, 231, 244], myInds)  % astro
clusterCat(b(a)) = 4; % astro CAT

[a, b] = ismember([65], myInds)  % micro
clusterCat(b(a)) = 5; % micro cat

[a, b] = ismember(252, myInds)  % oligo
clusterCat(b(a)) = 6; % oligo CAT

[a, b] = ismember(98, myInds)  % endo
clusterCat(b(a)) = 7; % endo CAT

[as, bs] = sort(clusterCat, 'descend')

% load GTExClusterRep_ctc and GTExClusterRep_itself

% now give me, for each cluster, what has happened to it : bs has
% the orders
% get the BS from above or get it from the clustering plot
bst = [244 126 231 57  84 45 107 35 53 92 161 250,... % astrocyte
     242 239 243 87 136 232 125 174 180 129 169 124 39 123 236 168,...
      170 246, 64, ... %pyramidal
      252 98 52 132 105 164 115 23 217 230 233 147 190,... % olygo
                                                           % dendro
      65 241 167 78 83 130 248 234 88 40 113 159 163 197, ...
      189, 253]

[a, bs] = ismember(bst, myInds);
myClusterIDs = myInds(bs);

lcDiff = ((GTExClusterRep_ctc.lc(bs) - GTExClusterRep_itself.lc(bs)) ...
         ./ GTExClusterRep_itself.lc(bs)) * 100;

h = figure
h = figure('units', 'centimeters', 'position', [0,0, 30, 8])
stem(lcDiff)
ylim([-120, 500])
xlim([0, 61])
xticks([1:60])
xticklabels(gtexCluster.clusterLabels(myInds(bs)));
xtickangle(90)
title('GTEx to CTC percent of links')
figFolder = ['~/resultsAndFigures/secondProject/clusterProfiling/']
file = sprintf('%sstemplot03_resize01_redo_Fixed', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% the R2 for the stem plot

load('~/resultsAndFigures/secondProject/gtexClusters.mat') %
load('~/data/GTEx/Brain_Cortex_expGenes.mat')
load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_AllenSCBasedMarker_redo.mat'])

clusterR2s = zeros(1, 253);
for i = 1:253
    inGenes = gtexCluster.syms(gtexCluster.cs1547 == i);
    [a, b] = ismember(dataSet.genes, inGenes);
    clusterR2s(i) = mean(result.r2(a));
end

%sR2 = gtexCluster.meanRs(myInds(bs));
%sR2 = myR2s(bs);
sR2 = clusterR2s(myClusterIDs)
h = figure
heatmap(gtexCluster.clusterLabels(myInds(bs)), 1, sR2)
title('r2')
colormap(1-bone)
figFolder = ['~/resultsAndFigures/secondProject/clusterProfiling/']
file = sprintf('%sstemplot03_r2heatmap_redo_Fixed', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% the cluster color for the stem plot
sColor = colors(bs);
h = figure
heatmap(gtexCluster.clusterLabels(myInds(bs)), 1, sColor)
colormap(jet)
title('the colors of clusters')
figFolder = ['~/resultsAndFigures/secondProject/clusterProfiling/']
file = sprintf('%sstemplot02_clusterColor_Fixed', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


% 14. getting the joint function of the genes in some clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following cluster IDs are connected 
sids = [40, 78, 83, 113 130 159 171 234 241]
sids = []
[a, b] = ismember(gtexCluster.cs1547, sids);
myGenes = gtexCluster.syms(a);
myCs = gtexCluster.cs1547(a);
[as, bs] = sort(myCs);
[fdrID, fdrTerms, inGenes, enrichMat, passedFDR] = ...
    fenrich_function(myGenes, .1);

h = figure
heatmap(enrichMat(bs, :)', 'GridVisible', 'off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat') 
        firstNet = GTExFiveNets.nets(2).net005;
        firstSyms = ...
            GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(t).expGenes);

inGenes = sum(firstNet + firstNet')>0;
inGenesSyms = dataSet.genes(inGenes);
inGenesR = result.r2(inGenes);
[a, b] = ismember(hkgSyms, inGenesSyms);
hkRs = inGenesR(b(a));

nonInds = ones(sum(inGenes), 1);
nonInds(b(a)) = 0;

otherRs = inGenesR(logical(nonInds));

[x, y] = ksdensity(inGenesR);
[x1, y2] = ksdensity(myRs);

h = figure
plot(y, x)
hold on
plot(y2, x1)

h = figure
sib = [inGenesR', myRs'];
kado = [zeros(1, length(inGenesR)), ones(1, length(myRs))];
boxplot(sib, kado)

figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sR2forHKversusAllInNetBox', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

sum(inGenesR < .2)
sum(myRs < .2)

