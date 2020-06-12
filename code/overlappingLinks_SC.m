% here I want to study the overlapping clusters in SC, finding
% possible missing interpretations in the GTEx dataset. 
% 1. general plots for the overlapping links in SC
% 2. looking at the cluster individually
% 3. implication of the functional enrichment between clusters 9
% and 10
% 4. they are sorted based on GTEx order now... I will pull out the
% cluster 9 and cluster 10. 


clear
% this file is made in gettingTheTSSDistForBiPrimer.m
load(['~/resultsAndFigures/secondProject/' ...
      'overlappingLinksBetweenSCNets/highRepLinks.mat'])
overlapRes

% load regression from mouse
load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_MarkersV1_cthr6.mat'])

resGTEx = result;
clear result

% load regression from SC
load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_AllenSCBasedMarker.mat'])
result
resSC = result;

load('~/resultsAndFigures/secondProject/gtexClusters.mat')
load('~/data/GTEx/Brain_Cortex_expGenes.mat')

h = figure
heatmap(log2(overlapRes.gtexBinNet01 + 1), 'GridVisible', 'off')

h = figure
book = hist(result.clusters, 10)

% 1. general plots for the overlapping links in SC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first important step:
% cluster inhSum and excSum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recluster the inh and exc 
% take it to clusterMyGenes.m

overlapRes.gcFilter = logical(ones(1,530)); % gene clusterfilter
overlapRes.gcFilter(23:62) = 0;

[x, inhOrder] = sort(overlapRes.inhClusters);
overlapRes.inhOrder = inhOrder;

[x, excOrder] = sort(overlapRes.excClusters);
overlapRes.excOrder = excOrder;

h = figure
heatmap(overlapRes.excSum(overlapRes.gcFilter, overlapRes.gcFilter), 'GridVisible', 'off')
title('exc')
h = figure
heatmap(overlapRes.inhSum(inhOrder, inhOrder), 'GridVisible', 'off')
title('inh')

% do the two big clusters of theirs overlap and if so, how

inds9 = [min(find(overlapRes.clusters == 9)):max(find(overlapRes.clusters ...
                                                  == 9))];
inds10 = [min(find(overlapRes.clusters == 10)):max(find(overlapRes.clusters ...
                                                  == 10))];

h = figure 
hist(overlapRes.inhClusters(inds9), ...
     unique(overlapRes.inhClusters(inds9)))

h = figure 
hist(overlapRes.excClusters(inds9), ...
     unique(overlapRes.excClusters(inds9)))

% plotting the figures: 

% the all sum figure
filter = overlapRes.gcFilter;
h = figure
heatmap(log2(overlapRes.sumNet(filter, filter)), 'GridVisible', 'off')
colormap(1-gray)
title('Clusters identified in SC sum network - gene count: 490 ')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%ssumNetHeatmap_gcFiltered', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% the inh sum 
h = figure
heatmap( log2(overlapRes.inhSum(filter, filter)), 'GridVisible', 'off')
colormap(1-gray)
title('Clusters identified in SC sum network - inh - gene count: 490 ')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%ssumNetHeatmap_inh_gcFiltered', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% the exci sum
h = figure
heatmap( log2(overlapRes.excSum(filter, filter)), 'GridVisible', 'off')
colormap(1-gray)
title('Clusters identified in SC sum network - exc - gene count: 490 ')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%ssumNetHeatmap_exc_gcFiltered', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% gtex net
h = figure
heatmap((overlapRes.gtexBinNet005(filter, filter)), 'GridVisible', 'off')
colormap(1-gray)
title('Clusters identified in SC sum network - gtexNetwork - gene count: 490 ')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%ssumNetHeatmap_gtex005_gcFiltered', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% all sum, each cluster sorted by GTEx
overlapRes
[a, b] = ismember(overlapRes.gtexGeneSyms(filter), gtexCluster.syms);
h = figure
book = overlapRes.gtexBinNet005(filter, filter);
heatmap((book(a,a)), 'GridVisible', 'off')
colormap(1-gray)

gclusterInds = a + 0;
gclusterInds(a) = gtexCluster.cs1547(b(a));
uniqueIDs = unique(gclusterInds);
book = hist(gclusterInds, uniqueIDs);
zerosIDs = uniqueIDs(book <5);

[a, b] = ismember(gclusterInds, zerosIDs);
gclusterInds(a) = 0;

[sa, sb] = sort(gclusterInds);

h = figure
%%% TAKE THIS 
sumFiltered = overlapRes.sumNet(filter, filter);
heatmap(log2(sumFiltered(sb, sb)), 'GridVisible', 'off')
colormap(1-gray)

title('Clusters identified in SC sum network - gtexOrdered ')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%ssumNetHeatmap_gtexOrdered_gcFiltered', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


% >>> sorting the gtex clusters based on the type (rather than the
% simple grouping)

template = gclusterInds;
unique(template)
order = [57, 252, 129 242 243 239, 64 168 253 163 169, 88, 203, 78]

for i = 1:length(order)
    template(template == order(i)) = i;
end

[sa, sb] = sort(template);

h = figure 
% I replaced overlapRes.sumNet with sumFiltered
heatmap(log2(sumFiltered(sb, sb)), 'GridVisible', 'off')
colormap(1-gray)
title('Clusters identified in SC sum network - gtexOrdered ')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%ssumNetHeatmap_gtexOrdered_grouped_gcFiltered', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% gtex sorted
gFiltered = overlapRes.gtexBinNet005(filter, filter);
h = figure
heatmap(gFiltered(sb, sb), 'GridVisible', 'off')
colormap(1-gray)


title('Clusters identified in SC sum network - gtexOrdered ')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%ssumNetHeatmap_gtexOrdered_gropued_gcFiltered', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% print the color bar 
h = figure
subplot(2, 1, 1)
heatmap(sa, 'GridVisible', 'off')
myCmat = [160, 160, 160;
          32, 102, 242;          
           237, 230, 105;
           140 191 136;
           140 191 136;
           140 191 136;
           140 191 136;
           221, 103, 51;
           221, 103, 51;
           221, 103, 51;
           221, 103, 51;
           221, 103, 51;
           219, 146, 117;
          130, 130, 130;
          160, 160, 160;
          ]./256;
colormap(myCmat)
order = [57, 252, 129 242 243 239, 64 168 253 163 169, 88, 203, 78]

mybar = hist(sa, unique(sa));
subplot(2, 1, 2)
bar([mybar; rand(1, 14)], 'Stacked')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%ssumNetHeatmap_gtexOrdered_grouped_clorbar_gcFiltered', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% get the color bar for the cluster IDs
h = figure
heatmap(sa, 'GridVisible', 'off')
colormap(lines(14))
h= figure
mybar = hist(sa, unique(sa))
bar([mybar; rand(1,14)], 'Stacked')
tl = unique(sa);
tlf = num2str(tl)
book = strsplit(tlf)
legend(book)
title('clusters ordered in GTEx ')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%sdraftPlot_clusterOrdersInGTEx_gcFiltered', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% get the color bar for the cluster IDs
myColors = [255, 127, 0; 
           255, 255, 51;
           228, 26, 28;
           247, 129, 191;
           255, 255, 51;
           255, 255, 51;           
           255, 255, 51;           
           255, 255, 51;           
           135, 135, 135;
            50, 50, 50;
           ]./255
filClusters = overlapRes.clusters(filter);
kado = filClusters(sb);
h = figure
heatmap(kado, 'GridVisible', 'off')
colormap(myColors)
title('clusters ordered in GTEx - SC labels')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%sdraftPlot_clusterOrdersInGTEx_SCLabels_grouped_gcFiltered', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% all sum, with no var 

% all sum, with high var

% 239, GTEX brain
 h = figure
heatmap(gFiltered(sb(sa == 6), sb(sa == 6)), 'GridVisible', 'off')
colormap(1 - gray)
title('clusters ordered in GTEx - SC labels')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%sgtexHeatmap_cluster239', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 239 SC no var

% 239, SC var

% 239, sc sum

% all sum, with GTEx colors

% all sum, GTEx sorted and colors

% 2. looking at the cluster individually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cluster 1 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% can we explain why they are correlated? 
% mitochondrial genes most of them
geneList = overlapRes.gtexGeneSyms(1:18)
geneList = overlapRes.geneSyms(1:18);
[a, b] = ismember(geneList, dataSet.genes);
gtexDSInds = b(a);
h = figure
plot(resSC.r2(b(a)), 'o')
hold on
plot(resGTEx.r2(b(a)), 'o')
legend('sc', 'gtex')
set(gca, 'XTick', [1:18])
set(gca, 'XTickLabels', geneList(a))
xtickangle(60)

resGTEx.coeff(gtexDSInds([13, 14]), :)

h = figure
heatmap(geneList, geneList, overlapRes.excSum(1:18, 1:18))
title('inhibitory')
h = figure
heatmap(geneList, geneList, overlapRes.inhSum(1:18, 1:18))
title('excitatory')
h = figure
heatmap(geneList, geneList, overlapRes.gtexBinNet01(1:18, 1:18))
title('GTEx')
h = figure
heatmap(geneList, geneList, overlapRes.stringExperimental(1:18, 1:18))
% are they similarly correlated in the bulk tissues? (same cluster)

% I got the expression levels and coexpression of CALM3 and NRGN in
% plots
% get teh box plot of expraession for CALM3 and NRGN
[a, b] = ismember('NRGN', dataSet.genes)
e1 = dataSet.mat(b,:);

[a, b] = ismember('CALM3', dataSet.genes)
e2 = dataSet.mat(b,:);

h = figure
boxplot(log2([e1' e2']))
ylim([7, 12])
set(gca, 'XTickLabels', {'NRGN', 'CALM3'})
title('expression levels in GTEx')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%sgtexExpression_NRGN_CALM3', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% cluster 2 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% can we explain why they are correlated? 
% mitochondrial genes most of them
geneList = overlapRes.gtexGeneSyms(1:18)
geneList = overlapRes.geneSyms(1:18);
[a, b] = ismember(geneList, dataSet.genes);
gtexDSInds = b(a);
h = figure
plot(resSC.r2(b(a)), 'o')
hold on
plot(resGTEx.r2(b(a)), 'o')
legend('sc', 'gtex')
set(gca, 'XTick', [1:18])
set(gca, 'XTickLabels', geneList(a))
xtickangle(60)

resGTEx.coeff(gtexDSInds([13, 14]), :)

h = figure
heatmap(geneList, geneList, overlapRes.excSum(1:18, 1:18))
title('inhibitory')
h = figure
heatmap(geneList, geneList, overlapRes.inhSum(1:18, 1:18))
title('excitatory')
h = figure
heatmap(geneList, geneList, overlapRes.gtexBinNet01(1:18, 1:18))
title('GTEx')
h = figure
heatmap(geneList, geneList, overlapRes.stringExperimental(1:18, 1:18))

% are they similarly correlated in the bulk tissues? (same cluster)

% if not, what is happening, why are some of them in a whole other
% cluster? 

% cluster 3 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% can we explain why they are correlated? 
% mitochondrial genes most of them
clusterIns = overlapRes.clusters == 3;
ggeneList = overlapRes.gtexGeneSyms(clusterIns)
scgeneList = overlapRes.geneSyms(clusterIns);
[a, b] = ismember(ggeneList, dataSet.genes);
gtexDSInds = b(a);
h = figure
plot(resSC.r2(b(a)), 'o')
hold on
plot(resGTEx.r2(b(a)), 'o')
legend('sc', 'gtex')
set(gca, 'XTick', [1:18])
set(gca, 'XTickLabels', ggeneList(a))
xtickangle(60)

resGTEx.coeff(gtexDSInds([13, 14]), :)

h = figure
heatmap(scgeneList, scgeneList, overlapRes.excSum(clusterIns, clusterIns))
title('inhibitory')
h = figure
heatmap(scgeneList, scgeneList, overlapRes.inhSum(clusterIns, clusterIns))
title('excitatory')
h = figure
heatmap(scgeneList, scgeneList, overlapRes.gtexBinNet01(clusterIns, ...
                                                  clusterIns))
title('GTEx')
h = figure
heatmap(scgeneList, scgeneList, overlapRes.stringExperimental(clusterIns, ...
                                                  clusterIns))

% are they similarly correlated in the bulk tissues? (same cluster)

% if not, what is happening, why are some of them in a whole other
% cluster? 

% cluster 4 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% can we explain why they are correlated? 
% can we explain why they are correlated? 
% mitochondrial genes most of them
clusterIns = overlapRes.clusters == 4;
ggeneList = overlapRes.gtexGeneSyms(clusterIns)
scgeneList = overlapRes.geneSyms(clusterIns);
[a, b] = ismember(ggeneList, dataSet.genes);
gtexDSInds = b(a);
h = figure
plot(resSC.r2(b(a)), 'o')
hold on
plot(resGTEx.r2(b(a)), 'o')
legend('sc', 'gtex')
set(gca, 'XTick', [1:18])
set(gca, 'XTickLabels', ggeneList(a))
xtickangle(60)

resGTEx.coeff(gtexDSInds([13, 14]), :)

h = figure
heatmap(scgeneList, scgeneList, overlapRes.excSum(clusterIns, clusterIns))
title('inhibitory')
h = figure
heatmap(scgeneList, scgeneList, overlapRes.inhSum(clusterIns, clusterIns))
title('excitatory')
h = figure
heatmap(scgeneList, scgeneList, overlapRes.gtexBinNet01(clusterIns, ...
                                                  clusterIns))
title('GTEx')
h = figure
heatmap(scgeneList, scgeneList, overlapRes.stringExperimental(clusterIns, ...
                                                  clusterIns))

% are they similarly correlated in the bulk tissues? (same cluster)

% if not, what is happening, why are some of them in a whole other
% cluster? 

% cluster 9
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% can we explain why they are correlated? 

clusterIns = overlapRes.clusters == 9;
ggeneList9 = overlapRes.gtexGeneSyms(clusterIns);
scgeneList9 = overlapRes.geneSyms(clusterIns);

% get the gtex cluster ids
[a, b] = ismember(ggeneList9, gtexCluster.syms);
gtexClusters9 = a + 0;
gtexClusters9(a) = gtexCluster.cs1547(b(a));

% get the member counts in each cluster
unGC9 = unique(gtexClusters9);
mc9 = hist(gtexClusters9, unique(gtexClusters9));
unGC9(mc9<5) = 0;

inIDs9 = unGC9(unGC9 > 0);

% sort them based on the GTEx cluster
counts9 = zeros(1, length(inIDs9));
thisSyms = gtexCluster.syms(gtexCluster.cs1547 == inIDs9(1));
[a, b] = ismember(ggeneList9, thisSyms);
counts9(1) = sum(a);
mySyms9 = ggeneList9(a);
for i = 2:length(inIDs)
    thisSyms = gtexCluster.syms(gtexCluster.cs1547 == inIDs9(i));
    [a, b] = ismember(ggeneList9, thisSyms);
    counts9(i) = sum(a);
    mySyms9 = [mySyms9 ggeneList9(a)];
end

% now mySyms9 are ordered: some of the genes are missing, since
% they don't belong to a cluster or etc

% now get their R2 values
[a, b] = ismember(mySyms9, dataSet.genes);
myR2s9 = resSC.r2(b(a));
h = figure
plot(myR2s9, 'o')

% boxplot grouping
groups9 = zeros(1, sum(counts9));
s = 1
for i = 1:length(counts9)
    e = s+counts9(i)-1;
    groups9(s:e) = inIDs9(i);
    s = e+1;
end

boxplot(myR2s9, groups9)

% >>>> get the SC symbols again
[a, b] = ismember(mySyms9, overlapRes.gtexGeneSyms);
mySyms9sc = overlapRes.geneSyms(b(a));

% get their avrg expression thing

% for R2
[a, b] = ismember(ggeneList9, dataSet.genes);
gtexDSInds = b(a);
h = figure
hist(resSC.r2(b(a)), 30)
plot(resSC.r2(b(a)), 'o')
hold on
plot(resGTEx.r2(b(a)), 'o')
legend('sc', 'gtex')
set(gca, 'XTick', [1:18])
set(gca, 'XTickLabels', ggeneList(a))
xtickangle(60)

% for cluster
[a, b] = ismember(ggeneList9, gtexCluster.syms);
gClusters9 = gtexCluster.cs1547(b(a));
allList9Clusters = a + 0;
allList9Clusters(a) = gtexCluster.cs1547(b(a));
exportGeneList = ggeneList9(allList9Clusters == 64);
h = figure
cIDs = unique(gClusters9);
cCounts = hist(gClusters9, unique(gClusters9))
inIDs = cCounts > 10;
idsAndCounts10 = [cIDs(inIDs)'; cCounts(inIDs)]

resGTEx.coeff(gtexDSInds([13, 14]), :)

h = figure
heatmap(scgeneList, scgeneList, overlapRes.excSum(clusterIns, clusterIns))
title('inhibitory')
h = figure
heatmap(scgeneList, scgeneList, overlapRes.inhSum(clusterIns, clusterIns))
title('excitatory')
h = figure
heatmap(scgeneList, scgeneList, overlapRes.gtexBinNet01(clusterIns, ...
                                                  clusterIns))
title('GTEx')
h = figure
heatmap(scgeneList, scgeneList, overlapRes.stringExperimental(clusterIns, ...
                                                  clusterIns))

% are they similarly correlated in the bulk tissues? (same cluster)

load('~/data/GTEx/GTExGeneIDfier.mat')

% get the NCBI IDs 

[a, b] = ismember(ggeneList9, gtexGenes.symbols);

myIDs = gtexGenes.NCBI(b(a));
book(1:10)

fileName = ['~/resultsAndFigures/secondProject/cluster9_SCOverlappingClusters_NCBI.txt']
fid = fopen(fileName, 'w')
for j = 1:length(myIDs)
    thisID = myIDs{j}
    if (length(thisID) > 0 && length(thisID) < 6)
        diff = 6 - length(thisID);
        for i = 1:diff
            thisID = ['0' thisID];
        end
        thisID = ['NM_', thisID]
        fprintf(fid, sprintf('%s\n', thisID))
    end
end
fclose(fid)

% in the bulk tissue, they are divided into other clusters. why is
% it happening 

% cluster 10
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% can we explain why they are correlated? 
clusterIns = overlapRes.clusters == 10;
ggeneList10 = overlapRes.gtexGeneSyms(clusterIns)
scgeneList10 = overlapRes.geneSyms(clusterIns);

% get the gtex cluster ids
[a, b] = ismember(ggeneList10, gtexCluster.syms);
gtexClusters10 = a + 0;
gtexClusters10(a) = gtexCluster.cs1547(b(a));

% get the member counts in each cluster
unGC10 = unique(gtexClusters10);
mc10 = hist(gtexClusters10, unique(gtexClusters10));
unGC10(mc10<5) = 0;

inIDs10 = unGC10(unGC10 > 0);

% sort them based on the GTEx cluster
counts10 = zeros(1, length(inIDs10));
thisSyms = gtexCluster.syms(gtexCluster.cs1547 == inIDs10(1));
[a, b] = ismember(ggeneList10, thisSyms);
counts10(1) = sum(a);
mySyms10 = ggeneList10(a);
for i = 2:length(inIDs10)
    thisSyms = gtexCluster.syms(gtexCluster.cs1547 == inIDs10(i));
    [a, b] = ismember(ggeneList10, thisSyms);
    counts10(i) = sum(a);
    mySyms10 = [mySyms10 ggeneList10(a)];
end

% now mySyms10 are ordered: some of the genes are missing, since
% they don't belong to a cluster or etc

% now get their R2 values
[a, b] = ismember(mySyms10, dataSet.genes);
myR2s10 = resSC.r2(b(a));
h = figure
plot(myR2s, 'o')

% boxplot grouping
groups10 = zeros(1, sum(counts10));
s = 1
for i = 1:length(counts10)
    e = s+counts10(i)-1;
    groups10(s:e) = inIDs10(i);
    s = e+1;
end

boxplot(myR2s10, groups10)

% >>>> get the SC symbols again
[a, b] = ismember(mySyms10, overlapRes.gtexGeneSyms);
mySyms10sc = overlapRes.geneSyms(b(a));

% sort them based on the GTEx cluster
thisSyms = gtexCluster.syms(gtexCluster.cs1547 == 78);
[a, b] = ismember(ggeneList10, thisSyms);
sum(a)
mySyms10 = ggeneList10(a);
mySyms10 = [mySyms10 ggeneList10(a)];

% for R2
[a, b] = ismember(ggeneList, dataSet.genes);
gtexDSInds = b(a);
h = figure
plot(resSC.r2(b(a)), 'o')
hold on
plot(resGTEx.r2(b(a)), 'o')
legend('sc', 'gtex')
set(gca, 'XTick', [1:18])
set(gca, 'XTickLabels', ggeneList(a))
xtickangle(60)

% for cluster
[a, b] = ismember(ggeneList10, gtexCluster.syms);
gClusters10 = gtexCluster.cs1547(b(a));
h = figure
cIDs = unique(gClusters);
cCounts = hist(gClusters, unique(gClusters))
inIDs = cCounts >= 5;
idsAndCounts10 = [cIDs(inIDs)'; cCounts(inIDs)]

% get the clusterIDs for sorting the expression patterns 
clusterIDs = zeros(size(a));
clusterIDs(a) = gClusters;
% now get the genes sorted based on the cluster IDs 
[cidsa, cidsb] = sort(clusterIDs)
sortscGeneList = scgeneList(cidsb);

[a, b] = ismember(ggeneList10, gtexGenes.symbols);

myIDs = gtexGenes.NCBI(b(a));
book(1:10)

fileName = ['~/resultsAndFigures/secondProject/cluster10_SCOverlappingClusters_NCBI.txt']
fid = fopen(fileName, 'w')
for j = 1:length(myIDs)
    thisID = myIDs{j}
    if (length(thisID) > 0 && length(thisID) < 6)
        diff = 6 - length(thisID);
        for i = 1:diff
            thisID = ['0' thisID];
        end
        thisID = ['NM_', thisID]
        fprintf(fid, sprintf('%s\n', thisID))
    end
end
fclose(fid)

% are they similarly correlated in the bulk tissues? (same cluster)

% if not, what is happening, why are some of them in a whole other
% cluster? 

% comparing clusters 9 and 10
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
myList = [mySyms9sc' mySyms10sc'];

[a, b] = ismember(myList, overlapRes.geneSyms);
gSubNet = overlapRes.gtexBinNet005(b, b);
h = figure
heatmap(gSubNet, 'GridVisible', 'off')
h = figure
sSubNet = log2(overlapRes.sumNet(b, b));
heatmap(sSubNet, 'GridVisible', 'off')

[a, b] = ismember(myList, filDataSet.geneSyms);
noVarSubNet = sumNetNoVar(b,b);
varSubNet = sumNetVar(b,b);

h = figure
heatmap(log2(noVarSubNet), 'GridVisible', 'off')
title('no Var')

h = figure
heatmap(log2(varSubNet), 'GridVisible', 'off')
title('Var')

[a, b] = ismember('CALM3', filDataSet.geneSyms)
% functional enrichment
% >>>>>>>>>>>>>>>>>>>>>>
[passedFDRIDs, passedFDRTerms, inGenes, enrichMat, fdrs] = ...
    fenrich_function(ggeneList9, .1);
c9.IDs = passedFDRIDs;
c9.terms = passedFDRTerms;
c9.genes = inGenes;
c9.enrichMat = enrichMat;
c9.fdrs = fdrs;

% getting the involved genes
sib9 = sum(c9.enrichMat');

h = figure
heatmap(c9.genes(sib9>0), c9.terms, c9.enrichMat(sib9>0,:)', 'GridVisible', 'off')

[passedFDRIDs, passedFDRTerms, inGenes, enrichMat, fdrs] = ...
    fenrich_function(ggeneList10, .1);
c10.IDs = passedFDRIDs;
c10.terms = passedFDRTerms;
c10.genes = inGenes;
c10.enrichMat = enrichMat;

sib10 = sum(c10.enrichMat');
h = figure
heatmap(c10.genes(sib10>0), c10.terms, c10.enrichMat(sib10>0,:)', 'GridVisible', 'off')

% functional terms do not overlap and both groups have functional
% term related to neurons, but with different natures.

% case of cluster 239
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

[a, b] = ismember(overlapRes.gtexGeneSyms, ...
                  gtexCluster.syms(gtexCluster.cs1547 == 239));
mySyms = overlapRes.gtexGeneSyms(a)

% 239, sc sum
h = figure
heatmap(mySyms, mySyms, ...
        log2(overlapRes.sumNet(a,a)+1), 'GridVisible', 'off')
colormap(1-gray)
title('cluster 239 - total sum')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%scluster239_totalSum', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 239, sc sum. inh
h = figure
heatmap(mySyms, mySyms, ...
        (overlapRes.inhSum(a,a)./38), 'GridVisible', 'off')
colormap(1-gray)
title('inh')
title('cluster 239 - inh')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%scluster239_inh', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 239, sc sum exc
h = figure
heatmap(mySyms, mySyms, ...
        (overlapRes.excSum(a,a))./23, 'GridVisible', 'off')
colormap(1 - gray)
title('cluster239 - exc')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%scluster239_exc', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 239, GTEX brain
[a, b] = ismember(mySyms, dataSet.genes);
h = figure
heatmap(GTExFiveNets.nets(2).net005(b, b)+0, 'GridVisible', 'off')
colormap(1-gray)
title('cluster239 in GTEx brain')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%scluster239_GTEx_brain', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 239 SC no var % get the network from generateBulkFromSC.m, 447
[a, b] = ismember(mySyms, bulkFromSC.geneSyms)
h = figure
heatmap((sumNet(b(a), b(a)))./15, 'GridVisible', 'off')
colormap(1-gray)
title('cluster 239 - 15 no Var synNets')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%scluster239_16NoVarSynNet', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 239, SC var % 
h = figure
heatmap((sumNet(b(a), b(a)))./15, 'GridVisible', 'off')
colormap(1-gray)
title('cluster 239 - 15 varSyn synNets')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%scluster239_16VarSynNet', figFolder);
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% final part of generate bulk from SC, you can get the no var and
% the high var or low var plots, which mimic (to some level) the
% bulk vs sc results. 
% now I want to do the funtional enrichment thing

list1 = mySyms(2:35);
list2 = mySyms(36:end);

[passedFDRIDs, passedFDRTerms, inGenes, enrichMat] = ...
    fenrich_function(list1, .1);
res1.IDs = passedFDRIDs;
res1.terms = passedFDRTerms;
res1.genes = inGenes;
res1.enrichMat = enrichMat;

h = figure
heatmap(res1.genes, res1.terms, res1.enrichMat')

[passedFDRIDs, passedFDRTerms, inGenes, enrichMat] = ...
    fenrich_function(list2, .1);
res2.IDs = passedFDRIDs;
res2.terms = passedFDRTerms;
res2.genes = inGenes;
res2.enrichMat = enrichMat;
h = figure
heatmap(res2.genes, res2.terms, res2.enrichMat')

% enrichment result for clusters
res1
res2
c9
c10

book = c9.enrichMat;
kado = sum(book');
sum(kado >0)
c9.subMat = book((kado>0), :);
c9.subGenes = c9.genes(kado>0);

book = c10.enrichMat;
kado = sum(book');
sum(kado >0)
c10.subMat = book((kado>0), :);
c10.subGenes = c10.genes(kado>0);

% enrichment results for 239
ttn
% genes involved
book = ttn.mat;
kado = sum(book');
sum(kado > 0)
ttn.subMat = book((kado >0), :);
ttn.subGenes = ttn.genes(kado > 0);

[a, b] = ismember(c9.IDs, ttn.IDs);
sum(a)

[a, b] = ismember(c10.IDs, ttn.IDs);
sum(a)

% get the GABA gene symbols 

[a, b] = sort(overlapRes.geneSyms);
gabGenes = a(158:162)

b(158:162)

h = figure 
heatmap(overlapRes.excSum(b(158:162), b(158:162)))
%uuuh, no they don't show high coexp in inh, they actually show
%high coexp in exc

% case of cluster 64 and the part that it is present in all thigns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[passedFDRIDs, passedFDRTerms, inGenes, enrichMat] = ...
    fenrich_function(exportGeneList, .1)
c64.ids = passedFDRIDs;
c64.terms = passedFDRTerms;
c64.inGenes = exportGeneList;
c64.mat = enrichMat;

h = figure
heatmap(c64.mat)

[a, b] = ismember(exportGeneList, dataSet.genes);

mean(result.r2(b))
h = figure
plot(sort(result.r2(b)))
hold on

% 3. implication of the functional enrichment between clusters 9
% and 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% based on the regression model, each cluster in the gtex bulk
% tissue has a vector of CTP variation. Most of the pyrimidal
% clusters have the same vector. 

% 4. the sum synth network between clusters 9 and 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bulkFromSC.geneSyms;
sumNet;
overlapRes
filSyms = overlapRes.geneSyms(filter);
filSymsOrdered = bulkFromSC.geneSyms(b(a));

[a, b] = ismember(overlapRes.geneSyms(filter), ...
                  bulkFromSC.geneSyms);

sumNetVar0 = sumNet(b(a), b(a));
sumNetVar33 = sumNet(b(a), b(a));

h = figure
heatmap(log2(sumNetVar0(sb, sb)), 'GridVisible', 'off')
colormap(1-gray)
title('var0 - GTExordered ')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%ssumNet_synthvar0_grouped', figFolder);
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


h = figure
heatmap(log2(sumNetVar33(sb, sb)), 'GridVisible', 'off')
colormap(1 - gray)
title('var33 - GTExordered ')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%ssumNet_synthvar33_grouped', figFolder);
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


% I want to plot here cluster 9, ordered by GTEx clusters
load(['~/resultsAndFigures/secondProject/' ...
      'overlappingLinksBetweenSCNets/highRepLinks.mat']) %

load('~/resultsAndFigures/secondProject/gtexClusters.mat')
load('~/data/GTEx/Brain_Cortex_expGenes.mat')

overlapRes
[a, b] = ismember(overlapRes.gtexGeneSyms, gtexCluster.syms);
gclusterInds = a + 0;
gclusterInds(a) = gtexCluster.cs1547(b(a));

template = gclusterInds;
unique(template)
order = [57, 252, 129 242 243 239, 64 168 253 163 169, 88, 203, 78]

for i = 1:length(order)
    template(template == order(i)) = i;
end
template(template > 14) = 0;
[sa, sb] = sort(template);

% 4. they are sorted based on GTEx order now... I will pull out the
% cluster 9 and cluster 10. 

% >>>>> cluster 9
overlapRes
[a, b] = ismember(overlapRes.gtexGeneSyms(overlapRes.clusters ==9), gtexCluster.syms);
gclusterInds = a + 0;
gclusterInds(a) = gtexCluster.cs1547(b(a));

template = gclusterInds;
unique(template)
order = [57, 252, 129 239 242 243, 64 168 253 163 169, 88, 203, 78]

for i = 1:length(order)
    template(template == order(i)) = i;
end
template(template > 14) = 0;
[sa, sb] = sort(template);

% now, plot the heatmap of cluster9, for sa >=0
subMat = overlapRes.sumNet(overlapRes.clusters == 9, overlapRes.clusters ...
                           == 9);
h = figure
heatmap(log2(subMat(sb(sa>0), sb(sa > 0))), 'GridVisible', 'off')
colormap(1 - gray)
title('GTEx Ordered cluster 9')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%scluster9_GTExOrdered_clusters', figFolder);
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% get the order of the genes: 
subList = overlapRes.geneSyms(overlapRes.clusters == 9);
subList9GTExOrdered = subList(sb(sa>0));
%% TAKE IT TO THE expressionInSC and give me the correlaion plot

% and plot the bar plot for it, so we know which cluster is where

h = figure
mybar = hist(sa(sa>0), unique(sa(sa>0)))
bar([mybar; rand(1,length(mybar))], 'Stacked')
legend(order)
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%scluster9_GTExOrdered_clusters_theBarPlot', figFolder);
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% >>>>> cluster 10
overlapRes
[a, b] = ismember(overlapRes.gtexGeneSyms(overlapRes.clusters ==10), gtexCluster.syms);
gclusterInds = a + 0;
gclusterInds(a) = gtexCluster.cs1547(b(a));

template = gclusterInds;
unique(template)
order = [57, 252, 129 239 242 243, 64 168 253 163 169, 88, 203, 78]

for i = 1:length(order)
    template(template == order(i)) = i;
end
template(template > 14) = 0;
[sa, sb] = sort(template);

% now, plot the heatmap of cluster10, for sa >=0
subMat = overlapRes.sumNet(overlapRes.clusters == 10, overlapRes.clusters ...
                           == 10);
h = figure
heatmap(log2(subMat(sb(sa>0), sb(sa > 0))), 'GridVisible', 'off')
colormap(1 - gray)
title('GTEx Ordered cluster 10')
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%scluster10_GTExOrdered_clusters', figFolder);
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% get the order of the genes: 
subList = overlapRes.geneSyms(overlapRes.clusters == 10);
subList10GTExOrdered = subList(sb(sa>0));
%% TAKE IT TO THE expressionInSC and give me the correlaion plot

% and plot the bar plot for it, so we know which cluster is where

h = figure
mybar = hist(sa(sa>0), unique(sa(sa>0)))
bar([mybar; rand(1,length(mybar))], 'Stacked')
legend(order)
figFolder = ['~/resultsAndFigures/secondProject/scSumNetClusters/']
file = sprintf('%scluster10_GTExOrdered_clusters_theBarPlot', figFolder);
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')
