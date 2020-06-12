% for all the genes, in each of the cell types, identify if they
% are markers or not, by expression. 
% 4. bar overlap diagaram for overlap of markers between GTEx and SC

geneList = thisMarkerGenes;

scExpGeneList = filDataSet.geneSyms;

expMat = normExp;
geneListExpMat = zeros(length(scExpGeneList), 69); 
geneListRatioMat = zeros(length(scExpGeneList), 69); 
avgExpMat = zeros(length(scExpGeneList), 69); 
for i = 1:75 % for each cluster - here I get the count of sc types
             % as 69
    i
    net(i).clusterName = clusterMeta.sortedClusterNames(i);
    [a, b] = ismember(clusterMeta.clusters, ...
                      clusterMeta.sortedClusterNames(i));
    %    thisExpMat = log2(expMat(:, a)+1);
    thisExpMat = log2(expMat(:, a)+1);    
    thisSCount = sum(a);
    thisExpGenes = logical(sum(thisExpMat>0, 2) > thisSCount/5);
    thisGeneSyms = filDataSet.geneSyms(thisExpGenes); 
        
    % get the corr and the networks
    smallMat = thisExpMat(thisExpGenes, :);
    [a, b] = ismember(scExpGeneList, thisGeneSyms);
    geneListExpMat(~a, i) = nan;
    book = smallMat(b(a), :);
    
    avgExpMat(a, i) = mean(book');
    
    % get the ratio of expressed
    book(book == 0) = nan;    
    geneListRatioMat(a, i) = sum(~isnan(book'))./ size(book, 2);
    
    % get the expression
     geneListExpMat(a, i) = nanmean(book');
end
geneListExpMat(isnan(geneListExpMat)) = 0;

whos geneListRatioMat
whos geneListExpMat
[a, ctSort] = sort(clusterMeta.sortedClusterNames(1:75));

geneListRatioMatSorted = geneListRatioMat(:, ctSort);
geneListExpMatSorted = geneListExpMat(:, ctSort);
sortedCTlabels = clusterMeta.sortedClusterNames(ctSort);
avgExpMatSorted = avgExpMat(:, ctSort);

whos geneListRatioMatSorted
whos geneListExpMatSorted
whos sortedCTlabels
whos avgExpMatSorted

% if a gene is present in less than 20%, I skip it
% thisRatio = geneListRatioMatSorted;
% thisRatio(geneListRatioMatSorted < .2) = nan;

% thisExp = geneListExpMatSorted;
% thisExp(geneListExpMatSorted < 3) = nan;

% sib = ~isnan(thisExp) + ~isnan(thisRatio);

h = figure
heatmap(geneListExpMatSorted)
colormap(jet)

h = figure
heatmap(avgExpMatSorted)
colormap(jet)

% find the marker genes which represent double expression : diff is
% greater than 1

scMarkerMat = zeros(16789, 7);

% >>>>>>>>>> for astro
kado1 = avgExpMatSorted(:, 1) - avgExpMatSorted(:, 3:end);
kado2 = avgExpMatSorted(:, 2) - avgExpMatSorted(:, 3:end);

sumk1 = sum(kado1' >= 1);
sumk2 = sum(kado2' >= 1);

halva = max([sumk1; sumk2]);

expFactor = max(avgExpMatSorted(:, 1), avgExpMatSorted(:, 2));

scMarkerMat(:, 1) = ((halva == 73) + (expFactor >= 3)') == 2;
scCellType{1} = 'Astrocyte'

% <<<<<<<<<<<<< end of astro

% >>>>>> for pyramidal

kado = mean(avgExpMatSorted(:, 4:27)')' - avgExpMatSorted(:, [1:3, 28:end]);
sumkado = sum(kado' >= 1);

halva = sumkado;

expFactor = mean(avgExpMatSorted(:, 4:27)');

scMarkerMat(:, 6) = ((expFactor>=3)+(halva > 48)) == 2; % I allowed for 5% error : 48+/51
scCellType{6} = 'Pyramidal'
% <<<<<< end of pyramidal

% Endothelial
kado = avgExpMatSorted(:, 3) - avgExpMatSorted(:, [1 2, 4:75]);
sumkado = sum(kado' >= 1);
halva = sumkado;

expFactor = avgExpMatSorted(:, 3) >= 3;

scMarkerMat(:, 2) = ((halva == 74) + expFactor') == 2; % I allowed for 5% error : 48+/51
scCellType{2} = 'Endothelial'

% Oligo
kado = avgExpMatSorted(:, 75) - avgExpMatSorted(:, [1:74]);
sumkado = sum(kado' >= 1);
halva = sumkado;

expFactor = avgExpMatSorted(:, 75) >= 3;

scMarkerMat(:, 5) = ((halva == 74) + expFactor') == 2; % I allowed for 5% error : 48+/51
scCellType{5} = 'Oligo'

% OPC
kado = avgExpMatSorted(:, 74) - avgExpMatSorted(:, [1:73, 75]);
sumkado = sum(kado' >= 1);
halva = sumkado;

expFactor = avgExpMatSorted(:, 74) >= 3;

scMarkerMat(:, 4) = ((halva == 74) + expFactor') == 2; % I allowed for 5% error : 48+/51
scCellType{4} = 'OPC'

% Micro
kado = avgExpMatSorted(:, 73) - avgExpMatSorted(:, [1:72, 74:75]);
sumkado = sum(kado' >= 1);
halva = sumkado;

expFactor = avgExpMatSorted(:, 73) >= 3;

scMarkerMat(:, 3) = ((halva == 74) + expFactor') == 2; % I allowed for 5% error : 48+/51
scCellType{3} = 'Micro'

% >>>>>> for Inh

kado = mean(avgExpMatSorted(:, 28:72)')' - avgExpMatSorted(:, [1:27 73:end]);
sumkado = sum(kado' >= 1);

halva = sumkado;

expFactor = mean(avgExpMatSorted(:, 28:72)');

scMarkerMat(:, 7) = ((expFactor>=3)+(halva > 28)) == 2; % I allowed for 5% error : 48+/51
scCellType{7} = 'Inhibitory'

% <<<<<< end of Inh

%%
scExpGeneList = filDataSet.geneSyms;

scBasedMarkers.markerMat = scMarkerMat;
scBasedMarkers.cellTypes = scCellType;
scBasedMarkers.genes = scExpGeneList;

%whos avgExpMatSorted
scMarkerMat = scBasedMarkers.markerMat;
sortMat = zeros(size(scMarkerMat));
for i = 1:7
    sortMat(:, i) = scMarkerMat(:, i) .* i;
end

book = sum(sortMat');

% removing OPC and Inhibitory
book(book == 4) = 0;
book(book == 7) = 0;

[a, b] = sort(book);
geneList = scBasedMarkers.genes(b(min(find(a>0)):end));
%selectedSortedMat = avgExpMatSorted(b(min(find(a>0)):end), :);
% h = figure
% heatmap(selectedSortedMat');
% colormap(jet)

scBasedMarkers.sortedGenes = scExpGeneList(b(min(find(a>0)):end));
scBasedMarkers.sortedExp = selectedSortedMat;

save('~/resultsAndFigures/secondProject/scBasedMarkers_withInh.mat', ...
     'scBasedMarkers')
load(['~/resultsAndFigures/secondProject/' ...
      'scBasedMarkers_withInh.mat'])

% saveing the bar plot
counts = sum(scBasedMarkers.markerMat);
h = figure
bar(counts, .75)
set(gca, 'XTickLabel', scBasedMarkers.cellTypes)
xtickangle(60)
ylim([0 190])
figFolder = ['~/resultsAndFigures/secondProject/generalFigures/']
file = sprintf('%sSCbasedCountOfMarkers_barplot', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

geneList = scBasedMarkers.sortedGenes;

% get the indexes 
counts = [ min(find(a > 0)) - min(find(a > 0)) + 1
           min(find(a > 1)) - min(find(a > 0))
           min(find(a > 2)) - min(find(a > 1))
           min(find(a > 3)) - min(find(a > 2))
%           min(find(a > 4)) - min(find(a > 3))
           min(find(a > 5)) - min(find(a > 4))
%           min(find(a > 6)) - min(find(a > 5))
%           length(a) - min(find(a > 6))]
                      length(a) - min(find(a > 5))]
cumCounts = cumsum(counts) % the indexes of the markers in the
                           % complete list

% the 'ins' comes from the clusterInNetworks.m
% given an ins, I want to get the count of each cell type

presentCounts = zeros(1, 5);
for i = 1:5 % the ins comes from the clustersInNetworks.m : the
            % genes in affy or GTEx
    presentCounts(i) = sum(ins((cumCounts(i)+1):cumCounts(i+1)));
end
cumPC = cumsum(presentCounts);
cumPC = [0 cumPC];

% get the density of fullNet in each part using the mCounts
fnDensity = zeros(5,5);
for i = 1:5 % on rows
    r1 = cumPC(i) + 1;
    r2 = cumPC(i+1);
    for j = 1:5 % on columns
        totalD = presentCounts(i) * presentCounts(j);
        c1 = cumPC(j) + 1;
        c2 = cumPC(j+1);
        littleFullNet = fullNet(r1:r2, c1:c2);
        fnDensity(i, j) = sum(littleFullNet(:))/totalD;
    end
end

widths = ceil((presentCounts./(sum(presentCounts)))*100)
cumWidths = [0 cumsum(widths)];

normHeatMap = zeros(5, 5);
for i = 1:5
    r1 = cumWidths(i) + 1;
    r2 = cumWidths(i+1);
    for j = 1:5
        c1 = cumWidths(j) + 1;
        c2 = cumWidths(j+1);
        normHeatMap(r1:r2, c1:c2) = fnDensity(i,j);
    end
end
plotMap = normHeatMap ./ .005;

h = figure
heatmap(plotMap, 'GridVisible', 'off')
heatmap(normHeatMap, 'GridVisible', 'off')
title('the Marker Coexpression - GTEx bulk - density')
title('the Marker Coexpression - sc bulkSim - density')

figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%smarkerCoexpression_scBulkSim_density', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 4. bar overlap diagaram for overlap of markers between GTEx and SC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% load the marker genes from SC
load(['~/resultsAndFigures/secondProject/' ...
      'scBasedMarkers_withInh.mat'])

load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])

% for each cell type get the individual counts and overlap, then
% ind - overlap
myBar = zeros(5, 3);

% astrocyte
gtexM = markers.genes( markers.finalMarkerGenes{1});
scM = scBasedMarkers.genes(logical(scBasedMarkers.markerMat(:, ...
                                                  1)));
[a, b] = ismember(scM, gtexM);
myBar(1, 1) = sum(a);
myBar(1, 2) = length(a) - sum(a);
myBar(1, 3) = length(gtexM) - sum(a);

% oligodendrocyte
gtexM = markers.genes( markers.finalMarkerGenes{8});
scM = scBasedMarkers.genes(logical(scBasedMarkers.markerMat(:, ...
                                                  5)));
[a, b] = ismember(scM, gtexM);
myBar(2, 1) = sum(a);
myBar(2, 2) = length(a) - sum(a);
myBar(2, 3) = length(gtexM) - sum(a);

% microglia
gtexM = markers.genes( markers.finalMarkerGenes{5});
scM = scBasedMarkers.genes(logical(scBasedMarkers.markerMat(:, ...
                                                  3)));
[a, b] = ismember(scM, gtexM);
myBar(3, 1) = sum(a);
myBar(3, 2) = length(a) - sum(a);
myBar(3, 3) = length(gtexM) - sum(a);

% pyramical
gtexM = markers.genes( markers.finalMarkerGenes{10});
scM = scBasedMarkers.genes(logical(scBasedMarkers.markerMat(:, ...
                                                  6)));
[a, b] = ismember(scM, gtexM);
myBar(4, 1) = sum(a);
myBar(4, 2) = length(a) - sum(a);
myBar(4, 3) = length(gtexM) - sum(a);

% endothelial
gtexM = markers.genes( markers.finalMarkerGenes{2});
scM = scBasedMarkers.genes(logical(scBasedMarkers.markerMat(:, ...
                                                  2)));
[a, b] = ismember(scM, gtexM);
myBar(5, 1) = sum(a);
myBar(5, 2) = length(a) - sum(a);
myBar(5, 3) = length(gtexM) - sum(a);
myBar1 = myBar;
% each cell type gets a count

load('~/data/GTEx/GTExGeneIDfier.mat')
load('~/data/brainSingleCell/scGeneIDfier.mat')

scEns = sc.ensemble;
for i = 1:length(scEns)
    if (length(scEns{i}) == 0)
        i
        scEns{i} = sprintf('sc%d', i);
    end
end
[a, b] = ismember(markers.genes, gtexGenes.symbols);
gtexEns = gtexGenes.ensemble(b);
for i = 1:length(gtexEns)
    if (length(gtexEns{i}) == 0)
        i
        gtexEns{i} = sprintf('gtex%d', i);
    end
end

% astrocyte
gtexM = gtexEns( markers.finalMarkerGenes{1});
scM = scEns(logical(scBasedMarkers.markerMat(:, ...
                                                  1)));
[a, b] = ismember(scM, gtexM);
myBar(1, 1) = sum(a);
myBar(1, 2) = length(a) - sum(a);
myBar(1, 3) = length(gtexM) - sum(a);

% oligodendrocyte
gtexM = gtexEns( markers.finalMarkerGenes{8});
scM = scEns(logical(scBasedMarkers.markerMat(:, ...
                                                  5)));
[a, b] = ismember(scM, gtexM);
myBar(2, 1) = sum(a);
myBar(2, 2) = length(a) - sum(a);
myBar(2, 3) = length(gtexM) - sum(a);

% microglia
gtexM = gtexEns( markers.finalMarkerGenes{5});
scM = scEns(logical(scBasedMarkers.markerMat(:, ...
                                                  3)));
[a, b] = ismember(scM, gtexM);
myBar(3, 1) = sum(a);
myBar(3, 2) = length(a) - sum(a);
myBar(3, 3) = length(gtexM) - sum(a);

% pyramical
gtexM = gtexEns( markers.finalMarkerGenes{10});
scM = scEns(logical(scBasedMarkers.markerMat(:, ...
                                                  6)));
[a, b] = ismember(scM, gtexM);
myBar(4, 1) = sum(a);
myBar(4, 2) = length(a) - sum(a);
myBar(4, 3) = length(gtexM) - sum(a);

% endothelial
gtexM = gtexEns( markers.finalMarkerGenes{2});
scM = scEns(logical(scBasedMarkers.markerMat(:, ...
                                                  2)));
[a, b] = ismember(scM, gtexM);
myBar(5, 1) = sum(a);
myBar(5, 2) = length(a) - sum(a);
myBar(5, 3) = length(gtexM) - sum(a);

h = figure
bar(myBar, 'stacked')
legend('overlap', 'sc only', 'GTEx only')
set(gca, 'XTickLabel', markers.cellTypes([1 8 5 10 2]))
xtickangle(60)

figFolder = ['~/resultsAndFigures/secondProject/figures/']
file = sprintf('%smarkerOverlaps', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


