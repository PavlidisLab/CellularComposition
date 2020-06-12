% In this file, I correct the cell type proportions for the
% GTEx-brain-cortex and save the dataset as the corrected one.
% 1. data preparation - getting the markers
% 2. do the regression
% 3. Some stuff that I didn't use: protein complex representation
% 4. test to see how does a random selection of count of marker
% genes explain the variation - it should be non-marker genes. 
% 5. correlation of the R2 values for different sets of markers in
% clusters
% 6. getting each clusters R2s variations

clear

% 0. obsolete - markers
%%%%%%%%%%%%%%%%%%%%%%%
% >>>
% markers in GTEx, final markers
load(['~/data/cellTypeVarianceFiles/' ...
      'markersFull.mat'])

% raw markers in Affy  - without filters
load(['~/data/cellTypeVarianceFiles/' ...
      'markerSet.mat'])

% filtered markers in affy
load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes.mat'])
 
% 1. data preparation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.1 loading the GTEx dataset and gettign the brain-cortex file
% 1.2 get the markers for cortex mouse

% 1.1 loading the GTEx dataset and gettign the brain-cortex file
load(['~/data/GTEx/' ...
      'GTExAllDataSetFromGeneLevel_v6p_newWithBlood_RPM.mat'])

uniqueGenes = unique(rpmGTExDS.genes);
myConMap = containers.Map(rpmGTExDS.genes, 1: ...
                          length(rpmGTExDS.genes));

% getting the brainCortex
myDS = rpmGTExDS.dataSets(13);
myTissue = rpmGTExDS.dataSets(13).tissue;

% collapsing the transcripts witnnh the same symbol
myMat = zeros(length(uniqueGenes), size(myDS.mat, 2));

for j = 1:length(myConMap)
    transInds = myConMap(uniqueGenes{j});
    tinyExpMat = myDS.mat(transInds, :);
    if length(transInds) > 1 
        myMat(j, :) = sum(tinyExpMat);
    else 
        myMat(j, :) = tinyExpMat;
    end
end

% the .05 is the count "3" reads divided by the per million
% factor of the average sample in GTEx (57). Basically, I am
% filtering the genes so that the selected genes have >=3
% reads in more than 20% of the samples 
expGenes = (sum(myMat' > .05) / size(myMat, 2)) > .2;
smallMat = myMat(expGenes, :);

dataSet.mat = smallMat;
dataSet.genes = uniqueGenes(expGenes);
save('~/data/GTEx/Brain_Cortex_expGenes.mat', 'dataSet')

load('~/data/GTEx/Brain_Cortex_expGenes.mat')
dataSet.mat = log2(dataSet.mat+1);

% I have my expression mat, now find the markers and correct the
% expression dataset based on them. 

% >>> obsolete
% markers in GTEx, final markers
% load(['~/data/cellTypeVarianceFiles/' ...
%       'markersFull.mat'])

% % raw markers in Affy  - without filters
% load(['~/data/cellTypeVarianceFiles/' ...
%       'markerSet.mat'])

% % filtered markers in affy
% load(['~/data/cellTypeVarianceFiles/' ...
%       'finalMarkerGenes.mat'])
% <<< obsolete

% 1.2 get the markers for cortex mouse

file = '~/data/cellTypeVarianceFiles/Cortex.txt'
fid = fopen(file)

headers = fgets(fid);
cellTypes = strsplit(headers, '\t');

markerListsRaw = textscan(fid, [repmat('%s', 1, 18)], ...
                'Headerlines', 1, 'Delimiter', '\t');

% >> get the homologue gene lists
file = '~/data/cellTypeVarianceFiles/humanAndMouse.txt'
fid = fopen(file)

linum = 0;
homoData = textscan(fid, [repmat('%d', 1, 3), repmat('%s', 1, 3)], ...
                'Headerlines', linum, 'Delimiter', '\t');

% >> get the markers for human 

% for a given set of markers, find the human set of genes 
% the human gene list 
uniqueExpGenes = uniqueGenes(expGenes);
markerSet = zeros(length(uniqueExpGenes), 18);
book = markerListsRaw{1};

% find the human homolog for each of the genes 
for cellType = 1:18
    [a, b] = ismember(markerListsRaw{cellType}, homoData{4});
    IDs = b(a);
    for i = 1:length(IDs)
        thisHGID = homoData{1}(IDs(i));
        homos = find(homoData{1} == thisHGID);
        [ta, tb] = ismember(homoData{2}(homos), 9606);
        if(sum(ta > 0))
            humanIDs = homos(ta);
            humanSyms = homoData{4}(humanIDs);
            % which of them is in the GTEx-brain-cortex genes
            [tta, ttb] = ismember(humanSyms, uniqueExpGenes);
            if(sum(tta))
                markerSet(ttb(tta), cellType) = 1;
            end
        end
    end
end

% raw markersets present in GTEx
save('~/data/cellTypeVarianceFiles/GTEx_brain_cortex_markerSet.mat', ...
     'markerSet')

load(['~/data/cellTypeVarianceFiles/' ...
      'GTEx_brain_cortex_markerSet.mat'])

addpath('~/codes/MATLAB/myCodes/cellTypeVariance/');
presentTypes = [1 2 3 5 10 11 12 13 14 18]; % cell types with
                                            % acceptable count of markers
finalMarkerGenes = cell(1, 10);
for i = 1:length(presentTypes)
    i
    finalMarkerGenes{i} = corrCluster_function_GTExExpGenes(markerSet(:,presentTypes(i)), .6);
end

markers.genes = uniqueExpGenes;
markers.finalMarkerGenes = finalMarkerGenes;
markers.cellTypes = cellTypes(presentTypes);
save('~/data/cellTypeVarianceFiles/finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat', ...
     'markers')

load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])

load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V0_cthr5.mat'])

% plotting the bar plot of the marker gene counts
kado = sum(markerSet);
pmCount = kado(presentTypes);
finalCounts = zeros(1, 10);
for i = 1:10
    finalCounts(i) = length(markers.finalMarkerGenes{i});
end

barMat = [finalCounts; pmCount-finalCounts]
h = figure
bar(barMat', 'stacked', 'EdgeColor', [1 1 1])
set(gca, 'XTickLabel', markers.cellTypes)
xtickangle(60)
legend('selected', 'filtered out')
title(['Count of mouse markers identified in GTEx-cortex dataset'])

figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%scountOfMarkers_barplot_V2', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% >>>>> testing: correcting the markers in the dataset
load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V2.mat'])

myFirstPCs = zeros(114, 10);
mySecondPCs = zeros(114, 10);
for i = 1:10
    markerExp = dataSet.mat(markers.finalMarkerGenes{i}, :);
    [w, ~, ~, ~, exp] = pca(markerExp');
    T = markerExp' * w;
    myFirstPCs(:, i) = T(:, 1);
    mySecondPCs(:, i) = T(:, 2);
    myThirdPCs(:, i) = T(:, 3);
    exps{i} = exp; % variation explained
end
kado1 = corr(myFirstPCs);
kado1 = kado1 - eye(size(kado1));
kado(kado < 0) = 0;

kado2 = corr(mySecondPCs);
kado2 = kado2 - eye(size(kado2));

h = figure
heatmap(markers.cellTypes, [markers.cellTypes], [kado1])
colormap(jet)

figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sMarkerSetCorr_V1', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

h = figure
heatmap(markers.cellTypes, [markers.cellTypes], [kado2])

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])

h = figure
orders = [10 3 4 1 2 5 6 7  8 9]
plotCellTypes = markers.cellTypes(orders)
cellTypesPlusCounts = cell(1, 10);
for i = 1:10
    d = length(markers.finalMarkerGenes{i});
    cellTypesPlusCounts{i} = sprintf('%d - %s', d, markers.cellTypes{i});
end

heatmap(plotCellTypes, cellTypesPlusCounts, kado(orders, orders))
title(['Correlation of the first PC of the marker genes - GTEx brain ' ...
       'cortex'])
set(gca, 'FontSize', 8)
figFolder = ['~/resultsAndFigures/secondProject/networkComplexes/']
file = sprintf('%sComplexesEnrichedInNetworks_1to40', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

modpcs = zeros(size(myFirstPCs));
sib = mean(myFirstPCs);
for i = 1:10
    modpcs(:, i) = myFirstPCs(:, i)./sib(i);
end

% 2. do the regression 

% >>>> Mancarci markers 
load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])


allMarkers = zeros(size(markers.genes));
for i = [1, 2, 5 8 10]
    allMarkers(markers.finalMarkerGenes{i}) = 1;
end

% >>>> snuc markers
load('~/resultsAndFigures/secondProject/scBasedMarkers.mat')
geneList = scBasedMarkers.sortedGenes;

sib = scBasedMarkers.markerMat(:, [1 2 3 5 6]);
geneList = scBasedMarkers.genes(logical(sum(sib')));

[a, b] = ismember(dataSet.genes, geneList);
sum(a)
markerExp = dataSet.mat(a, :);
allMarkers = a;
% >>>>

markerExp = dataSet.mat(logical(allMarkers), :);

[w, score, latent, tsquared, explained, mu] = pca(markerExp', ...
                                                  'NumComponents', 10);
T = markerExp' * w;

% correct it for all the genes
% featureMat = score(:, 1:5);
featureMat = score(:, 1:7);  % count of pcs matter
regOut = zeros(size(dataSet.mat));
%coeffs = zeros(size(dataSet.mat, 1), 6);
coeffs = zeros(size(dataSet.mat, 1), 8);

r2 = zeros(size(dataSet.mat , 1), 1);
for i = 1:size(dataSet.mat, 1)
    i
    y = (dataSet.mat(i, :));
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

% correlation of the weights for the markers. 
for i = 1:size(regOut, 1)
    regOut(i, :)  = regOut(i, :) + coeffs(i, 1);
end

correctedExp = dataSet.mat' - regOut';

result.coeff = coeffs;
result.featureMat = featureMat;
result.regOut = regOut;
result.r2 = r2;

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_5CellTypes_MarkersV1_cthr6.mat'],'result', '-v7.3')

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_AllenSCBasedMarker.mat'],'result', ...
     '-v7.3')

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_AllenSCBasedMarker.mat'],'result', ...
     '-v7.3')

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_AllenSCBasedMarker.mat'],'result', ...
     '-v7.3')

brainRes = load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_AllenSCBasedMarker.mat'])

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_MarkersV2_cthr7.mat'],'result', '-v7.3')

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_MarkersV2_cthr7.mat'])

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_MarkersV2_cthr7.mat'])

load('~/data/GTEx/Brain_Cortex_expGenes.mat')
% make the new binary network and save it 

% >>>>> obsolete
% correctedExp = dataSet.mat' - result.regOut';
% signs = -1*(correctedExp < 0) + (correctedExp>=0);
% mylog = log2(abs(correctedExp)+1);
% final = signs .* mylog;
% <<<<< obsolete

% just the residuals: 
final = correctedExp;

sib = corr(final); 
upCorr = sib(logical(triu(ones(size(sib)), 1)));
qs = quantile(upCorr, [999]);

binNet01 = triu(sib  > qs(990), 1);
ctc.net01 = sparse(binNet01);

binNet005 = triu(sib  > qs(995), 1);
ctc.net005 = sparse(binNet005);
ctc.geneSyms = dataSet.genes;

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'correctedBinNets_logCorrected_jusResiduals_scBasedMarkers' ...
      '.mat'],'ctc')

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'correctedBinNets_logCorrected_jusResiduals_scBasedMarkers' ...
      '_theQuantiles.mat'],'qs')


save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
 ...
 ...
      'correctedBinNets_logCorrected_jusResiduals_scBasedMarkers_7PCA_new.mat'],'ctc')

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
 ...
      'correctedBinNets_logCorrected_jusResiduals_scBasedMarkers_7PCA_new.mat'])


clear
load('~/resultsAndFigures/secondProject/GTExRegression/correctedBinNets.mat')
load('~/resultsAndFigures/secondProject/complexes_attributes.mat')

% 3. Some stuff that I didn't use: protein complex representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the complex representations...

cin = length(com);
myNet = ctc.net01;
thisNetSyms = ctc.geneSyms;
netLinkCount = sum(myNet(:));
netPotCount = length(thisNetSyms) * (length(thisNetSyms) -1) / ...
    2;
genesExpCount = zeros(1, cin);
linkCounts = zeros(1, cin);
for i = 1:cin
    [a, b] = ismember(thisNetSyms, com(i).atomGeneSyms);
    if (sum(a) >= 3)
        genesExpCount(i) = sum(a);
        linkCounts(i) = sum(sum(myNet(a, a)));
    end
end

% get the ps
ps = ones(1, 2004) .* -1;

potLinkCount = genesExpCount .* (genesExpCount - 1) ...
    ./2;
passedComplexes = genesExpCount >= 3;
passedComplexesIDs = find(passedComplexes);
passedCompCount = length(passedComplexesIDs);
thisNetLinkCount = netLinkCount;
thisNetPotCount = netPotCount;
for j = 1:passedCompCount
    myInd = passedComplexesIDs(j);
    ps(myInd) = 1 - hygecdf(linkCounts(myInd)-1, thisNetPotCount, ...
                            thisNetLinkCount, potLinkCount(myInd));
end

fdrFac = sum(ps >=0)

fdrPs = ps .* fdrFac;

complexRes_GTExCorrected.genesExpCount = genesExpCount;
complexRes_GTExCorrected.linkCounts = linkCounts;
complexRes_GTExCorrected.netLinkCount = netLinkCount;
complexRes_GTExCorrected.netPotCount = netPotCount;
complexRes_GTExCorrected.ps = ps;
complexRes_GTExCorrected.fdrs = fdrPs;

save('~/resultsAndFigures/secondProject/complexPresentation_GTExBrainCorrectedNet.mat', ...
     'complexRes_GTExCorrected')

res1 = fdrPs <= .1;
res2 = fdrPs >= 0;
tempRes = res1 + res2;
correctedComp = (tempRes == 2);

mypresentComps = book == 2;

halva = mypresentComps(, :);

sum((halva + correctedComp) == 2)

% TODO 1: this one
% add the GTEx corrected to this 
% 4. test to see how does a random selection of count of marker
% genes explain the variation - it should be non-marker genes. 
% TODO 2: get the pc from the first two pcs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('~/data/GTEx/Brain_Cortex_expGenes.mat')
dataSet.mat = log2(dataSet.mat+1);

% pick one marker set
% >>>> Mancarci markers 
load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])

allMarkers = zeros(size(markers.genes));
for i = [1, 2, 5 8 10]
    allMarkers(markers.finalMarkerGenes{i}) = 1;
end

markerExp = dataSet.mat(a, :);

% >>>> snuc markers
load('~/resultsAndFigures/secondProject/scBasedMarkers.mat')
geneList = scBasedMarkers.sortedGenes;

sib = scBasedMarkers.markerMat(:, [1 2 3 5 6]);
geneList = scBasedMarkers.genes(logical(sum(sib')));

[a, b] = ismember(dataSet.genes, geneList);
sum(a)
markerExp = dataSet.mat(a, :);
allMarkers = a;
% >>>>

% get the expression level of marker genes

meanExps = (mean(dataSet.mat'));
markerMeans = (mean(markerExp'));

h = figure
hist(markerMeans, 9)

expCounts = hist(markerMeans, [.2:1.2:11.2])
ecSum = cumsum(expCounts);

dataSetExp = dataSet.mat;

% I need the same count of inds for each group
borders = [-.4:1.2:11.8];

for j = 1:100
    j
    mockMarkerInds = zeros(1, length(markerMeans));

    tempList = find(((meanExps <= borders(2)) + (~allMarkers)) == 2);
    mockMarkerInds(1:ecSum(1)) = datasample(tempList, ecSum(1));
    for i = 2:(length(borders) - 1)
        tempList = find(((meanExps <= borders(i +1)) + (~allMarkers) + ...
                         (meanExps > borders(i))) == 3);
        mockMarkerInds((ecSum(i-1)+1):ecSum(i)) = datasample(tempList, ...
                                                          expCounts(i));
    end

    mockMarkerExp = dataSetExp(mockMarkerInds, :);
    [w, score, latent, tsquared, explained, mu] = pca(mockMarkerExp', ...
                                                      'NumComponents', ...
                                                      10);
    explained(1:10)
    T = markerExp' * w;

    % correct it for all the genes
    featureMat = score(:, 1:7);
    regOut = zeros(size(dataSetExp));
    coeffs = zeros(size(dataSetExp, 1), 8);
    r2 = zeros(size(dataSetExp , 1), 1);
    for i = 1:size(dataSetExp, 1)
        j
        i
        y = dataSetExp(i, :);
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


    % for i = 1:size(regOut, 1)
    %     regOut(i, :)  = regOut(i, :) + coeff(i, 1);
    % end


    %correctedExp = dataSetExp' - regOut';

    result(j).coeff = coeffs;
    result(j).featureMat = featureMat;
    %    result(j).regOut = regOut;
    result(j).r2 = r2;
    result(j).mmg = mockMarkerInds;
end

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_5CellTypes_100Mock_scBasedMarkers_redo' ...
      '.mat'],'result', '-v7.3')

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_5CellTypes_100Mock_MarkersV1_cthr6' ...
      '.mat'],'result', '-v7.3')

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_10CellTypes_100Mock_scBasedMarkers_redo' ...
      '.mat'],'result', '-v7.3')

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_100Mock_scBasedMarkers_redo' ...
      '.mat'],'result', '-v7.3')

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_100Mock_scBasedMarkers_redo' ...
      '.mat'])

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_100Mock_MarkersV0_cthr5' ...
      '.mat'],'result', '-v7.3')

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_100Mock_MarkersV2_cthr7' ...
      '.mat'],'result', '-v7.3')

% >>> do the plot
clear

load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])

markerInds = markers.finalMarkerGenes{1}';
for i = 2:10
    markerInds = [markerInds, markers.finalMarkerGenes{i}'];
end

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_100Mock_MarkersV1_ctrh6' ...
      '.mat'])
mockResGTEx = result;

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_100Mock_scBasedMarkers' ...
      '.mat'])
mockResSC = result;

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_AllenSCBasedMarker_redo.mat'])
mainResSC = result;

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_MarkersV2_cthr7.mat'])

load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])
mainResGTEx = result;

load('~/data/GTEx/Brain_Cortex_expGenes.mat')
[a, b] = ismember(scBasedMarkers.sortedGenes, dataSet.genes);
markerInds = b(a);

markerMainRs = mainResSC.r2(markerInds);

ps = zeros(1, 100);
means = zeros(1, 100);
for i = 1:100
    thisRs = mockResSC(i).r2(markerInds);
    means(i) = mean(thisRs);
    [h, p] = ttest(thisRs, markerMainRs);
    ps(i) = p;
end

% plot the distribution of R2 for markers 
xs = zeros(100, 100);
fs = zeros(100, 100);
for i = 1:100 
    i
    [f, x] = ksdensity(mockRes(i).r2(markerInds));
    xs(i, :) = x;
    fs(i, :) = f;
end

h = figure;
for i = 1:100
    plot(xs(i, :), fs(i,:), 'k');
    hold on
end

[f, x] = ksdensity(mainRes.r2(markerInds));
plot(x, f, 'r')

% the above plot doesn't show it right. I will do sth else. 

% SELECTED PLOT! % quantiles of the markers
%pData = zeros(100, length(markerInds));
pData = zeros(101, 20);
for i = 1:100 
    i
    %    pData(i, :) = sort(mockRes(i).r2(markerInds));
    pData(i, :) = quantile(mockRes(i).r2(markerInds), 20);
end
pData(101, :) = quantile(mainRes.r2(markerInds), 20);

% SELECTED PLOT! % quantiles of the markers
h = figure
for i = 1:100
    plot(pData(i, :), 'k')
    hold on
end
plot(pData(101, :), 'r')

title(['Qantiles for the Rsqrd values for marker genes - random vs ' ...
       'marker genes'])
       
set(gca, 'FontSize', 12)
xlabel('quantiles')
ylabel('Rsqrd')
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%smarkerGeneR2Dist_scBasedMarkers', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

plot(sort(mainRes.r2(markerInds)), 'r')

pData = [pData; mainRes.r2(markerInds)'];

spData = (sort(pData'));
csspData = cumsum(spData);

h = figure
for i = 1:100
    plot(csspData(:, i) , 'k')
    hold on
end
plot(csspData(:, 101),'r')

h = figure
hist(mainRes.r2(markerInds))

book = mean(csspData');
hold on
plot(book, 'b')

% plot the distribution of R2 for the genes in the network

% getting the genes in the network
load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat')
myNet = GTExFiveNets.nets(2).net005;
fullNet = myNet + myNet';
nds = sum(fullNet);
inGenes = nds >0; 
sum(inGenes)

qs = zeros(101, 100);
for i = 1:100
    qs(i, :) = quantile(mockRes(i).r2, 100);
end
qs(101, :) = quantile(mainRes.r2, 100);

h = figure
for i = 1:100
    plot(qs(i, :), 'k');
    hold on
end
plot(qs(101, :), 'r')
title(['Qantiles for the Rsqrd values for networkGenes - random vs ' ...
       'marker genes'])
       
set(gca, 'FontSize', 12)
xlabel('quantiles')
ylabel('Rsqrd')
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%snetworkGeneR2Dist_scBasedMarkers', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% plot the distribution of R2 for all the 26k expressed genes

mr = result;
h = figure
plot(sort(result.r2))
hold on 
i = 1
plot(sort(mr(i).r2))
hold on
i = i + 1

h = figure
for i = 1:9
    subplot(2, 5, i)
    % plot(result.r2(logical(allMarkers)), mr(i).r2(logical(allMarkers)), ...
    %      '.')
    hist(result(i).r2(logical(allMarkers)))
    xlim([0, 1])
    ylim([0 150])
    sum((result(i).r2(logical(allMarkers))) > .5)
end

h = figure
hist(mainRes.r2(logical(allMarkers)))
ylim([0 150])
sum(mainRes.r2(logical(allMarkers)) > .5)

xis = zeros(8, 100);
fis = zeros(8, 100);
for i = 1:8
    [f, x] = ksdensity(mr(i).r2(logical(allMarkers)));
    xis(i, :) = x;
    fis(i, :) = f;
end

h = figure
hold on
for i = 1:8
    plot(xis(i, :), fis(i,:))
end

h = figure
[f, x] = ksdensity(result.r2(logical(allMarkers)));
plot(x, f)

% 3. correlation of the R2 values for different sets of markers in
% clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

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

% load the gtex clusters ...
load('~/resultsAndFigures/secondProject/gtexClusters.mat') %


[a, b] = ismember(gtexCluster.syms, dataSet.genes);

sib = corr(resGTEx.r2(b), resSC.r2(b))

h = figure
scatter(resGTEx.r2(b), resSC.r2(b), 'filled', 'SizeData', 7)
alpha(.3)
xlabel('mouse markers')
ylabel('SC based markers')
set(gca, 'FontSize', 18)
set(gca, 'XTick', [0 .5, 1])
set(gca, 'YTick', [0 .5, 1])
xlim([-.1, 1])
ylim([-.1, 1])
title(sprintf('Pearson correlation %.2f', sib))
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%sscatter_R2_genesInGTExClusters_', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% get the distribution of R2 (both) - and the markers for each other
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('~/data/GTEx/Brain_Cortex_expGenes.mat')

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_AllenSCBasedMarker.mat'])
scRes = result;
clear result;

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_MarkersV1_cthr6.mat'])
gtexRes = result;
clear result

% get the null distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('~/resultsAndFigures/secondProject/gtexClusters.mat') %

% do the density plot for the resSC, resGTEx and also the genes in
% the network
[a, b] = ismember(gtexCluster.syms, dataSet.genes);

[yg, xg] = ksdensity(gtexRes.r2);
[ys, xs] = ksdensity(scRes.r2);
[ygc, xgc] = ksdensity(gtexRes.r2(b(a)));
[ysc, xsc] = ksdensity(scRes.r2(b(a)));

h = figure
plot(xg, yg)
hold on
plot(xs, ys)
plot(xgc, ygc)
plot(xsc, ysc)

legend('all genes, MouseMarkers', 'allGenes, SCmarkers', ['networkGenes, ' ...
                    'mouse markers'], 'networkGenes, scmarkers')
title('Distribution of r2 values')
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%sdistributionOfR2values', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])
gtexSet = [];
for i = 1:10
    gtexSet = [gtexSet markers.genes(markers.finalMarkerGenes{i})];
end

load('~/resultsAndFigures/secondProject/scBasedMarkers.mat')
scSet = scBasedMarkers.genes(sum(scBasedMarkers.markerMat(:, [1 2 3 ...
                   5 6])') > 0);

% gtex markers, sc R2
[a, b] = ismember(gtexSet, dataSet.genes);

bdatag = [gtexRes.r2(b(a))' , scRes.r2(b(a))'];
bgroupg = [ones(sum(a), 1)' ones(sum(a), 1)'.*2];

h = figure
boxplot(bdatag, bgroupg)
hold on
scatter(bgroupg-.25, bdatag)
line([.5, 1.25], [mean(gtexRes.r2), mean(gtexRes.r2)])
line([1.5, 2.25], [mean(scRes.r2), mean(scRes.r2)])
ylabel('r2')
title(sprintf('gtexMarkerBased r2 and SCMarkerBased r2 \n for gtexbased markers'))
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%sboxplotOfR2ForMarkerGenes_MouseMarkers', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% sc markers, gtex R2
[a, b] = ismember(scSet, dataSet.genes);

bdatas = [gtexRes.r2(b(a))' , scRes.r2(b(a))'];
bgroups = [ones(sum(a), 1)' ones(sum(a), 1)'.*2];
h = figure
boxplot(bdatas, bgroups)
hold on
scatter(bgroups-.25, bdatas)
line([.5, 1.25], [mean(gtexRes.r2), mean(gtexRes.r2)])
line([1.5, 2.25], [mean(scRes.r2), mean(scRes.r2)])

ylabel('r2')
title(sprintf('gtexMarkerBased r2 and SCMarkerBased r2 \n for sc markers'))
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%sboxplotOfR2ForMarkerGenes_scBasedMarkers', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% SELECTED PLOT! % quantiles of the markers
%pData = zeros(100, length(markerInds));
pData = zeros(101, 20);
for i = 1:100 
    i
    %    pData(i, :) = sort(mockRes(i).r2(markerInds));
    pData(i, :) = quantile(mockRes(i).r2(markerInds), 20);
end
pData(101, :) = quantile(scRes.r2(b(a)), 20);

% let's see at each quantile how many are bigger 
bigCounts = zeros(1, 20);
for i = 1:20
    bigCounts(i) = sum(pData(1:100, i) >= pData(101, i));
end

    
% plot the quantiles 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mockResSC
scSet
scRes

pData = zeros(101, 20);
[a, b] = ismember(gtexSet, dataSet.genes);
markerInds = b(a);
for i = 1:100 
    i
    %    pData(i, :) = sort(mockRes(i).r2(markerInds));
    pData(i, :) = quantile(mockResSC(i).r2(markerInds), 20);
end
pData(101, :) = quantile(scRes.r2(markerInds), 20);

bigCounts = zeros(1, 20);
for i = 1:20
    bigCounts(i) = sum(pData(1:100, i) >= pData(101, i));
end
bigCounts

% I should get the GTEx markers, from the SC based markers, and
% then comapre it with the mock...
h = figure
plot([1:20], pData(101,:))
hold on
ylim([-.2 1])
title('GTEx marker set predicted by SC marker set')
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%sgtexMarkersPredictedBySCMarkers', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

mockResGTEx
gtexSet
gtexRes

pData = zeros(101, 20);
[a, b] = ismember(scSet, dataSet.genes);
markerInds = b(a);
for i = 1:100 
    i
    %    pData(i, :) = sort(mockRes(i).r2(markerInds));
    pData(i, :) = quantile(mockResGTEx(i).r2(markerInds), 20);
end
pData(101, :) = quantile(gtexRes.r2(markerInds), 20);

bigCounts = zeros(1, 20);
for i = 1:20
    bigCounts(i) = sum(pData(1:100, i) >= pData(101, i));
end
bigCounts

h = figure
plot([1:20], pData(101,:))
hold on
ylim([0 .9])
title('SC marker set predicted by gtex marker set')
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%sscMarkersPredictedByGTExMarkers', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 4. getting each clusters R2s variations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
load('~/resultsAndFigures/secondProject/gtexClusters.mat') %

load('~/data/GTEx/Brain_Cortex_expGenes.mat')

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_AllenSCBasedMarker.mat'])
scRes = result;
clear result;

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_MarkersV1_cthr6.mat'])
gtexRes = result;
clear result

load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'GTExBloodClusterRep.mat']) % GTExClusterRep_Gblood.clusterInds forthe
myInds = GTExClusterRep_Gblood.clusterInds;

% getting the mat for a given cluster ID 
bs % ordered IDs from clusteringComparison

means = zeros(69, 7);
stds = zeros(69, 7);
for i = 1:length(myInds)
    mySyms = gtexCluster.syms(gtexCluster.cs1547 == myInds(i));
    [a, b] = ismember(mySyms, dataSet.genes);

    myValues = scRes.coeff(b(a), 2:end);
    means(i, :) = mean(myValues);
    stds(i, :) = std(myValues);
    % h = figure
    % errorbar(means(i, :), stds(i, :))
    % sib = gtexCluster.clusterLabels(myInds(bs(i)));
    % title(sprintf('%s', sib{1}))
    % ylim([-.03, .055])
    % figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
    % file = sprintf('%sR2plot_%s', figFolder, sib{1})
    % set(h, 'PaperOrientation', 'landscape')
    % print(h, '-deps', [file '.eps'])
    % print(h, '-dpdf', [file '.pdf'])
    % saveas(h, [file '.eps'], 'epsc')
end

sib = corr(means');
h = figure
heatmap(gtexCluster.clusterLabels(myInds(bs(1:32))), gtexCluster.clusterLabels(myInds(bs(1:32))),sib(1:32, 1:32))
colormap(hot)

boxplot(myValues)
mm = mean(myValues);
mv = std(myValues);
h = figure
errorbar(mm, mv)

% bs comes from the order of the figure, I give it in
% clusteringComparison.m

kado = means(bs(3:20), :);
sib = corr(kado');
h = figure
h = figure('units', 'centimeters', 'position', [0,0, 22, 22])
heatmap(gtexCluster.clusterLabels(myInds(bs(3:20))), ...
        gtexCluster.clusterLabels(myInds(bs(3:20))),sib, 'GridVisible', 'off')
colormap(bone)
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%sCTVCorr_AstroPyramidal_selected', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')



T = markerExp' * w;

% correct it for all the genes
featureMat = score(:, 1:7);
regOut = zeros(size(dataSet.mat));
coeffs = zeros(size(dataSet.mat, 1), 8);
r2 = zeros(size(dataSet.mat , 1), 1);
for i = 1:size(dataSet.mat, 1)
    i
    y = (dataSet.mat(i, :));
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

% correlation of the weights for the markers. 
for i = 1:size(regOut, 1)
    regOut(i, :)  = regOut(i, :) + coeffs(i, 1);
end

correctedExp = dataSet.mat' - regOut';

result.coeff = coeffs;
result.featureMat = featureMat;
result.regOut = regOut;
result.r2 = r2;
save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_AllenSCBasedMarker.mat'],'result', ...
     '-v7.3')

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_AllenSCBasedMarker.mat'],'result', ...
     '-v7.3')

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_AllenSCBasedMarker.mat'])

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_MarkersV2_cthr7.mat'],'result', '-v7.3')

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_MarkersV2_cthr7.mat'])

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_MarkersV2_cthr7.mat'])

load('~/data/GTEx/Brain_Cortex_expGenes.mat')
% make the new binary network and save it 
correctedExp = dataSet.mat' - result.regOut';
signs = -1*(correctedExp < 0) + (correctedExp>=0);
mylog = log2(abs(correctedExp)+1);
final = signs .* mylog;

% OR >> just the residuals: 
final = correctedExp;

sib = corr(final); 
upCorr = sib(logical(triu(ones(size(sib)), 1)));
qs = quantile(upCorr, [999]);

binNet01 = triu(sib  > qs(990), 1);
ctc.net01 = sparse(binNet01);

binNet005 = triu(sib  > qs(995), 1);
ctc.net005 = sparse(binNet005);
ctc.geneSyms = dataSet.genes;

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'correctedBinNets_logCorrected_jusResiduals_scBasedMarkers' ...
      '.mat'],'ctc')

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'correctedBinNets_logCorrected_jusResiduals_scBasedMarkers' ...
      '_theQuantiles.mat'],'qs')


save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
 ...
 ...
      'correctedBinNets_logCorrected_jusResiduals_scBasedMarkers_7PCA_new.mat'],'ctc')

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
 ...
      'correctedBinNets_logCorrected_jusResiduals_scBasedMarkers_7PCA_new.mat'])


clear
load('~/resultsAndFigures/secondProject/GTExRegression/correctedBinNets.mat')
load('~/resultsAndFigures/secondProject/complexes_attributes.mat')

% get the complex representations...
cin = length(com);
myNet = ctc.net01;
thisNetSyms = ctc.geneSyms;
netLinkCount = sum(myNet(:));
netPotCount = length(thisNetSyms) * (length(thisNetSyms) -1) / ...
    2;
genesExpCount = zeros(1, cin);
linkCounts = zeros(1, cin);
for i = 1:cin
    [a, b] = ismember(thisNetSyms, com(i).atomGeneSyms);
    if (sum(a) >= 3)
        genesExpCount(i) = sum(a);
        linkCounts(i) = sum(sum(myNet(a, a)));
    end
end


% get the ps
ps = ones(1, 2004) .* -1;

potLinkCount = genesExpCount .* (genesExpCount - 1) ...
    ./2;
passedComplexes = genesExpCount >= 3;
passedComplexesIDs = find(passedComplexes);
passedCompCount = length(passedComplexesIDs);
thisNetLinkCount = netLinkCount;
thisNetPotCount = netPotCount;
for j = 1:passedCompCount
    myInd = passedComplexesIDs(j);
    ps(myInd) = 1 - hygecdf(linkCounts(myInd)-1, thisNetPotCount, ...
                            thisNetLinkCount, potLinkCount(myInd));
end

fdrFac = sum(ps >=0)

fdrPs = ps .* fdrFac;

complexRes_GTExCorrected.genesExpCount = genesExpCount;
complexRes_GTExCorrected.linkCounts = linkCounts;
complexRes_GTExCorrected.netLinkCount = netLinkCount;
complexRes_GTExCorrected.netPotCount = netPotCount;
complexRes_GTExCorrected.ps = ps;
complexRes_GTExCorrected.fdrs = fdrPs;

save('~/resultsAndFigures/secondProject/complexPresentation_GTExBrainCorrectedNet.mat', ...
     'complexRes_GTExCorrected')

res1 = fdrPs <= .1;
res2 = fdrPs >= 0;
tempRes = res1 + res2;
correctedComp = (tempRes == 2);

mypresentComps = book == 2;

halva = mypresentComps(, :);

sum((halva + correctedComp) == 2)

% TODO 1: this one
% add the GTEx corrected to this 
% 2. test to see how does a random selection of count of marker
% genes explain the variation - it should be non-marker genes. 
% TODO 2: get the pc from the first two pcs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I need similar expression levels: Does the R2 change?

% get the expression level of marker genes

meanExps = (mean(log2(dataSet.mat' + 1)));
markerMeans = (mean(markerExp'));

dataSetExp = log2(dataSet.mat + 1);

h = figure
hist(markerMeans, 9)
expCounts = hist(markerMeans, [.2:1.2:11.2])
ecSum = cumsum(expCounts);

% I need the same count of inds for each group
borders = [-.4:1.2:11.8];

for j = 1:100
    j
    mockMarkerInds = zeros(1, length(markerMeans));

    tempList = find(((meanExps <= borders(2)) + (~allMarkers)) == 2);
    mockMarkerInds(1:ecSum(1)) = datasample(tempList, ecSum(1));
    for i = 2:(length(borders) - 1)
        tempList = find(((meanExps <= borders(i +1)) + (~allMarkers) + ...
                         (meanExps > borders(i))) == 3);
        mockMarkerInds((ecSum(i-1)+1):ecSum(i)) = datasample(tempList, ...
                                                          expCounts(i));
    end

    mockMarkerExp = dataSetExp(mockMarkerInds, :);
    [w, score, latent, tsquared, explained, mu] = pca(mockMarkerExp', ...
                                                      'NumComponents', ...
                                                      10);
    explained(1:10)
    T = markerExp' * w;

    % correct it for all the genes
    featureMat = score(:, 1:5);
    regOut = zeros(size(dataSetExp));
    coeffs = zeros(size(dataSetExp, 1), 6);
    r2 = zeros(size(dataSetExp , 1), 1);
    for i = 1:size(dataSetExp, 1)
        j
        i
        y = dataSetExp(i, :);
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


    % for i = 1:size(regOut, 1)
    %     regOut(i, :)  = regOut(i, :) + coeff(i, 1);
    % end


    %correctedExp = dataSetExp' - regOut';

    result(j).coeff = coeffs;
    result(j).featureMat = featureMat;
    %    result(j).regOut = regOut;
    result(j).r2 = r2;
    result(j).mmg = mockMarkerInds;
end


save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_100Mock_scBasedMarkers_redo' ...
      '.mat'],'result', '-v7.3')

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_100Mock_scBasedMarkers_redo' ...
      '.mat'])

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_100Mock_MarkersV0_cthr5' ...
      '.mat'],'result', '-v7.3')

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_100Mock_MarkersV2_cthr7' ...
      '.mat'],'result', '-v7.3')

% >>> do the plot
clear

load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])

markerInds = markers.finalMarkerGenes{1}';
for i = 2:10
    markerInds = [markerInds, markers.finalMarkerGenes{i}'];
end

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_100Mock_MarkersV1_ctrh6' ...
      '.mat'])
mockResGTEx = result;

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_100Mock_scBasedMarkers' ...
      '.mat'])
mockResSC = result;

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_AllenSCBasedMarker_redo.mat'])
mainResSC = result;

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_MarkersV2_cthr7.mat'])

load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])
mainResGTEx = result;

load('~/data/GTEx/Brain_Cortex_expGenes.mat')
[a, b] = ismember(scBasedMarkers.sortedGenes, dataSet.genes);
markerInds = b(a);

markerMainRs = mainResSC.r2(markerInds);

ps = zeros(1, 100);
means = zeros(1, 100);
for i = 1:100
    thisRs = mockResSC(i).r2(markerInds);
    means(i) = mean(thisRs);
    [h, p] = ttest(thisRs, markerMainRs);
    ps(i) = p;
end

% plot the distribution of R2 for markers 
xs = zeros(100, 100);
fs = zeros(100, 100);
for i = 1:100 
    i
    [f, x] = ksdensity(mockRes(i).r2(markerInds));
    xs(i, :) = x;
    fs(i, :) = f;
end

h = figure;
for i = 1:100
    plot(xs(i, :), fs(i,:), 'k');
    hold on
end

[f, x] = ksdensity(mainRes.r2(markerInds));
plot(x, f, 'r')

% the above plot doesn't show it right. I will do sth else. 

% SELECTED PLOT! % quantiles of the markers
%pData = zeros(100, length(markerInds));
pData = zeros(101, 20);
for i = 1:100 
    i
    %    pData(i, :) = sort(mockRes(i).r2(markerInds));
    pData(i, :) = quantile(mockRes(i).r2(markerInds), 20);
end
pData(101, :) = quantile(mainRes.r2(markerInds), 20);

% SELECTED PLOT! % quantiles of the markers
h = figure
for i = 1:100
    plot(pData(i, :), 'k')
    hold on
end
plot(pData(101, :), 'r')

title(['Qantiles for the Rsqrd values for marker genes - random vs ' ...
       'marker genes'])
       
set(gca, 'FontSize', 12)
xlabel('quantiles')
ylabel('Rsqrd')
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%smarkerGeneR2Dist_scBasedMarkers', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

plot(sort(mainRes.r2(markerInds)), 'r')

pData = [pData; mainRes.r2(markerInds)'];

spData = (sort(pData'));
csspData = cumsum(spData);

h = figure
for i = 1:100
    plot(csspData(:, i) , 'k')
    hold on
end
plot(csspData(:, 101),'r')

h = figure
hist(mainRes.r2(markerInds))

book = mean(csspData');
hold on
plot(book, 'b')

% plot the distribution of R2 for the genes in the network

% getting the genes in the network
load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat')
myNet = GTExFiveNets.nets(2).net005;
fullNet = myNet + myNet';
nds = sum(fullNet);
inGenes = nds >0; 
sum(inGenes)

qs = zeros(101, 100);
for i = 1:100
    qs(i, :) = quantile(mockRes(i).r2, 100);
end
qs(101, :) = quantile(mainRes.r2, 100);

h = figure
for i = 1:100
    plot(qs(i, :), 'k');
    hold on
end
plot(qs(101, :), 'r')
title(['Qantiles for the Rsqrd values for networkGenes - random vs ' ...
       'marker genes'])
       
set(gca, 'FontSize', 12)
xlabel('quantiles')
ylabel('Rsqrd')
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%snetworkGeneR2Dist_scBasedMarkers', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% plot the distribution of R2 for all the 26k expressed genes

mr = result;
h = figure
plot(sort(result.r2))
hold on 
i = 1
plot(sort(mr(i).r2))
hold on
i = i + 1

h = figure
for i = 1:9
    subplot(2, 5, i)
    % plot(result.r2(logical(allMarkers)), mr(i).r2(logical(allMarkers)), ...
    %      '.')
    hist(result(i).r2(logical(allMarkers)))
    xlim([0, 1])
    ylim([0 150])
    sum((result(i).r2(logical(allMarkers))) > .5)
end

h = figure
hist(mainRes.r2(logical(allMarkers)))
ylim([0 150])
sum(mainRes.r2(logical(allMarkers)) > .5)

xis = zeros(8, 100);
fis = zeros(8, 100);
for i = 1:8
    [f, x] = ksdensity(mr(i).r2(logical(allMarkers)));
    xis(i, :) = x;
    fis(i, :) = f;
end

h = figure
hold on
for i = 1:8
    plot(xis(i, :), fis(i,:))
end

h = figure
[f, x] = ksdensity(result.r2(logical(allMarkers)));
plot(x, f)

% 3. correlation of the R2 values for different sets of markers in
% clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

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

% load the gtex clusters ...
load('~/resultsAndFigures/secondProject/gtexClusters.mat') %


[a, b] = ismember(gtexCluster.syms, dataSet.genes);

sib = corr(resGTEx.r2(b), resSC.r2(b))

h = figure
scatter(resGTEx.r2(b), resSC.r2(b), 'filled', 'SizeData', 7)
alpha(.3)
xlabel('mouse markers')
ylabel('SC based markers')
set(gca, 'FontSize', 18)
set(gca, 'XTick', [0 .5, 1])
set(gca, 'YTick', [0 .5, 1])
xlim([-.1, 1])
ylim([-.1, 1])
title(sprintf('Pearson correlation %.2f', sib))
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%sscatter_R2_genesInGTExClusters_', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% get the distribution of R2 (both) - and the markers for each other
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('~/data/GTEx/Brain_Cortex_expGenes.mat')

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_AllenSCBasedMarker.mat'])
scRes = result;
clear result;

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_MarkersV1_cthr6.mat'])
gtexRes = result;
clear result

% get the null distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('~/resultsAndFigures/secondProject/gtexClusters.mat') %

% do the density plot for the resSC, resGTEx and also the genes in
% the network
[a, b] = ismember(gtexCluster.syms, dataSet.genes);

[yg, xg] = ksdensity(gtexRes.r2);
[ys, xs] = ksdensity(scRes.r2);
[ygc, xgc] = ksdensity(gtexRes.r2(b(a)));
[ysc, xsc] = ksdensity(scRes.r2(b(a)));

h = figure
plot(xg, yg)
hold on
plot(xs, ys)
plot(xgc, ygc)
plot(xsc, ysc)

legend('all genes, MouseMarkers', 'allGenes, SCmarkers', ['networkGenes, ' ...
                    'mouse markers'], 'networkGenes, scmarkers')
title('Distribution of r2 values')
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%sdistributionOfR2values', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


load(['~/data/cellTypeVarianceFiles/' ...
      'finalMarkerGenes_GTExBrainCortex_V1_cthr6.mat'])
gtexSet = [];
for i = 1:10
    gtexSet = [gtexSet markers.genes(markers.finalMarkerGenes{i})];
end

load('~/resultsAndFigures/secondProject/scBasedMarkers.mat')
scSet = scBasedMarkers.genes(sum(scBasedMarkers.markerMat(:, [1 2 3 ...
                   5 6])') > 0);

% gtex markers, sc R2
[a, b] = ismember(gtexSet, dataSet.genes);

bdatag = [gtexRes.r2(b(a))' , scRes.r2(b(a))'];
bgroupg = [ones(sum(a), 1)' ones(sum(a), 1)'.*2];

h = figure
boxplot(bdatag, bgroupg)
hold on
scatter(bgroupg-.25, bdatag)
line([.5, 1.25], [mean(gtexRes.r2), mean(gtexRes.r2)])
line([1.5, 2.25], [mean(scRes.r2), mean(scRes.r2)])
ylabel('r2')
title(sprintf('gtexMarkerBased r2 and SCMarkerBased r2 \n for gtexbased markers'))
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%sboxplotOfR2ForMarkerGenes_MouseMarkers', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% sc markers, gtex R2
[a, b] = ismember(scSet, dataSet.genes);

bdatas = [gtexRes.r2(b(a))' , scRes.r2(b(a))'];
bgroups = [ones(sum(a), 1)' ones(sum(a), 1)'.*2];
h = figure
boxplot(bdatas, bgroups)
hold on
scatter(bgroups-.25, bdatas)
line([.5, 1.25], [mean(gtexRes.r2), mean(gtexRes.r2)])
line([1.5, 2.25], [mean(scRes.r2), mean(scRes.r2)])

ylabel('r2')
title(sprintf('gtexMarkerBased r2 and SCMarkerBased r2 \n for sc markers'))
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%sboxplotOfR2ForMarkerGenes_scBasedMarkers', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% SELECTED PLOT! % quantiles of the markers
%pData = zeros(100, length(markerInds));
pData = zeros(101, 20);
for i = 1:100 
    i
    %    pData(i, :) = sort(mockRes(i).r2(markerInds));
    pData(i, :) = quantile(mockRes(i).r2(markerInds), 20);
end
pData(101, :) = quantile(scRes.r2(b(a)), 20);

% let's see at each quantile how many are bigger 
bigCounts = zeros(1, 20);
for i = 1:20
    bigCounts(i) = sum(pData(1:100, i) >= pData(101, i));
end

    
% plot the quantiles 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mockResSC
scSet
scRes

pData = zeros(101, 20);
[a, b] = ismember(gtexSet, dataSet.genes);
markerInds = b(a);
for i = 1:100 
    i
    %    pData(i, :) = sort(mockRes(i).r2(markerInds));
    pData(i, :) = quantile(mockResSC(i).r2(markerInds), 20);
end
pData(101, :) = quantile(scRes.r2(markerInds), 20);

bigCounts = zeros(1, 20);
for i = 1:20
    bigCounts(i) = sum(pData(1:100, i) >= pData(101, i));
end
bigCounts

% I should get the GTEx markers, from the SC based markers, and
% then comapre it with the mock...
h = figure
plot([1:20], pData(101,:))
hold on
ylim([-.2 1])
title('GTEx marker set predicted by SC marker set')
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%sgtexMarkersPredictedBySCMarkers', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

mockResGTEx
gtexSet
gtexRes

pData = zeros(101, 20);
[a, b] = ismember(scSet, dataSet.genes);
markerInds = b(a);
for i = 1:100 
    i
    %    pData(i, :) = sort(mockRes(i).r2(markerInds));
    pData(i, :) = quantile(mockResGTEx(i).r2(markerInds), 20);
end
pData(101, :) = quantile(gtexRes.r2(markerInds), 20);

bigCounts = zeros(1, 20);
for i = 1:20
    bigCounts(i) = sum(pData(1:100, i) >= pData(101, i));
end
bigCounts

h = figure
plot([1:20], pData(101,:))
hold on
ylim([0 .9])
title('SC marker set predicted by gtex marker set')
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%sscMarkersPredictedByGTExMarkers', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 4. getting each clusters R2s variations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
load('~/resultsAndFigures/secondProject/gtexClusters.mat') %

load('~/data/GTEx/Brain_Cortex_expGenes.mat')

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_7PCA_AllenSCBasedMarker.mat'])
scRes = result;
clear result;

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_10CellTypes_MarkersV1_cthr6.mat'])
gtexRes = result;
clear result

load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'GTExBloodClusterRep.mat']) % GTExClusterRep_Gblood.clusterInds forthe
myInds = GTExClusterRep_Gblood.clusterInds;

% getting the mat for a given cluster ID 
bs % ordered IDs from clusteringComparison

means = zeros(69, 7);
stds = zeros(69, 7);
for i = 1:length(myInds)
    mySyms = gtexCluster.syms(gtexCluster.cs1547 == myInds(i));
    [a, b] = ismember(mySyms, dataSet.genes);

    myValues = scRes.coeff(b(a), 2:end);
    means(i, :) = mean(myValues);
    stds(i, :) = std(myValues);
    % h = figure
    % errorbar(means(i, :), stds(i, :))
    % sib = gtexCluster.clusterLabels(myInds(bs(i)));
    % title(sprintf('%s', sib{1}))
    % ylim([-.03, .055])
    % figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
    % file = sprintf('%sR2plot_%s', figFolder, sib{1})
    % set(h, 'PaperOrientation', 'landscape')
    % print(h, '-deps', [file '.eps'])
    % print(h, '-dpdf', [file '.pdf'])
    % saveas(h, [file '.eps'], 'epsc')
end

sib = corr(means');
h = figure
heatmap(gtexCluster.clusterLabels(myInds(bs(1:32))), gtexCluster.clusterLabels(myInds(bs(1:32))),sib(1:32, 1:32))
colormap(hot)

boxplot(myValues)
mm = mean(myValues);
mv = std(myValues);
h = figure
errorbar(mm, mv)

% bs comes from the order of the figure, I give it in
% clusteringComparison.m

kado = means(bs(3:20), :);
sib = corr(kado');
h = figure
h = figure('units', 'centimeters', 'position', [0,0, 22, 22])
heatmap(gtexCluster.clusterLabels(myInds(bs(3:20))), ...
        gtexCluster.clusterLabels(myInds(bs(3:20))),sib, 'GridVisible', 'off')
colormap(bone)
figFolder = ['~/resultsAndFigures/secondProject/GTExRegression/']
file = sprintf('%sCTVCorr_AstroPyramidal_selected', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')



