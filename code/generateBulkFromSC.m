% this is the file to simulate bulk tissue dataset, build the
% network from it and compare it with the bulk tissue coexpression.

% 3. get the sum simulated network

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

% get the inds for each cell type:
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
[a, b] = sort(clusterMeta.sortedClusterNames);
b(1:2) % Astro
b(3) % endo
b(4:27) % Exc - pyra
b(28:72) % Inh
b(73) % micro
b(74) % OPC
b(75) % oligo

% Pyra
[c, d] = ismember(clusterMeta.clusters, ...
                  clusterMeta.sortedClusterNames(b(4:27)));
sum(c)
inds{1} = find(c);

% Inh
[c, d] = ismember(clusterMeta.clusters, ...
                  clusterMeta.sortedClusterNames(b(28:72)));
sum(c)
inds{2} = find(c);

% Oligo
[c, d] = ismember(clusterMeta.clusters, ...
                  clusterMeta.sortedClusterNames(b(75)));
sum(c)
inds{3} = find(c);

% Astro
[c, d] = ismember(clusterMeta.clusters, ...
                  clusterMeta.sortedClusterNames(b(1:2)));
sum(c)
inds{4} = find(c);

% Micro
[c, d] = ismember(clusterMeta.clusters, ...
                  clusterMeta.sortedClusterNames(b(73)));
sum(c)
inds{5} = find(c);

% OPC
[c, d] = ismember(clusterMeta.clusters, ...
                  clusterMeta.sortedClusterNames(b(74)));
sum(c)
inds{6} = find(c);


% Endo
[c, d] = ismember(clusterMeta.clusters, ...
                  clusterMeta.sortedClusterNames(b(3)));
sum(c)
inds{7} = find(c);

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Building a dataset of 100 samples, where each sample is sum of 200
% samples. The ratio between non neuron and neuron should have the
% average of 50 50, the ratio between oligo, astro and micro should
% be 75, 20, 5. So in total: 20% Pyra, 20% Inh, 45%Olig 12% As,
% 3%micro

%% first way of making networks
for f = 1:65
    f
    % with 100 samples with sigma 3
    counts = zeros(100, 5);
    simBulkExpMat = zeros(size(expMat, 1), 100);
    selectedSamples = cell(1, 100);
    
    for i = 1:100
        
        % >>>>> with normal distribution
        % thisCounts = normrnd([40.2 40.2 90.2 24.2 6.2], 3);
        %thisCounts = normrnd([40 40 80 24 10], [40 40 80 24,10]./3);
        
        %        thisCounts = floor(max([thisCounts; 0 0 0 0 0 ]))
        
        % >>>>> with no Variation
        
        thisCounts = [40 40 80 24 10];
        % >>>>>>>>>>
        
        counts(i, :) = thisCounts;

        tempInds = zeros(1, size(normExp, 2));
        thisCounts(thisCounts < 0) = 0;
        for j = 1:5
            localCount = floor(min(thisCounts(j), length(inds{j})));
            dso = datasample(inds{j}, localCount);
            tempMat = normExp(:, dso);
            simBulkExpMat(:, i) = simBulkExpMat(:, i) + sum(tempMat, 2);
            tempInds(dso) = 1;
        end
        selectedSamples{i} = find(tempInds);
    end

    % >>> normalizing simBulkExpMat
    sumExp = sum(simBulkExpMat);
    milFac = sumExp ./ 1000000;
    normSimMat = zeros(size(simBulkExpMat));
    sampleCount = size(simBulkExpMat, 2);
    for i = 1:sampleCount
        normSimMat(:, i) = simBulkExpMat(:, i)./milFac(i);
    end

    % >>>> corr with p
    % kado = sum(normSimMat' > 0);
    % fixedMat = normSimMat;
    % fixedMat(normSimMat == 0) = nan;
    % tic
    % [sib1, p1] = corr(fixedMat', 'rows', ...
    %                   'pairwise', 'Tail', 'right');
    % toc
    % upP1 = p1(logical(triu(ones(size(p1)), 1)));
    % qs = quantile(upP1, 1000);

    % myBinNet = zeros(size(p1));
    % myBinNet(p1 < qs(10)) = 1;
    % sum(myBinNet(:))
    % <<<<<<<< 
    
    sib = corr(normSimMat');
    upSib = sib(logical(triu(ones(size(sib)), 1)));
    qs = quantile(upSib, 1000);
    myBinNet01 = (sib > qs(985)) - eye(size(sib));
    myBinNet005 = (sib > qs(995)) - eye(size(sib));

    bulkFromSC.counts = counts;
    bulkFromSC.simBulkExpMat = simBulkExpMat;
    bulkFromSC.selectedSamples = selectedSamples;
    bulkFromSC.binNet005 = sparse(myBinNet005);
    bulkFromSC.binNet01 = sparse(myBinNet01);
    %    bulkFromSC.normrndInput = [40 40 40 40 10; 15 15 15 15 5];
        bulkFromSC.normrndInput = [[40 40 80 24 10], [0 0 0 0 0]./3];
    bulkFromSC.geneSyms = filDataSet.geneSyms;
    save(sprintf('~/resultsAndFigures/secondProject/SimBulkNetworksFromSC/bulkFromSC_%d_Var0.mat', f), ...
         'bulkFromSC')
end

%% second way of making networks
for f = 1:20
    % with 100 samples with sigma 3
    counts = zeros(100, 5);
    simBulkExpMat = zeros(size(expMat, 1), 100);
    selectedSamples = cell(1, 100);
    %  baseCounts = [80 80 160 48 20];
    baseCounts = [40 40 40 40 20];
    for i = 1:100
        i
        baseInds = zeros(1, sum(baseCounts));
        
        % dso = datasample(inds{1}, baseCounts(1));
        % baseInds(1:80) = dso;
        
        % dso = datasample(inds{2}, baseCounts(2));
        % baseInds(81:160) = dso;
        
        % dso = datasample(inds{3}, baseCounts(3));
        % baseInds(161:320) = dso;
        
        % dso = datasample(inds{4}, baseCounts(4));
        % baseInds(321:368) = dso;
        
        % baseInds(369:388) = inds{5};
        
        dso = datasample(inds{1}, baseCounts(1));
        baseInds(1:40) = dso;
        
        dso = datasample(inds{2}, baseCounts(2));
        baseInds(41:80) = dso;
        
        dso = datasample(inds{3}, baseCounts(3));
        baseInds(81:120) = dso;
        
        dso = datasample(inds{4}, baseCounts(4));
        baseInds(121:160) = dso;
        
        baseInds(161:180) = inds{5};

        
        tempInds = datasample(baseInds, 100);
        
        selectedSamples{i} = tempInds;
        simBulkExpMat(:, i) = sum(normExp(:, tempInds)');
    end

    % >>> normalizing simBulkExpMat
    sumExp = sum(simBulkExpMat);
    milFac = sumExp ./ 1000000;
    normSimMat = zeros(size(simBulkExpMat));
    sampleCount = size(simBulkExpMat, 2);
    for i = 1:sampleCount
        normSimMat(:, i) = simBulkExpMat(:, i)./milFac(i);
    end

    % >>>> corr with p
    % kado = sum(normSimMat' > 0);
    % fixedMat = normSimMat;
    % fixedMat(normSimMat == 0) = nan;
    % tic
    % [sib1, p1] = corr(fixedMat', 'rows', ...
    %                   'pairwise', 'Tail', 'right');
    % toc
    % upP1 = p1(logical(triu(ones(size(p1)), 1)));
    % qs = quantile(upP1, 1000);

    % myBinNet = zeros(size(p1));
    % myBinNet(p1 < qs(10)) = 1;
    % sum(myBinNet(:))
    % <<<<<<<< 
    
    sib = corr(normSimMat');
    upSib = sib(logical(triu(ones(size(sib)), 1)));
    qs = quantile(upSib, 1000);
    myBinNet01 = (sib > qs(985)) - eye(size(sib));
    myBinNet005 = (sib > qs(995)) - eye(size(sib));

    %bulkFromSC.counts = counts;
    bulkFromSC.simBulkExpMat = simBulkExpMat;
    bulkFromSC.selectedSamples = selectedSamples;
    bulkFromSC.binNet005 = sparse(myBinNet005);
    bulkFromSC.binNet01 = sparse(myBinNet01);
    %    bulkFromSC.normrndInput = [40 40 40 40 10; 15 15 15 15 5];
    %       bulkFromSC.normrndInput = [[40 40 80 24 10], [40 40 80 24 10]/1.2];
    bulkFromSC.geneSyms = filDataSet.geneSyms;
    save(sprintf('~/resultsAndFigures/secondProject/SimBulkNetworksFromSC/bulkFromSC_%d_newComb6_from180.mat', f), ...
         'bulkFromSC')
end


sumNet = zeros(16789, 16789);
olcs = zeros(1, 10);
for f = 1:10
    load(sprintf('~/resultsAndFigures/secondProject/SimBulkNetworksFromSC/bulkFromSC_%d.mat', ...
                 f))
    sumNet = sumNet + bulkFromSC.binNet005;
    %    thison = gtexFullNet + bulkFromSC.binNet005(a, a);
    %olcs(f) = sum(thison(:) == 2);
end

aggNet = sumNet > 5;

scSyms = bulkFromSC.geneSyms;
[a, b] = ismember(scSyms, thisNetSyms); % thisNetSyms = gtex
                                        % network symbols
aggNetGTEx = aggNet(a, a);
binNetGTEx = myBinNet005(a, a);

net1 = aggNetGTEx;
net1 = binNetGTEx;
net1 = scNet;
net2 = gtexFullNet;
overlapNet = net1 + net2;
n = 13807;
n = 9544
d1 = sum(net1(:)) / (n*(n-1))
d2 = sum(net2(:)) / (n*(n-1))
do = sum(overlapNet(:) == 2) / (n*(n-1))

do / (d1 * d2)

nds = sum(myBinNet005);

% how does the correlation of markers look? 
load('~/resultsAndFigures/secondProject/scBasedMarkers.mat')

[a, b] = ismember(scBasedMarkers.sortedGenes, filDataSet.geneSyms);
ctInds = [1:179, 311:437 471:792]; % astor, micro, oligo , pyramidal
[a, b] = ismember(scBasedMarkers.sortedGenes(ctInds), ...
                  filDataSet.geneSyms);
[a, b] = ismember(scBasedMarkers.sortedGenes(ctInds), ...
                  testGeneList);

whos a
sum(a)
geneList = scBasedMarkers.sortedGenes(ctInds);
geneList = filDataSet.geneSyms;
testGeneList = geneList(a);
kado = sib(a, a); % the a comes from the clustersInNetworks.m, a
                  % are the genes which are present in GTEx network

kado = sib(b, b);
kado = bulkFromSC.binNet005(b, b);
kado = aggNet(b, b);
kado = kado - eye(length(kado));
h = figure
heatmap(kado+0)
hist(kado(:))

qs = quantile(kado(:), 1000);

myBinNet = zeros(size(kado));
myBinNet(kado > qs(985)) = 1;
sum(myBinNet(:))

nds = sum(myBinNet);
sum(nds > 0)

h = figure
heatmap(myBinNet(nds>0, nds>0))

h = figure
subplot(2, 3, 1)
plot(sort(nds(1:179)))
title('Astro')

subplot(2, 3, 2)
plot(sort(nds(179:311)))
title('Endo')

subplot(2, 3, 3)
plot(sort(nds(311:437)))
title('Micro')

subplot(2, 3, 4)
plot(sort(nds(437:471)))
title('OPC')

subplot(2, 3, 5)
plot(sort(nds(471:612)))
title('Oligo')

subplot(2, 3, 6)
plot(sort(nds(612:792))) 
title('Pyramidal')

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% here I continue: I found a set of genes are correlated wit heach
% other in the myBinNet005, disproportionally. They are masking the
% rest of the observed coexpression. Based on the ND, I found
% them. 

% the network is almost empty except for a very few sets of genes. 
nds = sum(myBinNet005);
geneList = filDataSet.geneSyms;
subGeneList = geneList(nds > 1);
subBin = myBinNet005(nds >1, nds>1);
subNDs = sum(subBin);
hubInds = find(subNDs > 1600);
h = figure
heatmap(subGeneList(hubInds), subGeneList(hubInds), subBin(hubInds, ...
                                                  hubInds))
localGeneList = subGeneList(hubInds);

whos nds
invGenes = find(nds > 1600);

% normSimMat is the simulated matrix (where coexpression of
% myBinNet005 has come from)

h = figure
heatmap(log2(normSimMat(invGenes, :)+1))


h = figure
heatmap(log2(simBulkExpMat(invGenes, :)+1))


h = figure
heatmap(log2(normExp(invGenes, selectedSamples{9})+1))

mySum9 = sum(normExp(invGenes, selectedSamples{9})');
h = figure
heatmap(mySum9)

h = figure
heatmap(normExp(invGenes, selectedSamples{11}))


% it is clear that samples 60 and 93 and these bunch of genes are
% problematic. For now, I forget where they are coming from (which
% samples) and I just remove the samples 60 and 93 from my
% simulated data. 

normSimMatMinus = normSimMat;

normSimMatMinus(:, [60, 93]) = [];
h = figure
heatmap(log2(normSimMatMinus(invGenes, :)+1))

% now I will build the networks without them
sibm = corr(normSimMatMinus');
upSibm = sibm(logical(triu(ones(size(sibm)), 1)));
qsm = quantile(upSibm, 1000);
myBinNet01m = (sibm > qsm(990)) - eye(size(sibm));
myBinNet005m = (sibm > qsm(995)) - eye(size(sibm));

sum(sum(myBinNet005m(a, a)))

subBinNet005m = myBinNet005m(a, a);
subBinNet01m = myBinNet01m(a, a);

halva = subBinNet01m + fullNet;

% now that we got the networks going, let's find those genes!

whos normExp

invExp = normExp(invGenes, :);
mySum = sum(invExp);
myCells = find(mySum > 800);

plotList=clusterMeta.clusters(myCells);
[pa, pb] = sort(plotList);
modPlotLabs = cell(1, length(pa));
for i = 1:length(pa)
    modPlotLab{i} = sprintf('%d %s', i, pa{i});
end

heatmap(modPlotLab, geneList(invGenes), log2(normExp(invGenes, pb)+1))

trouble.simExpMat = normSimMat;
trouble.problemExp = invGenes;
save('~/resultsAndFigures/secondProject/normSimMat_Trouble.mat', 'trouble')

% the added network is about 7.5 times more reproducible in GTEx
% compared to 3.2 times.

% Does it meet the reproducibility between Affy and GTEx?
% the rep between affy and GTEx is 29

% the marker block stays perfect. 

% get the sum networks of simulated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

files = dir(['~/resultsAndFigures/secondProject/' ...
             'SimBulkNetworksFromSC/*Var0.mat'])

files = dir(['~/resultsAndFigures/secondProject/' ...
             'SimBulkNetworksFromSC/*Var33.mat'])

files = dir(['~/resultsAndFigures/secondProject/' ...
             'SimBulkNetworksFromSC/*newComb8*.mat'])
        load([files(1).folder '/' files(1).name])
        bulkFromSC

files = dir(['~/resultsAndFigures/secondProject/' ...
             'SimBulkNetworksFromSC/withCounts/*.mat'])

sumNet = zeros(16789, 16789);
for j = 1:length(files)
    j
    %    load(sprintf('%s/bulkFromSC_%d.mat', files(1).folder, j))
        load([files(j).folder '/' files(j).name])
        sumNet = sumNet + bulkFromSC.binNet005;
    % [aj, bj] = ismember({'NRGN', 'CALM3'}, bulkFromSC.geneSyms);
    % tinyNet = bulkFromSC.binNet005(bj, bj);
    % tsNet = tsNet + tinyNet;
end

myNet = sparse(sumNet);
sumNet = myNet;
save(['~/resultsAndFigures/secondProject/SimBulkNetworksFromSC/' ...
      'sumVarNet_33Var_65Net.mat'], 'sumNet')

[a, b] = ismember(mySyms, bulkFromSC.geneSyms)
h = figure
heatmap(log2(sumNet(b(a), b(a))))

[a, b] = ismember({'NRGN', 'CALM3'}, bulkFromSC.geneSyms)

sumNet(b, b)
[a, b] = ismember({'NRGN', 'CALM3'}, dataSet.genes)
result.r2(b)
