% Here I do the QC, preprocessing and exploratory plots. 
% there are 2 files, one is read for exons, and one is read for
% introns. I do everything for both of them so sections in
% this file will have exon, intron and then "sum". 
% I do all the downstream analyses on the exon files only. 
% >>> version 3 it is 5 and 10 % of cells 
% >>> version 4 it is 2 and 5% of cells

% 0. reading the attribute file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
dataFolder = '~/data/brainSingleCell/'
file = [dataFolder 'columns-nuclei.csv']
fid = fopen(file)
linum = 1 % file has one header. 

myLine = fgets(fid);
book = strsplit(myLine, ',');

metaData = textscan(fid, [repmat('%s', 1, 56)], ...
                'Headerlines', linum, 'Delimiter', ',','Bufsize', ...
                1000000095);

for i = 3:10
    book{i}
    this = metaData{i};
    this{1}
end

% the metadata that I need: neun positive and negative, the
% sampleID 3 (8 samples), regions 6 (all MGT, skip) 7, neuron 8 (positive or
% negative), sampletype 9 (nuclei, I don't need)

mData.cellUniqueID = metaData{1};
mData.sampleID = metaData{3};
mData.regions = metaData{6};
mData.NeuN = metaData{8};
[a, b] = ismember(metaData{8}, 'NeuN-positive');
mData.NeuNPosBinary = a;

save('~/data/brainSingleCell/dataSet_meta_new.mat', 'mData')
load('~/data/brainSingleCell/dataSet_meta_new.mat')

% 1. reading the Exon file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFolder = '~/data/brainSingleCell/'
file = [dataFolder 'exons-table.csv']
fid = fopen(file)
linum = 1 % file has one header. 

myLine = fgets(fid);
book = strsplit(myLine, ',');

expData = textscan(fid, ['%s', repmat('%f', 1, 15928)], ...
                'Headerlines', linum, 'Delimiter', ',','Bufsize', ...
                1000000095);
whos expData
firstCol = expData{1};
book = expData(2:end); % first column is labels.
myExpMat = cell2mat(book);

dataSet.expMat = myExpMat;
dataSet.firstCol = expData{1};
save('~/data/brainSingleCell/dataSet_exon.mat', 'dataSet', '-v7.3')

save('~/data/brainSingleCell/dataSet_exon_corrected.mat', 'dataSet', '-v7.3')

load('~/data/brainSingleCell/dataSet_exon.mat')

% 2. reading the Intron file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFolder = '~/data/brainSingleCell/'
file = [dataFolder 'introns-table.csv']
fid = fopen(file)
linum = 1 % file has one header. 

book = strsplit(myLine, ',');

expData = textscan(fid, ['%s', repmat('%f', 1, 15928)], ...
                'Headerlines', linum, 'Delimiter', ',','Bufsize', ...
                1000000095);
whos expData
firstCol = expData{1};
book = expData(2:end); % first column is labels.
myExpMat = cell2mat(book);

dataSet.expMat = myExpMat;
dataSet.firstCol = expData{1};
geneSyms = cell(1, 50281);
for i = 1:length(geneSyms)
    temps = dataSet.firstCol{i};
    temps = temps(2:(end-1));
    geneSyms{i} = temps;
end
save('~/data/brainSingleCell/dataSet_intron.mat', 'dataSet', ['-' ...
                    'v7.3'])

% 3. Reading the read count and expressed gene count
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
figFolder = '~/data/brainSingleCell/figures/QCFigures/'
load('~/data/brainSingleCell/dataSet_meta_new.mat')

% >>> 3.1 exon
load('~/data/brainSingleCell/dataSet_exon.mat')
expMat = dataSet.expMat;

% >>>>>> 3.1.1 exon checking the Read Count (RC) per cell
rawCellCount = size(expMat, 2);

% getting the sum of reads
readSum = sum(expMat);

h = figure;
hist(readSum, 20)
title('Distribution of sample read counts')
xlabel('read count')
ylabel('cell count')
fileName = ['exon_readCount_raw'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

h = figure;
hist(readSum(mData.NeuNPosBinary), 20)
title('Distribution of sample read counts')
xlabel('read count')
ylabel('cell count')
fileName = ['exon_readCount_raw_NeuNPositive'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

h = figure;
hist(readSum(~mData.NeuNPosBinary), 20)
title('Distribution of sample read counts')
xlabel('read count')
ylabel('cell count')
fileName = ['exon_readCount_raw_NeuNNegative'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

% just checking
sum(readSum < 4e5)
sum(readSum > 2e6)

% getting MAD (for NEUN-positive samples)
med = median(readSum)
% checking for outlier cells  
MAD3 = readSum./med;
% removing samples which are not within 1/3*MAD3 and 3*MAD3
inSamples1 = MAD3 < 3;
inSamples2 = MAD3 > 1/3;
inSamples = (inSamples1 & inSamples2);
sum(inSamples)
rawCellCount - sum(inSamples1)
rawCellCount - sum(inSamples2)
rawCellCount - sum(inSamples)
exonInCellsRCPos = inSamples;

% >>>>>> 3.1.1 exon checking the gene count (GC) per cell
temp = expMat > 0;
expGeneCount = sum(temp);

h = figure;
hist(expGeneCount, 20)
title('Distribution of sample exp genes')
xlabel('exp gene count')
ylabel('cell count')
fileName = ['exon_expGeneCount_raw'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

med = median(expGeneCount)
% the 2000 thrshold is random. I do what I want! 
inSamples = expGeneCount > 2000;
exonInCellsGC = inSamples;

% >>>>>> 3.1.4 exon checking for expressed genes (EG)

geneMean = mean(expMat');
geneVar = var(expMat');

h = figure
scatter(log10(geneMean + 1), log10(geneVar+1), '.')
title('mean vs variance plot for the genes')
xlabel('log10 mean')
ylabel('log10 var')
fileName = ['exon_mean_var'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

h = figure
logGeneMean = log10(geneMean + .0003);
hist(logGeneMean, 20)
xlabel('log10 mean')
ylabel('gene count')
fileName = ['exon_log10mean'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);


% >> plotting the one gene with highest var
[a, b] = max(geneVar);
sib = expMat(b, :);
h = figure
hist(sib)

% >> plotting a few genes with high variance extremely small mean
gs1 = (geneMean < 1);
gs2 = geneVar > 2;
g3 = gs1 & gs2;
g3inds = find(g3);

% since count of cells is so high in this dataset, some genes
% expressed at small cells and low levels might miss the
% expression. I will keep genes with minimum of 10 reads in minimum
% of 100 cells as expressed
temp = expMat > 0;
cellCount = sum(temp');
hist(cellCount, 20)

%% NOTE: Filter of > 100 cells%%
pass1 = cellCount > (2 * size(expMat, 2)/100);
% now for the pass1 genes, get the mean in the samples they are
% expressed
pass1IDs = find(pass1);
pass1Means = zeros(size(pass1IDs));
for i = 1:length(pass1IDs)
    expCells = temp(pass1IDs(i), :);
    pass1Means(i) = mean(expMat(pass1IDs(i), expCells));
end
hist(pass1Means)

%% NOTE: Filter of average > 20%%
finalList = pass1IDs(pass1Means > 20);

% >>>>>>>>>> bulding final dataSet: 
geneSyms = dataSet.firstCol(finalList);
samples = exonInCellsGC & exonInCellsRC;
finalExpMat = expMat(finalList, samples);
filDataSet.expMat = finalExpMat;
filDataSet.geneSyms = geneSyms;

save('~/data/brainSingleCell/filDataSet_exon_100Cell_20Mean.mat', 'filDataSet', ['-' ...
                    'v7.3'])

% a whole other set of filters where i distinguish for NeuN pos and Neg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

neup = ismember(mData.NeuN, 'NeuN-positive');
neupExpMat = expMat(:, neup);
neunExpMat = expMat(:, ~neup);

% >>>>> FIRST: neunP : neup
% genes expressed in 2% of neup samples
temp = neupExpMat > 0;
cellCount = sum(temp');
pass1 = cellCount > (2 * size(neupExpMat, 2)/100);
% now for the pass1 genes, get the mean in the samples they are
% expressed
pass1IDs = find(pass1);
pass1Means = zeros(size(pass1IDs));
pass1LowCells = zeros(size(pass1IDs));
for i = 1:length(pass1IDs)
    expCells = temp(pass1IDs(i), :);
    thisCellCount = sum(expCells);
    if(thisCellCount < (5 * size(neupExpMat, 2)/100))
        pass1LowCells(i) = 1;
    end
    pass1Means(i) = mean(neupExpMat(pass1IDs(i), expCells));
end
hist(log10(pass1Means), 40)
[ma, mb] = max(bars)
log10(median(pass1Means))
qs = quantile(pass1Means, 3)

lowCellsFilter = pass1Means >=  qs(2);
generalCellsFilter = pass1Means >= qs(1);

% genes expressed in low count of cells but highly expressed
passedLowCells = lowCellsFilter & pass1LowCells;
passedAvgCells = generalCellsFilter & ~pass1LowCells;;
finalSet = (passedLowCells + passedAvgCells);
neupFinalGenes = pass1IDs(logical(finalSet));

% >>>>> SECOND: neunN : neun
% genes expressed in 2% of neup samples
temp = neunExpMat > 0;
cellCount = sum(temp');
pass1 = cellCount > (2 * size(neunExpMat, 2)/100);
% now for the pass1 genes, get the mean in the samples they are
% expressed
pass1IDs = find(pass1);
pass1Means = zeros(size(pass1IDs));
pass1LowCells = zeros(size(pass1IDs));
for i = 1:length(pass1IDs)
    expCells = temp(pass1IDs(i), :);
    thisCellCount = sum(expCells);
    if(thisCellCount < (5 * size(neunExpMat, 2)/100))
        pass1LowCells(i) = 1;
    end
    pass1Means(i) = mean(neunExpMat(pass1IDs(i), expCells));
end
hist(log10(pass1Means), 30)
qs = quantile(pass1Means, 3)

lowCellsFilter = pass1Means >= qs(2);
generalCellsFilter = pass1Means >= qs(1);

% genes expressed in low count of cells but highly expressed
passedLowCells = lowCellsFilter & pass1LowCells;
passedAvgCells = generalCellsFilter & ~pass1LowCells;;
finalSet = (passedLowCells + passedAvgCells);
neunFinalGenes = pass1IDs(logical(finalSet));

finalList = unique([neunFinalGenes neupFinalGenes]);

% >>>>>>>>>> bulding final dataSet: 
geneSyms = cell(size(dataSet.firstCol(finalList)));
for i = 1:length(geneSyms)
    book = dataSet.firstCol{finalList(i)};
    geneSyms{i} = book(2:end-1);
end
samples = exonInCellsGC & exonInCellsRC;
finalExpMat = expMat(finalList, samples);
filDataSet.expMat = finalExpMat;
filDataSet.geneSyms = geneSyms;

% V4 is filtered by the first and second tertiles as the filters
% for low cells filter and genreal cells filter
save('~/data/brainSingleCell/filDataSet_exon_V4.mat', 'filDataSet', ['-' ...
                    'v7.3'])

filMeta.uniqueIDs = mData.cellUniqueID(samples);
filMeta.sampleID = mData.sampleID(samples);
filMeta.regions = mData.regions(samples);
filMeta.NeuN = mData.NeuN(samples);
filMeta.NeuNPosBinary = mData.NeuNPosBinary(samples);
save('~/data/brainSingleCell/dataSet_meta_filtered_exon_V4_new.mat', ...
     'filMeta')
% >>> 3.2 intron 
clear
figFolder = '~/data/brainSingleCell/figures/QCFigures/'
load('~/data/brainSingleCell/dataSet_meta.mat')

load('~/data/brainSingleCell/dataSet_intron.mat')
expMat = dataSet.expMat;

% >>>>>> 3.1.1 intron checking the Read Count (RC) per cell
rawCellCount = size(expMat, 2);

% getting the sum of reads
readSum = sum(expMat);

h = figure;
hist(readSum, 20)
title('Distribution of sample read counts')
xlabel('read count')
ylabel('cell count')
fileName = ['intron_readCount_raw'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

% just checking
sum(readSum < 4e5)
sum(readSum > 2e6)

% getting MAD
med = median(readSum)
% checking for outlier cells  
MAD3 = readSum./med;
% removing samples which are not within 1/3*MAD3 and 3*MAD3
inSamples1 = MAD3 < 3;
inSamples2 = MAD3 > 1/3;
inSamples = (inSamples1 & inSamples2);
sum(inSamples)
rawCellCount - sum(inSamples1)
rawCellCount - sum(inSamples2)
rawCellCount - sum(inSamples)
intronInCellsRC = inSamples;

% >>>>>> 3.1.1 intron checking the gene count (GC) per cell
temp = expMat > 0;
expGeneCount = sum(temp);

h = figure;
hist(expGeneCount, 20)
title('Distribution of sample exp genes')
xlabel('exp gene count')
ylabel('cell count')
fileName = ['intron_expGeneCount_raw'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

med = median(expGeneCount)
% the 2000 thrshold is random. I do what I want! 
inSamples = expGeneCount > 2000;
intronInCellsGC = inSamples;

% >>>>>> 3.1.4 exon checking for expressed genes (EG)

geneMean = mean(expMat');
geneVar = var(expMat');

h = figure
scatter(log10(geneMean + 1), log10(geneVar+1), '.')
title('mean vs variance plot for the genes')
xlabel('log10 mean')
ylabel('log10 var')
fileName = ['intron_mean_var'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

h = figure
logGeneMean = log10(geneMean + .0003);
hist(logGeneMean, 20)
xlabel('log10 mean')
ylabel('gene count')
fileName = ['intron_log10mean'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

% >> plotting the one gene with highest var
[a, b] = max(geneVar);
sib = expMat(b, :);
h = figure
hist(sib)

% >> plotting a few genes with high variance extremely small mean
gs1 = (geneMean < 1);
gs2 = geneVar > 2;
g3 = gs1 & gs2;
g3inds = find(g3);

% since count of cells is so high in this dataset, some genes
% expressed at small cells and low levels might miss the
% expression. I will keep genes with minimum of 10 reads in minimum
% of 100 cells as expressed
temp = expMat > 0;
cellCount = sum(temp');
hist(cellCount, 20)

% skipping this for intron
%% NOTE: Filter of > 100 cells%%
pass1 = cellCount > (5 * size(expMat, 2)/100);
% now for the pass1 genes, get the mean in the samples they are
% expressed
pass1IDs = find(pass1);
pass1Means = zeros(size(pass1IDs));
for i = 1:length(pass1IDs)
    expCells = temp(pass1IDs(i), :);
    pass1Means(i) = mean(expMat(pass1IDs(i), expCells));
end
hist(pass1Means)

%% NOTE: Filter of average > 20%%
finalList = pass1IDs(pass1Means > 20);

% >>>>>>>>>> bulding final dataSet: 
geneSyms = dataSet.firstCol(finalList);
samples = intronInCellsGC & intronInCellsRC;
finalExpMat = expMat(finalList, samples);
filDataSet.expMat = finalExpMat;
filDataSet.geneSyms = geneSyms;

save('~/data/brainSingleCell/filDataSet_intron_100Cell_20Mean.mat', 'filDataSet', ['-' ...
                    'v7.3'])

% a whole other set of filters where i distinguish for NeuN pos and Neg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

neup = ismember(mData.NeuN, 'NeuN-positive');
neupExpMat = expMat(:, neup);
neunExpMat = expMat(:, ~neup);

% >>>>> FIRST: neunP : neup
% genes expressed in 2% of neup samples
temp = neupExpMat > 0;
cellCount = sum(temp');
pass1 = cellCount > (2 * size(neupExpMat, 2)/100);
% now for the pass1 genes, get the mean in the samples they are
% expressed
pass1IDs = find(pass1);
pass1Means = zeros(size(pass1IDs));
pass1LowCells = zeros(size(pass1IDs));
for i = 1:length(pass1IDs)
    expCells = temp(pass1IDs(i), :);
    thisCellCount = sum(expCells);
    if(thisCellCount < (5 * size(neupExpMat, 2)/100))
        pass1LowCells(i) = 1;
    end
    pass1Means(i) = mean(neupExpMat(pass1IDs(i), expCells));
end
hist(log10(pass1Means), 30)
qs = quantile(pass1Means, 3)

lowCellsFilter = pass1Means >= qs(2);
generalCellsFilter = pass1Means >= qs(1);

% genes expressed in low count of cells but highly expressed
passedLowCells = lowCellsFilter & pass1LowCells;
passedAvgCells = generalCellsFilter & ~pass1LowCells;;
finalSet = (passedLowCells + passedAvgCells);
neupFinalGenes = pass1IDs(logical(finalSet));

% >>>>> SECOND: neunN : neun
% genes expressed in 2% of neup samples
temp = neunExpMat > 0;
cellCount = sum(temp');
pass1 = cellCount > (2 * size(neunExpMat, 2)/100);
% now for the pass1 genes, get the mean in the samples they are
% expressed
pass1IDs = find(pass1);
pass1Means = zeros(size(pass1IDs));
pass1LowCells = zeros(size(pass1IDs));
for i = 1:length(pass1IDs)
    expCells = temp(pass1IDs(i), :);
    thisCellCount = sum(expCells);
    if(thisCellCount < (5 * size(neunExpMat, 2)/100))
        pass1LowCells(i) = 1;
    end
    pass1Means(i) = mean(neunExpMat(pass1IDs(i), expCells));
end
hist(log10(pass1Means), 30)
qs = quantile(pass1Means , 3)

lowCellsFilter = pass1Means >= qs(2);
generalCellsFilter = pass1Means >= qs(1);

% genes expressed in low count of cells but highly expressed
passedLowCells = lowCellsFilter & pass1LowCells;
passedAvgCells = generalCellsFilter & ~pass1LowCells;;
finalSet = (passedLowCells + passedAvgCells);
neunFinalGenes = pass1IDs(logical(finalSet));

finalList = unique([neunFinalGenes neupFinalGenes]);

% >>>>>>>>>> bulding final dataSet: 
geneSyms = dataSet.firstCol(finalList);
samples = intronInCellsGC & intronInCellsRC;
finalExpMat = expMat(finalList, samples);
filDataSet.expMat = finalExpMat;
filDataSet.geneSyms = geneSyms;

save('~/data/brainSingleCell/filDataSet_intron_V4.mat', 'filDataSet', ['-' ...
                    'v7.3'])

filMeta.sampleID = mData.sampleID(samples);
filMeta.regions = mData.regions(samples);
filMeta.NeuN = mData.NeuN(samples);
filMeta.NeuNPosBinary = mData.NeuNPosBinary(samples);
save('~/data/brainSingleCell/dataSet_meta_filtered_intron_V4.mat', ...
     'filMeta')

% >>> 3.3 sum 

clear
figFolder = '~/data/brainSingleCell/figures/QCFigures/'
load('~/data/brainSingleCell/dataSet_meta.mat')

load('~/data/brainSingleCell/dataSet_intron.mat')
expMat1 = dataSet.expMat;

load('~/data/brainSingleCell/dataSet_exon.mat')
expMat2 = dataSet.expMat;

expMat = expMat1 + expMat2;
clear expMat1 expMat2

% >>>>>> 3.1.1 sum - checking the Read Count (RC) per cell
rawCellCount = size(expMat, 2);

% getting the sum of reads
readSum = sum(expMat);

h = figure;
hist(readSum, 20)
title('Distribution of sample read counts')
xlabel('read count')
ylabel('cell count')
fileName = ['sum_readCount_raw'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

% just checking
sum(readSum < 10e5)
sum(readSum > 4e6)

% getting MAD
med = median(readSum)
% checking for outlier cells  
MAD3 = readSum./med;
% removing samples which are not within 1/3*MAD3 and 3*MAD3
inSamples1 = MAD3 < 3;
inSamples2 = MAD3 > 1/3;
inSamples = (inSamples1 & inSamples2);
sum(inSamples)
rawCellCount - sum(inSamples1)
rawCellCount - sum(inSamples2)
rawCellCount - sum(inSamples)
sumInCellsRC = inSamples;

% >>>>>> 3.1.1 sum - checking the gene count (GC) per cell
temp = expMat > 0;
expGeneCount = sum(temp);

h = figure;
hist(expGeneCount, 20)
title('Distribution of sample exp genes')
xlabel('exp gene count')
ylabel('cell count')
fileName = ['sum_expGeneCount_raw'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

med = median(expGeneCount)
% the 2000 thrshold is random. I do what I want! 
inSamples = expGeneCount > 4000;
sumInCellsGC = inSamples;

% >>>>>> 3.1.4 intron checking for expressed genes (EG)
geneMean = mean(expMat');
geneVar = var(expMat');

h = figure
scatter(log10(geneMean + 1), log10(geneVar+1), '.')
title('mean vs variance plot for the genes')
xlabel('log10 mean')
ylabel('log10 var')
fileName = ['sum_mean_var'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

h = figure
logGeneMean = log10(geneMean + .0003);
hist(logGeneMean, 20)
xlabel('log10 mean')
ylabel('gene count')
fileName = ['sum_log10mean'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

% >> plotting the one gene with highest var
[a, b] = max(geneVar);
sib = expMat(b, :);
h = figure
hist(sib)

% >> plotting a few genes with high variance extremely small mean
gs1 = (geneMean < 1);
gs2 = geneVar > 2;
g3 = gs1 & gs2;
g3inds = find(g3);

% since count of cells is so high in this dataset, some genes
% expressed at small cells and low levels might miss the
% expression. I will keep genes with minimum of 10 reads in minimum
% of 100 cells as expressed
temp = expMat > 0;
cellCount = sum(temp');
hist(cellCount, 20)

% skipping this
% %% NOTE: Filter of > 100 cells%%
% pass1 = cellCount > (5 * size(expMat, 2)/100);
% % now for the pass1 genes, get the mean in the samples they are
% % expressed
% pass1IDs = find(pass1);
% pass1Means = zeros(size(pass1IDs));
% for i = 1:length(pass1IDs)
%     expCells = temp(pass1IDs(i), :);
%     pass1Means(i) = mean(expMat(pass1IDs(i), expCells));
% end
% hist(pass1Means)

% %% NOTE: Filter of average > 20%%
% finalList = pass1IDs(pass1Means > 20);

% % >>>>>>>>>> bulding final dataSet: 
% geneSyms = dataSet.firstCol(finalList);
% samples = intronInCellsGC & intronInCellsRC;
% finalExpMat = expMat(finalList, samples);
% filDataSet.expMat = finalExpMat;
% filDataSet.geneSyms = geneSyms;

% save('~/data/brainSingleCell/filDataSet_intron_100Cell_20Mean.mat', 'filDataSet', ['-' ...
%                     'v7.3'])

% a whole other set of filters where i distinguish for NeuN pos and Neg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

neup = ismember(mData.NeuN, 'NeuN-positive');
neupExpMat = expMat(:, neup);
neunExpMat = expMat(:, ~neup);

% >>>>> FIRST: neunP : neup
% genes expressed in 2% of neup samples
temp = neupExpMat > 0;
cellCount = sum(temp');
pass1 = cellCount > (2 * size(neupExpMat, 2)/100);
% now for the pass1 genes, get the mean in the samples they are
% expressed
pass1IDs = find(pass1);
pass1Means = zeros(size(pass1IDs));
pass1LowCells = zeros(size(pass1IDs));
for i = 1:length(pass1IDs)
    expCells = temp(pass1IDs(i), :);
    thisCellCount = sum(expCells);
    if(thisCellCount < (5 * size(neupExpMat, 2)/100))
        pass1LowCells(i) = 1;
    end
    pass1Means(i) = mean(neupExpMat(pass1IDs(i), expCells));
end
h = figure
hist(log10(pass1Means), 30)
qs = quantile(pass1Means, 3)
 
lowCellsFilter = pass1Means >= qs(2);
generalCellsFilter = pass1Means >= qs(1);

% genes expressed in low count of cells but highly expressed
passedLowCells = lowCellsFilter & pass1LowCells;
passedAvgCells = generalCellsFilter & ~pass1LowCells;
finalSet = (passedLowCells + passedAvgCells);
neupFinalGenes = pass1IDs(logical(finalSet));

% >>>>> SECOND: neunN : neun
% genes expressed in 2% of neun samples
temp = neunExpMat > 0;
cellCount = sum(temp');
pass1 = cellCount > (2 * size(neunExpMat, 2)/100);
% now for the pass1 genes, get the mean in the samples they are
% expressed
pass1IDs = find(pass1);
pass1Means = zeros(size(pass1IDs));
pass1LowCells = zeros(size(pass1IDs));
for i = 1:length(pass1IDs)
    expCells = temp(pass1IDs(i), :);
    thisCellCount = sum(expCells);
    if(thisCellCount < (5 * size(neunExpMat, 2)/100))
        pass1LowCells(i) = 1;
    end
    pass1Means(i) = mean(neunExpMat(pass1IDs(i), expCells));
end
h = figure
hist(log10(pass1Means), 30)
qs = quantile(pass1Means, 3)

lowCellsFilter = pass1Means >= qs(2);
generalCellsFilter = pass1Means >= qs(1);

% genes expressed in low count of cells but highly expressed
passedLowCells = lowCellsFilter & pass1LowCells;
passedAvgCells = generalCellsFilter & ~pass1LowCells;
finalSet = (passedLowCells + passedAvgCells);
neunFinalGenes = pass1IDs(logical(finalSet));

finalList = unique([neunFinalGenes neupFinalGenes]);

% >>>>>>>>>> bulding final dataSet: 
geneSyms = dataSet.firstCol(finalList);
samples = sumInCellsGC & sumInCellsRC;
finalExpMat = expMat(finalList, samples);
filDataSet.expMat = finalExpMat;
tempSym;
for i = 1:length(geneSyms)
    tempStr = geneSyms{i};
    tempSym{i} = tempStr(2:end-1);
end

filDataSet.geneSyms = tempSym;

save('~/data/brainSingleCell/filDataSet_intronAndExon_V4.mat', 'filDataSet', ['-' ...
                    'v7.3'])

load('~/data/brainSingleCell/filDataSet_intronAndExon_V4.mat')

filMeta.cellUniqueID = mData.cellUniqueID(samples);
filMeta.sampleID = mData.sampleID(samples);
filMeta.regions = mData.regions(samples);
filMeta.NeuN = mData.NeuN(samples);
filMeta.NeuNPosBinary = mData.NeuNPosBinary(samples);
save('~/data/brainSingleCell/dataSet_meta_filtered_exonAndIntron_V4.mat', ...
     'filMeta')

% Difference between V2 and V3 is the average count used for the
% filter. This could be changed. 
