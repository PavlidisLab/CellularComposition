% this file has different parts. First is processing the microarray
% datasets
% 00. main files and roots
% 0. fix colormap
% 1. microarray processing
% 2. getting list of DEG from the paper
% 3. get the bulk tissue data
% 4. identify the clusters (code from clusteringGTExNetworks.m)
% 5. expression of the markers in different clusters 
% 6. regression
% 7. Result 2: study of the R2 values for different clusters
% 8. cluster functional enrichment - brain comparison
% 9. gene expression levels in different clusters versus celltypes
% 10. print clusters for cytoscape
% 11. cluster the clusters and sort them out and finish
% 11.5 get the coorelation of CT profiles and expression in celltypes
% 12. working with diff genes and their expression AND filtering
% markers based on coexpression
% 13. Blood TAN work - Cluster bloodtan network
% 14. expression level of genes in datasets
% 15. blood mouse data for Erythrocyte markers 
% 16. first report plot
% 17. affy blood dataset human bulk
% 18. hk genes in the clusters
% 19. that regression thing : are the Rs from markers more
% correlated than Rs from random groups of genes. 
% 20. affy blood dataset human bulk GSE27562
% 21. correlation of different celltypes for clusters
% 22. correlation of R2, CT profiles and variance for affy for
% blood genes
% 23. plot for distribution of R2
% 24. R2 for clusters
% 25. writing supplemental files

% D1. testing exp level of diff genes

% TODO: add the Erythrocyte markers and then cluster
% do the average variance and average R2 for clusters
% do the r2 for the new markers for each cluster
% do the functional enrichment - report etc.

%  00. main files and roots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataFolder = '/home/mfarahbod/data/affyArray/bloodCellTypeDataset/'

% load GTEx blood network, clustering and gene symbols, outperm and
% clusterOverlap
load('~/resultsAndFigures/secondProject/bloodDS.mat') % bloodDS

load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat')
bloodNet = GTExFiveNets.nets(1);
bloodExpGenes = bloodDS.genes(bloodNet.expGenes);

load('~/resultsAndFigures/secondProject/gtexClusters_blood.mat') % gtexCluster
cs = gtexCluster.cs1547;
clusterMemberCount = hist(cs, unique(cs))
inClusterIds = find(clusterMemberCount >= 20);

% cluster blood cluster
% cinfo.inIDs = inClusterIds;
% cinfo.overlapDensity = overlapDensity;
% cinfo.outperm = outperm;
% cinfo.clusterLabels = clusterLabels;
% save('~/resultsAndFigures/secondProject/bloodClusterInfo.mat', 'cinfo')
load('~/resultsAndFigures/secondProject/bloodClusterInfo.mat') %
                                                               % cinfo 
selectedCs = cinfo.outperm([1:54, 57, 59, 61]); % removing clusters
                                                % with low
                                                % expression levels

load('~/resultsAndFigures/secondProject/bloodDiffGenesInfo.mat') 
diff

% load GTEx brain cluster enrichment 

brainRes = load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_AllenSCBasedMarker.mat'])

load('~/data/GTEx/Brain_Cortex_expGenes.mat')
brainDS.mat = log2(dataSet.mat+1);
brainDS.genes = dataSet.genes;

[a, b] = ismember(gtexCluster.syms, brainDS.genes);
brr = brainRes.result.r2(b(a));
[a, b] = ismember(bloodExpGenes, gtexCluster.syms(a));
blr = regRes.r2(a);
corr(blr, brr) % this corr is .11 - compared to blood datasets
               % which is .37 it is interesting

% 0. fix colormap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = colormap(1- [1-gray(100)' copper(100)']');

book = temp(1:(length(temp)/2), :)';
book = [book ones(3, 40) temp(100:end, :)']';

colormap(book)

% >>>> the four colormaps
bcmap = [247,247,247 % black
204,204,204
150,150,150
82,82,82]

lcmap = [239,243,255 % blue
189,215,231
107,174,214
33,113,181]

rcmap = [254,229,217 % red
252,174,145
251,106,74
203,24,29] 

pcmap = [242,240,247 % purple
203,201,226
158,154,200
106,81,163]

ocmap = [254,237,222  % orange
253,190,133
253,141,60
217,71,1]


% 1. microarray processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1.1 preprocess the .cel files
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

celFolder = '/home/mfarahbod/data/affyArray/bloodCellTypeDataset/celFiles'
celFolder = '/space/grp/marjan/data/GSE22886/B'
tarList = dir(celFolder)
dataFolder = '/home/mfarahbod/data/affyArray/bloodCellTypeDataset/'

name = 'GSE22886_RAW'
system(['gunzip ' celFolder '/*.gz'], '-echo');

% MATLAB can not pass tilde right, it has problem
libFile = '~/data/affyArray/bloodCellTypeDataset/HG_U133A.CDF'
libFile = ['~/data/affyArray/bloodCellTypeDataset/HG-' ...
           'U133_Plus_2.CDF'];
libFile = '/space/grp/marjan/data/HG-U133_Plus_2.CDF';
libFile = '/space/grp/marjan/data/HG-U133B.CDF';

% 1.2 getting list of the cel files here
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

tempCelFileList = dir(fullfile(celFolder ,'*.CEL'));
celFileList = cell(1, length(tempCelFileList));
for j = 1:length(tempCelFileList)
    celFileList{j} = tempCelFileList(j).name;
end

expr = affyrma(celFileList, libFile, 'CELPath', ...
               celFolder);
save('~/data/affyArray/bloodCellTypeDataset/GSE22886_A.mat', 'expr')
save('~/data/affyArray/bloodCellTypeDataset/GSE22886_B.mat', ...
     'expr')

clear expr
load(['/home/mfarahbod/data/affyArray/bloodCellTypeDataset/' ...
      'GSE22886_A.mat'])
dmwrite(expr, ['/home/mfarahbod/data/affyArray/' ...
               'bloodCellTypeDataset/GSE22886_A.txt'])

clear expr
load(['/home/mfarahbod/data/affyArray/bloodCellTypeDataset/' ...
      'GSE22886_B.mat'])
dmwrite(expr, ['/home/mfarahbod/data/affyArray/' ...
               'bloodCellTypeDataset/GSE22886_B.txt'])

% 1.3 getting the platform files
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fid = fopen(['/home/mfarahbod/data/affyArray/bloodCellTypeDataset/' ...
             'GPL96_noParents.an.txt'])
gpl = textscan(fid, [repmat('%s', 1, 11)], 'Delimiter', '\t', ...
               'Headerlines',8);

%GPL96 probe ID's
gpl96.probeID = gpl{1};

%GPL96 gene Symbols
gpl96.geneSymbol = gpl{2};
gpl96.uniqueGeneSymbols = unique(gpl{2});

%GPL96 gene names
gpl96.geneNames = gpl{3};

save([dataFolder, 'gpl96.mat'], 'gpl96')

fid = fopen(['/home/mfarahbod/data/affyArray/bloodCellTypeDataset/' ...
             'GPL97_noParents.an.txt'])
gpl = textscan(fid, [repmat('%s', 1, 11)], 'Delimiter', '\t', ...
               'Headerlines',8);

%GPL97 probe ID's
gpl97.probeID = gpl{1};

%GPL97 gene Symbols
gpl97.geneSymbol = gpl{2};
gpl97.uniqueGeneSymbols = unique(gpl{2});

%GPL97 gene names
gpl97.geneNames = gpl{3};

save([dataFolder, 'gpl97.mat'], 'gpl97')

% 1.4 read the .txt files and build the gene datasets
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% >>>>>>>>> for samples A
fileName = ['/home/mfarahbod/data/affyArray/' ...
               'bloodCellTypeDataset/GSE22886_A.txt'];
fid = fopen(fileName)
linum = 1
header = textscan(fid, '%s', 1, 'delimiter', '\n', 'headerlines', ...
                  linum -1)

%getting the header:
tempHeader = header{1}{1};
headers = strsplit(tempHeader, '\t');
sampleCount = length(headers);

%getting the GOEID from a string. 
GEOID = cell(1, sampleCount);
expression = 'GSM[0123456789]+';
for j = 1:sampleCount
    GEOID(j) = regexp(headers{j}, expression, 'match');
end

%I need length of header for the sampleCount
fileName = ['/home/mfarahbod/data/affyArray/' ...
               'bloodCellTypeDataset/GSE22886_A.txt'];
fid = fopen(fileName)
data = textscan(fid, [repmat('%s', 1,1) repmat('%f', 1, ...
                                               sampleCount)], ...
                'Headerlines', linum, 'Delimiter', '\t');
fclose(fid);

gCount = length(gpl96.uniqueGeneSymbols)
pCount = length(data{1});
dataExpr = cell2mat(data(2:end));

for j = 1:length(data{1})
    if(strcmp(data{1}{j}(1), '"'))
        data{1}{j} = data{1}{j}(2:end);
    end
    if(strcmp(data{1}{j}(end), '"'))
        data{1}{j} = data{1}{j}(1:(end-1));
    end
end

[a, b] = ismember(data{1}, gpl96.probeID);
geneSymbols = gpl96.geneSymbol(b(a)); %gene ID 22k
myProbes = data{1}(a); %probe names - 22k
myDataExpr = dataExpr(a ,:); %probe values
pCount = length(myProbes); 

geneMap = containers.Map(gpl96.uniqueGeneSymbols, [1:gCount]);
probeMap = containers.Map(myProbes, geneSymbols);%  NCBI geneID

probeDataSet.mat = myDataExpr;
probeDataSet.GEOID = GEOID;

% % on data, but on the platform
divMat = zeros(gCount, sampleCount);

%    adding the maps
tic
for j = 1:pCount
    sib = values(geneMap, values(probeMap, myProbes(j)));
    divMat(sib{1}, :) = divMat(sib{1}, :) + 1;
end
toc

%getting the porbe list for each gene
%max probe count for each gene
probeInd = zeros(gCount, max(max(divMat))); 
for j = 1:pCount
    sib = values(geneMap, values(probeMap, myProbes(j)));
    temp = find(probeInd(sib{1}, :) == 0);
    pInd = temp(1);
    probeInd(sib{1}, pInd) = j;
end

probeMean = mean(probeDataSet.mat, 2);
a = mean(probeMean)
newGeneExpr = zeros(gCount, sampleCount);
affectedGeneCount = 0;

q = quantile(probeDataSet.mat(:), 5)
thr = q(1);

% for each gene
for j = 1: length(probeInd)
    pCount = sum(probeInd(j,:)>0);
    %if there are any with exp > thr, get the mean of bigger ones
    largeProbes = find(probeMean(probeInd(j,1:pCount))> thr);

    if(length(largeProbes) > 0)
        largeProbesInd = probeInd(j, largeProbes);
        newGeneExpr(j,:) = mean(probeDataSet.mat(largeProbesInd, ...
                                                   :), 1);
        
        if(pCount > length(largeProbes))
            affectedGeneCount = affectedGeneCount + 1;
        end
    else
        newGeneExpr(j, :) = mean(probeDataSet.mat(probeInd(j, ...
                                                          1:pCount), ...
                                                    :), 1);
    end
end

finalMat = newGeneExpr;
dataSet.mat = finalMat;
dataSet.GEOID = probeDataSet.GEOID;
dataSet.geneSyms = gpl96.uniqueGeneSymbols;

save([dataFolder 'GSE22886A_DS.mat'], 'dataSet');

% >>>>>>>>> for samples B
fileName = ['/home/mfarahbod/data/affyArray/' ...
               'bloodCellTypeDataset/GSE22886_B.txt'];
fid = fopen(fileName)
linum = 1
header = textscan(fid, '%s', 1, 'delimiter', '\n', 'headerlines', ...
                  linum -1)

%getting the header:
tempHeader = header{1}{1};
headers = strsplit(tempHeader, '\t');
sampleCount = length(headers);

%getting the GOEID from a string. 
GEOID = cell(1, sampleCount);
expression = 'GSM[0123456789]+';
for j = 1:sampleCount
    GEOID(j) = regexp(headers{j}, expression, 'match');
end

%I need length of header for the sampleCount
fileName = ['/home/mfarahbod/data/affyArray/' ...
               'bloodCellTypeDataset/GSE22886_B.txt'];
fid = fopen(fileName)
data = textscan(fid, [repmat('%s', 1,1) repmat('%f', 1, ...
                                               sampleCount)], ...
                'Headerlines', linum, 'Delimiter', '\t');
fclose(fid);

gCount = length(gpl97.uniqueGeneSymbols)
pCount = length(data{1});
dataExpr = cell2mat(data(2:end));

for j = 1:length(data{1})
    if(strcmp(data{1}{j}(1), '"'))
        data{1}{j} = data{1}{j}(2:end);
    end
    if(strcmp(data{1}{j}(end), '"'))
        data{1}{j} = data{1}{j}(1:(end-1));
    end
end

[a, b] = ismember(data{1}, gpl97.probeID);
geneSymbols = gpl97.geneSymbol(b(a)); %gene ID 22k
myProbes = data{1}(a); %probe names - 22k
myDataExpr = dataExpr(a ,:); %probe values
pCount = length(myProbes); 

geneMap = containers.Map(gpl97.uniqueGeneSymbols, [1:gCount]);
probeMap = containers.Map(myProbes, geneSymbols);%  NCBI geneID

probeDataSet.mat = myDataExpr;
probeDataSet.GEOID = GEOID;

% % on data, but on the platform
divMat = zeros(gCount, sampleCount);

%    adding the maps
tic
for j = 1:pCount
    sib = values(geneMap, values(probeMap, myProbes(j)));
    divMat(sib{1}, :) = divMat(sib{1}, :) + 1;
end
toc

%getting the porbe list for each gene
%max probe count for each gene
probeInd = zeros(gCount, max(max(divMat))); 
for j = 1:pCount
    sib = values(geneMap, values(probeMap, myProbes(j)));
    temp = find(probeInd(sib{1}, :) == 0);
    pInd = temp(1);
    probeInd(sib{1}, pInd) = j;
end

probeMean = mean(probeDataSet.mat, 2);
a = mean(probeMean)
newGeneExpr = zeros(gCount, sampleCount);
affectedGeneCount = 0;

q = quantile(probeDataSet.mat(:), 5)
thr = q(1);

% for each gene
for j = 1: size(probeInd, 1)
    pCount = sum(probeInd(j,:)>0);
    %if there are any with exp > thr, get the mean of bigger ones
    largeProbes = find(probeMean(probeInd(j,1:pCount))> thr);

    if(length(largeProbes) > 0)
        largeProbesInd = probeInd(j, largeProbes);
        newGeneExpr(j,:) = mean(probeDataSet.mat(largeProbesInd, ...
                                                   :), 1);
        
        if(pCount > length(largeProbes))
            affectedGeneCount = affectedGeneCount + 1;
        end
    else
        newGeneExpr(j, :) = mean(probeDataSet.mat(probeInd(j, ...
                                                          1:pCount), ...
                                                    :), 1);
    end
end

finalMat = newGeneExpr;
dataSet.mat = finalMat;
dataSet.GEOID = probeDataSet.GEOID;
dataSet.geneSyms = gpl97.uniqueGeneSymbols;

save([dataFolder 'GSE22886B_DS.mat'], 'dataSet');

dataFolder = '/home/mfarahbod/data/affyArray/bloodCellTypeDataset/'
% check out the heatmap with labels and sample correltaions - QC
load([dataFolder 'GSE22886A_DS.mat']);
Ads = dataSet;

load([dataFolder 'GSE22886B_DS.mat']);
Bds = dataSet;

% identification of samples and cell types - seems OK 
sib = corr(Ads.mat);

h = figure
heatmap(sib)

[a, b] = ismember(Ads.geneSyms, Bds.geneSyms);

h = figure
sib = corr(Bds.mat)
heatmap(sib)

halva = [Bds.mat; Ads.mat];
sib = corr(halva);
h = figure
heatmap(sib)

% - label samples

groups = {'B', 'CD14+', 'CD4+Te', 'CD4+T', 'CD8+T', 'IgGAmemB', ...
          'IgMmemB', 'Monocyte', 'NK', 'Neutrophils', 'pbm', 'ppb'}

labels = [repmat(5, 1, 4), repmat(4, 1, 14), repmat(3, 1, 6), ...
          repmat(9, 1, 15), repmat(1, 1, 7), repmat(6, 1, 4), ...
          repmat(7, 1, 4), repmat(12, 1, 3), repmat(11, 1, 4), ...
          repmat(8, 1, 36), repmat(2, 1, 12), repmat(10, 1, 5)]

% summary Ads 
samat = zeros(length(Ads.mat), 12);
for i = 1:12
    samat(:, i) = mean(Ads.mat(:, labels == i)');
end
Ads.samat = samat;

% summary Bds
sbmat = zeros(length(Bds.mat), 12);
for i = 1:12
    sbmat(:, i) = mean(Bds.mat(:, labels == i)');
end
Bds.samat = sbmat;

% >> obsolete
% extract the markers 
kado = Ads.mat(1, :);
h = figure
plot(labels, Ads.mat(1, :), 'o')
h = figure
plot(labels, Bds.mat(2, :), 'o')

sib = corr(Ads.mat);

h = figure
heatmap( sib, 'GridVisible', 'off')

h = figure
heatmap(labels)
colormap(jet)

% 2. getting list of DEG from the paper
%%%%%%%%%%%%%%%%%%%%%%%  read the table file for diff exp genes 

fileName = ['/home/mfarahbod/data/affyArray/bloodCellTypeDataset/' ...
            'SuppTable1_DEGs-Table1.csv']
fid = fopen(fileName)
diffExpMat = textscan(fid, ['%s' repmat('%f', 1, 22)], 'Delimiter', ',', ...
               'Headerlines',3)

diffGenes = diffExpMat{1};
diffMat = cell2mat(diffExpMat(2:end));

fileName = ['/home/mfarahbod/data/affyArray/bloodCellTypeDataset/' ...
            'SuppTable1_DEGs-Table1.csv']
fid = fopen(fileName)
header = fgetl(fid) % it is the third line

headers = strsplit(header, ',');
cellTypes = headers(2:23)

% count of marker for each cellType (they do overlap)
mCounts = sum(diffMat);
overlap = diffMat' * diffMat;
h = figure('units', 'centimeters', 'position', [0,0, 25,20])
heatmap(cellTypes, cellTypes, overlap)
colormap(1-gray)
title(['Count and overlap of marker genes idendified in different ' ...
       'cell types'])
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('heatmap_bloodMarkerGenes_countAndOverlap', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 3. get the bulk tissue data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% code from (GTEx_cellTypeCorrection.m)
load(['~/data/GTEx/' ...
      'GTExAllDataSetFromGeneLevel_v6p_newWithBlood_RPM.mat'])

uniqueGenes = unique(rpmGTExDS.genes);
myConMap = containers.Map(rpmGTExDS.genes, 1: ...
                          length(rpmGTExDS.genes));

GTExBloodDS = rpmGTExDS.dataSets(54);
GTExGenes = rpmGTExDS.genes;

myMat = zeros(length(uniqueGenes), size(GTExBloodDS.mat, 2));

for j = 1:length(myConMap)
    transInds = myConMap(uniqueGenes{j});
    tinyExpMat = GTExBloodDS.mat(transInds, :);
    if length(transInds) > 1 
        myMat(j, :) = sum(tinyExpMat);
    else 
        myMat(j, :) = tinyExpMat;
    end
end

bloodDS.genes = uniqueGenes;
bloodDS.mat = myMat;

save('~/resultsAndFigures/secondProject/bloodDS.mat', 'bloodDS')

[a, b] = ismember(diffGenes, uniqueGenes);

% loading blood network (code in buildingNetworks_GTEx_SC.m)
load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat')
bloodNet = GTExFiveNets.nets(1);
bloodExpGenes = bloodDS.genes(bloodNet.expGenes);

% 4. identify the clusters (code from clusteringGTExNetworks.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 1
gtexGeneSyms = GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(t).expGenes);
gtexNet = GTExFiveNets.nets(t).net005;
gtexND = sum(gtexNet + gtexNet');

gtexSNet = gtexNet(gtexND > 2, gtexND > 2);
gtexsfNet = gtexSNet + gtexSNet';
gtexFSyms = gtexGeneSyms(gtexND >2);

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

% save the two clusters
gtexCluster.cs1547 = t;
gtexCluster.cs1546 = t;
gtexCluster.syms = gtexFSyms;
save('~/resultsAndFigures/secondProject/gtexClusters_blood.mat', ...
     'gtexCluster')
load('~/resultsAndFigures/secondProject/gtexClusters_blood.mat')
% general cluster Info
cs = gtexCluster.cs1547;
clusterMemberCount = hist(cs, unique(cs))
inClusterIds = find(clusterMemberCount >= 20);
inClusterIds2 = find(clusterMemberCount >= 3);

inClusterIds = [find(inClusterIds), 229 ]

% 5. expression of the markers in different clusters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.1 for each of the clusters, how many marker genes we have, what are
% they. 

inDiffGenes = diffGenes(sum(diffMat') == 1);

for i = 1: 9%length(inClusterIds)
    i
    % get the genes
    cid = inClusterIds(outperm(i));
    
    % cgenes = gtexCluster.syms(ismember(gtexCluster.cs1547, [14, 26, ...
    %                216, 22, 228, 63]));
    
    % cgenes = gtexCluster.syms(ismember(gtexCluster.cs1547, [200, ...
    %                     215, 32, 28, 208, 210]));
    
    % cgenes = gtexCluster.syms(ismember(gtexCluster.cs1547, [232]));
    
    
    cgenes = gtexCluster.syms(gtexCluster.cs1547 == cid);
    cgc = length(cgenes)
    
    % see if there are more than 3 marker genes in them
    [a, b] = ismember(diffGenes, cgenes);
    
    if (sum(a) >= 3)
        myBar = sum(diffMat(a, :));
        h = figure
        bar(myBar)        
        set(gca, 'XTick', [1:22], 'XTickLabel', cellTypes)
        xtickangle(90)
        title(sprintf('clusterID %d, geneCount %d', cid, cgc))
        
        thisGenes = diffGenes(a);
        [aa, bb] = ismember(thisGenes, Ads.geneSyms);
        plotMat = zscore(Ads.samat(bb(aa), :)');
        h = figure
        heatmap(groups, 1:sum(aa), plotMat');
        title(sprintf('clusterID %d, geneCount %d', cid, cgc))
    end
    
    % plot the histogram of celltypes for that cluster with cluster
    % ID on top
end

% >>>>>>> 5.2 where are the markers: which clusters % getting it from
% diffMat
fileName = ['/home/mfarahbod/data/affyArray/bloodCellTypeDataset/' ...
            'SuppTable1_DEGs-Table1.csv']
fid = fopen(fileName)
diffExpMat = textscan(fid, ['%s' repmat('%f', 1, 22)], 'Delimiter', ',', ...
               'Headerlines',3)

diffGenes = diffExpMat{1};
diffMat = cell2mat(diffExpMat(2:end));

fileName = ['/home/mfarahbod/data/affyArray/bloodCellTypeDataset/' ...
            'SuppTable1_DEGs-Table1.csv']
fid = fopen(fileName)
header = fgetl(fid) % it is the third line

headers = strsplit(header, ',');
cellTypes = headers(2:23)

diff.genes = diffGenes;
diff.mat = diffMat;
diff.cellTypes = cellTypes;

save('~/resultsAndFigures/secondProject/bloodDiffGenesInfo.mat', 'diff') % gtexCluster

for i = 1:size(diffMat, 2)
    i

    % thisDiffGenes = diffGenes(logical((diffMat(:, 13) - diffMat(:, ...
    %                                                   22))>0));
    %thisDiffGenes = inM(13).list;
    % thisDiffGenes = markerListsReO(243:272);
    % for i = 1:9
    % thisDiffGenes = homoMarkerList{i};
    % for i = 1:length(IDs)
    %     cid = IDs(i)
    %      thisDiffGenes = myNetGenes(c == cid);
    
    % >>> first, getting it from that thing. 
    thisDiffGenes = diffGenes(logical(diffMat(:, i)));
    [a, b] = ismember(thisDiffGenes, gtexCluster.syms);
    sum(a)
    
    % get the genes
    cids = gtexCluster.cs1547(b(a));
    uniqueCids = unique(cids);
    myBar = hist(cids, uniqueCids);
    tpcs = find(myBar > 2); % three plus clusters
    if length(tpcs > 0)
        h = figure
        bar(myBar(tpcs))        
        set(gca, 'XTick', [1:length(uniqueCids(tpcs))], 'XTickLabel', ...
                 uniqueCids(tpcs))
        title(sprintf('%s - total: %d', cellTypes{i}, sum(a)))
        figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/MarkersInClusters/']
        file = sprintf('%sDiffGenes_%d', figFolder, i);
        set(h, 'PaperOrientation', 'landscape')
        %print(h, '-dpng', [file '.png'])
        print(h, '-deps', [file '.eps'])
        print(h, '-dpdf', [file '.pdf'])
        saveas(h, [file '.eps'], 'epsc')
    end
end

% >>>> 5.3 where are the markers: which clusters % getting it from
% inM

for i = 1:length(inM)
    i
    thisDiffGenes = inM(i).list;
    [a, b] = ismember(thisDiffGenes, gtexCluster.syms);
    sum(a)
    
    % get the genes
    cids = gtexCluster.cs1547(b(a));
    uniqueCids = unique(cids);
    myBar = hist(cids, uniqueCids);
    tpcs = find(myBar > 2); % three plus clusters
    
    if (length(uniqueCids) == 1)
        myBar = sum(a);
        tpcs = 1;
    end

    if length(tpcs > 0)
        h = figure
        bar(myBar(tpcs))        
        set(gca, 'XTick', [1:length(uniqueCids(tpcs))], 'XTickLabel', ...
                 uniqueCids(tpcs))
        title(sprintf('%d - total: %d', i, sum(a)))
        figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/MarkersInClusters/']
        file = sprintf('%sinM_filteredMarkers_%d', figFolder, i);
        set(h, 'PaperOrientation', 'landscape')
        %print(h, '-dpng', [file '.png'])
        print(h, '-deps', [file '.eps'])
        print(h, '-dpdf', [file '.pdf'])
        saveas(h, [file '.eps'], 'epsc')
    end
end

% 6. regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bloodDS
bloodNet
whos bloodExpGenes

% get bloodDS 
[a, b] = ismember(bloodExpGenes, bloodDS.genes);
bloodExpMat = log2(bloodDS.mat(b(a), :)+1);

% >> getting all diff genes: 
[a, b] = ismember(diffGenes, bloodExpGenes);
GTExDiffExpMat = bloodExpMat(b(a), :);

% >> 6.1 getting only the marker genes - getting the R2s
load(['~/resultsAndFigures/secondProject/' ...
      'bloodCorrelatedMarkers.mat'])
whos inM
% count of markers in inM
 
mCount = length(inM(1).list);
mCounts = zeros(1,15);
mCounts(1) = mCount;
markerLists = inM(1).list;
for i = 2:13
    mCount = mCount + length(inM(i).list);
    mCounts(i) = length(inM(i).list);
    markerLists = [markerLists; inM(i).list];
end

markerLists = [markerLists; Bmarkers'; ERmarkers'];
markerLists = [markerLists; m177List'; m163List']; % this is what I
                                                   % am plotting
mCounts(14) = length(m177List);
mCounts(15) = length(m163List);
inMlistOrder = [1 2 3 4 12 6 7 9 10 5 8 11 13]
 
markerListsReO = inM(1).list;
for i = 2:13
    order = inMlistOrder(i);
    markerListsReO = [markerListsReO; inM(order).list];
end
markerListsReO = [markerListsReO; m177List'; m163List']; % this is what I

[a, b] = ismember(markerListsReO, bloodExpGenes);
markerExpMat = bloodExpMat(b(a), :);
sib = corr(markerExpMat');

h = figure
heatmap(sib, 'GridVisible', 'off')
colormap(parula)

% OR get the markers as this and that 
load('~/resultsAndFigures/secondProject/bloodDiffGenesInfo.mat') 
diff.cellTypes{:}
% cellTypes we are getting: 
tempSelected = [1, 2, 4:13, 22];

inGenes = diff.genes((sum(diff.mat(:, tempSelected)') > 0));
inGenes = diff.genes;
length(inGenes)

[a, b]  = ismember(inGenes, bloodExpGenes);
markerExpMat = bloodExpMat(b(a), :);

% getting the prediction
[w, score, latent, tsquared, explained, mu] = pca(markerExpMat', ...
                                                  'NumComponents', 10);
T = GTExDiffExpMat' * w;

% correct it for all the genes
featureMat = score(:, 1:7);
regOut = zeros(size(bloodExpMat));
coeffs = zeros(size(bloodExpMat, 1), 8);
r2 = zeros(size(bloodExpMat , 1), 1);
for i = 1:size(bloodExpMat, 1)
    i
    y = (bloodExpMat(i, :));
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

%result.markers = markerLists;
result.markers = inGenes;
result.coeff = coeffs;
result.featureMat = featureMat;
result.regOut = regOut;
result.r2 = r2;
save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_diffGenes_all.mat'],'result', ...
     '-v7.3')

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_diffGenes_selectedTypes.mat'],'result', ...
     '-v7.3')

[a, b] = ismember(gtexCluster.syms, bloodExpGenes);

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240_Eryth177_Eryth163_5feature.mat'],'result', ...
     '-v7.3')

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240_Eryth177_Eryth163.mat'])
regResM = result;
clear result

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240.mat'])
regRes = result;

h = figure
hist(r2)

% >> 6.2 getting the random MgeneR2s

meanExps = (mean(bloodExpMat'));
markerMeans = (mean(markerExpMat'));

dataSetExp = bloodExpMat;
[allMarkers, b] = ismember(bloodExpGenes, markerLists);

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

    mockMarkerExp = bloodExpMat(mockMarkerInds, :);
    [w, score, latent, tsquared, explained, mu] = pca(mockMarkerExp', ...
                                                      'NumComponents', ...
                                                      10);
    explained(1:10)

    % correct it for all the genes
    featureMat = score(:, 1:5);
    regOut = zeros(size(bloodExpMat));
    coeffs = zeros(size(bloodExpMat, 1), 6);
    r2 = zeros(size(bloodExpMat , 1), 1);
    for i = 1:size(bloodExpMat, 1)
        j
        i
        y = bloodExpMat(i, :);
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

    resultR(j).coeff = coeffs;
    resultR(j).featureMat = featureMat;
    %    result(j).regOut = regOut;
    resultR(j).r2 = r2;
    resultR(j).mmg = mockMarkerInds;
end

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionRegression_bloodData_240RandomGenes_100trial' ...
      '.mat'],'resultR', '-v7.3')

% get the featuer vector for two 177 and 163 - 

% [a, b] = ismember(gtexCluster.syms, bloodExpGenes);
% cr2 = r2(b(a));
% cmeanExp = mean(bloodExpMat(b(a), :)');
% cbvar = var(bloodExpMat(b(a), :)');
% [a, b] = ismember(gtexCluster.syms, Ads.geneSyms);

% ctp = Ads.samat(b(a), :);
% ctpVar = var(ctp');

[a, b] = ismember(bloodExpGenes, Ads.geneSyms);
[a, b] = ismember(bloodExpGenes, dataSet.hhomoGeneSyms);
bvar = var(bloodExpMat(a, :)');
meanExp = mean(bloodExpMat(a, :)');

corr(bvar', result.r2(a))
corr(bvar', ctpVar')

corr(meanExp', bvar')
corr(meanExp', ctpVar')

corr(ra, ctpVar')

meanSelect = meanExp > 5;
corr(meanExp(meanSelect)', bvar(meanSelect)')
corr(bvar(meanSelect)', ctpVar(meanSelect)')
ra = result.r2(a);
halva = [meanExp(meanSelect)', bvar(meanSelect)', ctpVar(meanSelect)', ...
         ra(meanSelect)];
corr(halva)
h = figure
heatmap({'meanExp', 'bulkVar', 'ctpVar', 'R2'}, {'meanExp', 'bulkVar', 'ctpVar', 'R2'},corr(halva) - eye(4))

ctp = Ads.samat(b(a), :);
%ctp = Ads.mat(b(a), :);
ctpVar = var(ctp');

ctpVar = var(dataSet.hMat(b(a), :)');

% >>> 6.4 get the average feature from R2 for the clusters, are
% they orthogonal? At all?

pMat = zeros(length(inClusterIds), 7);
geneInCount = zeros(1, length(inClusterIds));
for i = 1:length(inClusterIds)
    % get genes in the cluster
    i
    geneSet = gtexCluster.syms(gtexCluster.cs1547 == inClusterIds(i));
              
    %>>> doing it with just Ads
    % [a, b] = ismember(geneSet, Ads.geneSyms);
    % % get them in the expression data
    % smallExp = (Ads.samat(b(a), :));
    % whos smallExp
    % %inGenes = sum(smallExp' > 8) > 3;
    
    % % the average 
    % %avgTP(i, :) = mean(smallExp(inGenes,:));
    % avgCTP(i, :) = mean(zscore(smallExp')');
    
    %>>> doing it with just Ads and Bds
    [a, b] = ismember(geneSet, bloodExpGenes);
    pMat(i, :) = mean(regRes.coeff(b, 2:end));
end
h= figure
kado = corr(pMat(outperm, :)');
myGTExFMat = pMat(outperm, :);
heatmap(clusterLabels(outperm), clusterLabels(outperm), kado)
colormap(1- [1-gray' copper']')
% get R2s for marker
[a, b] = ismember(markerLists, bloodExpGenes);
markerR2s = regRes.r2(b);
qsa = quantile(markerR2s, 30);

[a, b] = ismember(gtexCluster.syms, bloodExpGenes);
clusterGeneR2s = regRes.r2(b);
qsb = quantile(clusterGeneR2s, 30)
h = figure
plot(qsa, qsb, '.')
xlabel('markers')
ylabel('cluster')

xlim([.5 1])
ylim([.5 1])
hold on
plot([.5 1], [.5 1])


% 7. Result 2: study of the R2 values for different clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the genes
meanClusterRs = zeros(1, 238);
for i = 1:length(meanClusterRs)
    %cid = inClusterIds(i);
    cid = i
    cgenes = gtexCluster.syms(gtexCluster.cs1547 == cid);
    cgc = length(cgenes)
    
    [a, b] = ismember(bloodExpGenes, cgenes);
    %    myrs = result.r2(a);
    myrs = regResM.r2(a);
    meanClusterRs(i) = mean(myrs);
end

inClusterIds = find(inClusterIds);
incr2 = meanClusterRs(inClusterIds)
incmc = clusterMemberCount(inClusterIds)

temp = find(incr2 < .5)
inClusterIds(temp)
gtexCluster.syms(gtexCluster.cs1547 == 175)

h = figure
plot(sort(incr2), 'o')

% 8. cluster functional enrichment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the functional enrichment for all the clusters 
% the code is fenrich_function.m

load('~/data/general/GOdata_GTExV6_nonComp.mat')
GTExGO = GOdata;

[a, b] = ismember(bloodExpGenes, GTExGO.geneSymbols);
wholeFMat = GTExGO.matP(b(a), :);
sib = sum(wholeFMat);
inTerms = ((sib >=3) + (sib<=200)) == 2;

fMat = wholeFMat(:, inTerms);
inTermNames = GTExGO.GOTerms(inTerms);
inTermsGOID = GTExGO.GOID(inTerms);
rawGeneCount = sib(inTerms);
expGeneCounts = sum(fMat);
fdrThr = .1
for i = 1:length(inClusterIds)
    fenrichList = gtexCluster.syms(gtexCluster.cs1547 == inClusterIds(i));
    [a, b] = ismember(fenrichList, bloodExpGenes);
    cfMat = sum(fMat(b(a), :));
    inFs = cfMat >= 3;
    if sum(inFs) > 1
        i
        sum(inFs)
        fMatGeneS = fMat(b(a), inFs);

        myFunctionMat = fMat(b(a), inFs);

        % adjust the last parameter accordingly
        % ps = 1 - hygecdf(cfMat(inFs) -1, length(dataSet.genes), expGeneCounts(inFs), ...
        %                  sum(a));
        ps = 1 - hygecdf(cfMat(inFs) -1, length(bloodExpGenes), expGeneCounts(inFs), ...
                         sum(a));


        fdrs = ps * sum(inFs);
        myNames = inTermNames(inFs);
        passedFDRTerms = myNames(fdrs < fdrThr);
        myIDs = inTermsGOID(inFs);
        passedFDRIDs = myIDs(fdrs < fdrThr);
        enrichMat = fMatGeneS(:, fdrs<fdrThr);
        passedFDRs = fdrs(:, fdrs<fdrThr);
        clusterFE(i).myNames = inTermNames(inFs);
        clusterFE(i).passedFDRTerms = myNames(fdrs < fdrThr);
        clusterFE(i).myIDs = inTermsGOID(inFs);
        clusterFE(i).passedFDRIDs = myIDs(fdrs < fdrThr);
        clusterFE(i).enrichMat = fMatGeneS(:, fdrs<fdrThr);
        clusterFE(i).passedFDRs = fdrs(:, fdrs<fdrThr);
        clusterFE(i).enrichMat = fMatGeneS(:, fdrs<fdrThr);
        clusterFE(i).clusterID = inClusterIds(i);
        clusterFE(i).geneList = fenrichList;
    end
end

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodCluster_fenrich.mat'],'clusterFE', ...
     '-v7.3')
fe = load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodCluster_fenrich.mat'])

% 8.2 >>> which brain clusters are they
brainc = load(['~/resultsAndFigures/secondProject/' ...
               'gtexClusters.mat']);
i = i
clusterFE(outperm(i)).passedFDRTerms(:)
i = i + 1

for i = 1:length(outperm)
    myID = inClusterIds(outperm(i))
    geneList = gtexCluster.syms(gtexCluster.cs1547 == myID);

    [a, b] = ismember(brainc.gtexCluster.syms, geneList);
    inds = brainc.gtexCluster.cs1547(a);
    myBar = hist(inds, unique(inds));
    h = figure
    bar(myBar)        
    set(gca, 'XTick', [1:length(unique(inds))], 'XTickLabel', ...
             unique(inds))
    title(sprintf('%d - %d', myID, sum(a)))
    sum(a)
    myID
    figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/bloodBrainClusterDist/']
    file = sprintf('%sblood_ID%d', figFolder, myID);
    set(h, 'PaperOrientation', 'landscape')
    %print(h, '-dpng', [file '.png'])
    print(h, '-deps', [file '.eps'])
    print(h, '-dpdf', [file '.pdf'])
    saveas(h, [file '.eps'], 'epsc')
end

% 8.3 >>> identification of blood specific functional terms
load(['~/resultsAndFigures/secondProject/moduleRepresentation/' ...
      'funcRes_withAffyBloodLiver_FDR05.mat'])

sib = funcRes.funPresence * funcRes.funPresence';
bs = sum(funcRes.funPresence([5,8], :)) > 0 ;

neg = sum(funcRes.funPresence([1:4], :)) >0; % removing brain, but
                                             % not liver as much
                                             % liver is inflamation

%passing the terms in brain clusters 52, 65, 105, 164, 190, 233
load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'functionalEnrichment.mat'])

[a, b] = ismember(clusterFRes.clusterIDs, [52, 65, 105, 164, 190, ...
                    233]);
myMat = clusterFRes.ps * length(clusterFRes.fTerms);
passFDR = ((myMat > 0) + (myMat < .1)) == 2;

myFs = passFDR(a, :);
sib = sum(myFs);
myTermIDs = clusterFRes.fIDs(sib > 0); % fterms that belong to
                                       % microglia clusters
% I should also add to this, cluster IDs 
[a, b] = ismember(myTermIDs, funcRes.uniqueTermIDs);
neg(b(a)) = false;

% get the blood specific terms
kado = (bs - neg) > 0;
sum(kado)

bloodSterms = funcRes.uniqueTermNames(kado);
bloodSterms(:)

% 9. gene expression levels in different clusters versus celltypes
%%%%%%%%%%%%%% plotting the expression level of diff genes in
%%%%%%%%%%%%%% clusters

IDs = inClusterIds;
IDs = affyInClusters;
for i = 1:length(IDs)
    %for i = 1:length(inMB)    
    i
     cid = IDs(i);
     % cgenes = myNetGenes(c == cid);
     cgenes = gtexCluster.syms(gtexCluster.cs1547 == cid);
    %cgenes = affyCluster.syms(affyCluster.cs1547 == cid);
    
    % cid = i
    %         cgenes = inMB(1).list;
    % cgenes = inM(11).list;
    
    cgc = length(cgenes)
    
    % >>> this is for the affy human cell types
    [a, b] = ismember(cgenes, [Ads.geneSyms', Bds.geneSyms']);
    sum(a)
    %wholeMat = [Ads.mat' Bds.mat']';
    wholeASMat = [Ads.samat' Bds.samat']';
    mybox = wholeSMat(b(a), :);
    h = figure
        subplot(1, 2, 1)
    boxplot(zscore(mybox')')
       set(gca, 'XTick', [1:12], 'XTickLabels', groups)
    xtickangle(90)
    title(sprintf('clusterId: cs1547-%d, geneCount: %d', cid, ...
                      cgc))
    subplot(1, 2, 2)        
    boxplot(mybox)
       set(gca, 'XTick', [1:12], 'XTickLabels', groups)
    xtickangle(90)
    
figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/clusterCTProfile/']
file = sprintf('%sclusterCTP_%d', figFolder, cid);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


    % >> this is for mouse cell types
    [a, b] = ismember(cgenes, ctDataSet.dataSet.hhomoGeneSyms);
    sum(a)
    mybox = ctDataSet.dataSet.hMat(b(a), :);
    h = figure
    subplot(1, 2, 1)
    boxplot(zscore(mybox')')
        title(sprintf('clusterId: cs1547-%d, geneCount: %d', cid, ...
                      cgc))

        %       set(gca, 'XTick', [1:12], 'XTickLabels', groups)
    % xtickangle(90)
    % title(sprintf('clusterId: cs1547-%d, geneCount: %d', cid, ...
    %                   cgc))
    subplot(1, 2, 2)        
    boxplot(mybox)
    %    set(gca, 'XTick', [1:12], 'XTickLabels', groups)
    % xtickangle(90)
    
    
figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/clusterCTProfile/']
file = sprintf('%sclusterCTP_%d_mm', figFolder, cid);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

end

cid = 180;
cgenes = gtexCluster.syms(gtexCluster.cs1546 == cid);
cgenes = gtexCluster.syms(gtexCluster.cs1547 == cid);

cgenes = gtexCluster.syms(ismember(gtexCluster.cs1547, [93, 229, 44, ...
                   51, 27, 108, 204, 77]));

cgenes = gtexCluster.syms(ismember(gtexCluster.cs1547, [43, ...
                   12, 21, 38]));

    cgenes = gtexCluster.syms(ismember(gtexCluster.cs1547, [14, 26, ...
                   216, 22, 228, 63]));
    
    cgenes = gtexCluster.syms(ismember(gtexCluster.cs1547, [200, ...
                        215, 32, 28, 208, 210]));
    
    cgenes = inMB(1).list;
    cgenes= markerLists(81:89);

cid = 181;
cgenes = gtexCluster.syms(gtexCluster.cs1547 == cid);

cgc = length(cgenes)
[a, b] = ismember(cgenes, [Ads.geneSyms', Bds.geneSyms']);
sum(a)
wholeMat = [Ads.mat' Bds.mat']';
wholeSMat = [Ads.samat' Bds.samat']';
mybox = wholeSMat(b(a), :);
mybox = wholeMat(b(a), :);
mybox = Ads.samat((b(a)), :);
mybox = Ads.mat((b(a)), :);
mg = mean(mybox');
mcmat = zeros(size(mybox));
for i = 1:length(mybox)
    mcmat(i, :) = mybox(i, :) - mg(i);
end

h = figure
hold on
for i = 1:9
    plot(zscore(mybox(i, :)))
end

h = figure
boxplot(mybox)
plot(median(mybox))
hold on
boxplot(mcmat)
boxplot(zscore(mybox')')
set(gca, 'XTick', [1:12], 'XTickLabels', groups)
xtickangle(90)

set(gca, 'XTick', [1:114], 'XTickLabels', Ads.GEOID)
xtickangle(90)

% I want genes which are expressed at high levels at least in some
% samples. 
book = mybox > 8;
sib = sum(book') > 3;
h = figure
heatmap(mybox(sib, :))

med1 = median(mybox);
med2 = median(mybox);
med3 = median(mybox);
med4 = median(mybox);
h = figure
plot(med1)
hold on
plot(med2)
hold on
plot(med3)
hold on
plot(med4)

legend('m1', 'm2', 'm3')

whos bloodExpGenes

% get bloodDS 
[a, b] = ismember(cgenes, bloodExpGenes);
smallMat = bloodExpMat(b(a), :);

h = figure
heatmap(smallMat, 'GridVisible', 'off')
mean(smallMat')

h = figure
hist(mean(bloodExpMat'))

cgenes = markerLists(62:70);
    
        [a, b] = ismember(cgenes, dataSet.hhomoGeneSyms);
    sum(a)
    mybox = dataSet.hMat(b(a), :);
    h = figure
    subplot(1, 2, 1)
    boxplot(zscore(mybox')')
        title(sprintf('clusterId: cs1547-%d, geneCount: %d', cid, ...
                      cgc))

        %       set(gca, 'XTick', [1:12], 'XTickLabels', groups)
    % xtickangle(90)
    % title(sprintf(clusterId: cs1547-%d, geneCount: %d', cid, ...
    %                   cgc))
    subplot(1, 2, 2)        
    boxplot(mybox)

% 10. print clusters for cytoscape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 10.1 getting the density of clusters
gtexCluster
bloodNet
whos bloodExpGenes
gtexFullNet = bloodNet.net005 + bloodNet.net005';

ccount = length(inClusterIds);
overlapDensity = zeros(ccount, ccount);
for i = 1:ccount
    i
    geneSet1 = gtexCluster.syms(gtexCluster.cs1547 == inClusterIds(i)); 
    [a, b1] = ismember(geneSet1, bloodExpGenes);
    overlapDensity(i, i) = sum(sum(gtexFullNet(b1, b1)))/(length(a) *(length(a)-1));
    for j = i+1:ccount
        geneSet2 = gtexCluster.syms(gtexCluster.cs1547 == inClusterIds(j)); 
        [a, b2] = ismember(geneSet2, bloodExpGenes);
        chunk = gtexFullNet(b1, b2);
        overlapDensity(i, j) = max(sum(sum(chunk))/(length(b1) * ...
                                                (length(b2))), 0.001);

    end
end

overlapDensity = full(overlapDensity);
fullod = (overlapDensity + triu(overlapDensity, ...
                                                  1)');

h = figure
heatmap(fullod(outperm, outperm), 'GridVisible', 'off')

h = figure('units', 'centimeters', 'position', [0,0, 25, 25])
heatmap(clusterLabels(outperm), clusterLabels(outperm), fullod(outperm, ...
                                                  outperm), ...
        'GridVisible', 'off')

heatmap(cinfo.clusterLabels(selectedCs), cinfo.clusterLabels(selectedCs), fullod(selectedCs, ...
                                                  selectedCs), ...
        'GridVisible', 'off')

colormap(1-copper)
figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%sGTExBloodNetoverlap', figFolder);
file = sprintf('%sGTExBloodNetoverlap_selectedCs', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


cc = length(unique(gtexCluster.cs1547))
firstCs = [27, 44, 51, 93, 229];
odtest = zeros(5, cc)
for i = 1:5
    i
    geneSet1 = gtexCluster.syms(gtexCluster.cs1547 == firstCs(i)); 
    [a, b1] = ismember(geneSet1, bloodExpGenes);
    for j = 1:cc
        geneSet2 = gtexCluster.syms(gtexCluster.cs1547 == j); 
        [a, b2] = ismember(geneSet2, bloodExpGenes);
        chunk = gtexFullNet(b1, b2);
        odtest(i, j) = sum(sum(chunk))/(length(b1) *(length(b2)));
    end
end


fileName = 'clusterInfo_blood_full.txt'
file = ['~/resultsAndFigures/secondProject/cytoFilesForClusters/' fileName]
fid = fopen(file, 'w')
fprintf(fid, ['clusterID\t geneCount\t density \t R2\n'])
for i = 1:length(inClusterIds)
    ID = inClusterIds(i);
    gCount = incmc(i);
    density = overlapDensity(i, i);
    R2 = incr2(i);
    fprintf(fid, '%d\t%d\t%.2f\t%.2f\n', ID, gCount, density, R2);
end
fclose(fid)


fileName = 'clusterInfo_blood_full_R2From240M_ERplus.txt'
file = ['~/resultsAndFigures/secondProject/cytoFilesForClusters/' fileName]
fid = fopen(file, 'w')
fprintf(fid, ['clusterID \t R2_ERp\n'])
for i = 1:length(inClusterIds)
    ID = inClusterIds(i);
    %gCount = incmc(i);
    %density = overlapDensity(i, i);
    R2 = incr2(i);
    fprintf(fid, '%d\t%.2f\n', ID, R2);
end
fclose(fid)

% >>>>> clusternetwork file
% for each two clusters, I need to have their link density - that's
% all.
fileName = 'interClusterInfo_blood_full.txt'
file = ['~/resultsAndFigures/secondProject/cytoFilesForClusters/' fileName]
fid = fopen(file, 'w')
[a, b, c] = find(triu(overlapDensity, 1));
fprintf(fid, ['clusterID_1\t clusterID_2\t density\n'])
for i = 1:length(c)
    i
    ID1 = inClusterIds(a(i));
    ID2 = inClusterIds(b(i));
    d = c(i);
    fprintf(fid, '%d\t %d\t %.2f\n', ID1, ID2, d);
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

% 11. cluster the clusters and sort them out and finish
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure
heatmap(inClusterIds(outperm), inClusterIds(outperm), fullod(outperm, ...
                                                  outperm))
colormap(1- gray)
IDlist = [228, 22, 216, 26, 14, 226, 48, 78, 236, 68, 90, 235, 230, ...
         217, 209, 131,...
         63, 200, 215, 32, 218, 28, 59, 122, 210, 96, 227, 75, 208, ...
         212, 221, 232, 113, 220, 61, 116, 38, 83, 21, 55, 211, ...
         153, 124, 149, 120, 157, 167, 177, 163, 142, 130, 151, 238, ...
         175, 237, 141]

[a, b] = ismember(IDlist, inClusterIds);

h = figure
heatmap(IDlist, IDlist, fullod(b(a), b(a)))

clusterDistMat = squareform(tril(1-fullod, -1));
 
z = linkage(clusterDistMat, 'average');
h = figure
[H, T, outperm2] = dendrogram(z, 140);
[H, T, outperm] = dendrogram(z, 63);

h = figure
heatmap(inClusterIds(outperm), inClusterIds(outperm) , fullod(outperm, ...
                                                  outperm), ...
        'GridVisible', 'off')
colormap([cool' spring']')

h = figure
heatmap(inClusterIds2(outperm2), inClusterIds2(outperm2), fullod(outperm2, outperm2), 'GridVisible', 'off')

% 11.5 get the coorelation of CT profiles and expression in
% celltypes - both mouse and human
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

avgCTP = zeros(length(inClusterIds), 12);
medCTP = zeros(length(inClusterIds), 12);
avgCTPwhole = zeros(length(inClusterIds), 114);
geneInCount = zeros(1, length(inClusterIds));
avgCTPmouse = zeros(length(inClusterIds), 20);
geneInCountMouse = zeros(1, length(inClusterIds));
for i = 1:length(inClusterIds)
    % get genes in the cluster
    i
    geneSet = gtexCluster.syms(gtexCluster.cs1547 == inClusterIds(i));
              
    %>>> doing it with just Ads
    % [a, b] = ismember(geneSet, Ads.geneSyms);
    % % get them in the expression data
    % smallExp = (Ads.samat(b(a), :));
    % whos smallExp
    % %inGenes = sum(smallExp' > 8) > 3;
    
    % % the average 
    % %avgCTP(i, :) = mean(smallExp(inGenes,:));
    % avgCTP(i, :) = mean(zscore(smallExp')');
    
    %>>> doing it with just Ads and Bds
    [a, b] = ismember(geneSet, [Ads.geneSyms', Bds.geneSyms']);
    geneInCount(i) = sum(a);
    wholeSMat = [Ads.samat' Bds.samat']';
    smallExp = wholeSMat(b(a), :);
    avgCTP(i, :) = mean(zscore(smallExp')');
    medCTP(i, :) = median(smallExp);
    
    % % %>>> doing it with mouse
    % [a, b] = ismember(geneSet, dataSet.hhomoGeneSyms);
    % geneInCountMouse(i) = sum(a);
    % smallExp = dataSet.hMat(b(a), :);
    % avgCTPmouse(i, :) = mean(zscore(smallExp')');
    
    %>>> doing it with just Ads and Bds
    [a, b] = ismember(geneSet, [Ads.geneSyms', Bds.geneSyms']);
    geneInCount(i) = sum(a);
    wholeMat = [Ads.mat' Bds.mat']';
    smallExp = wholeMat(b(a), :);
    avgCTPwhole(i, :) = mean(zscore(smallExp')');
end

grouporder = [1 6 7, 2 8, 3 4 5 9, 10 11, 12 ]
grouporder = [1 6 7, 2 8, 3 4 5 9, 10]
grouporder = [8, 10, 1 3 4 5 9 32]

plotMat = [avgCTP, avgCTPmouse];
h = figure
heatmap(groups(grouporder), [clusterLabels(outperm)], zscore(plotMat(outperm, grouporder ...
                                                  )')');

heatmap(zscore(plotMat(outperm, grouporder)')')
colormap(1- [1-gray' copper']')
h = figure
heatmap([1:20],clusterLabels(outperm), avgCTPmouse(outperm, :))
heatmap([1:114], inClusterIds(outperm), zscore(avgCTP(outperm, : ...
                                                  )')');
colormap(1- [1-gray' copper']')
title('mean per cluster on exp zscore')
title('mean per cluster on exp zscore - exp filtered')

sib = corr((avgCTP(outperm, grouporder)'));
h = figure
heatmap(inClusterIds(outperm), inClusterIds(outperm), sib)
colormap(1- [1-gray' copper']')
title('correlation of mean per cluster on exp zscore')

h = figure
heatmap([1:114], inClusterIds(outperm), avgCTP(outperm,:))
colormap(1- [1-gray' copper']')
title('exp  not filtered')
 
clusterOrder = outperm;
clusterOrder(62) = [];

sib = corr(avgCTP(clusterOrder, :)', 'type', 'Pearson');

h = figure
heatmap(inClusterIds(clusterOrder(1:54)), inClusterIds(clusterOrder(1:54)), ...
        sib(1:54, 1:54))
colormap([cool' spring']')


sib = corr(Ads.mat');
halva = sib(:);
h = figure
hist(halva, 40)

% 12. working with diff genes and their expression AND filtering
% markers based on coexpression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
whos diffGenes
whos diffMat
whos bloodExpMat
whos bloodExpGenes
whos ERmarkers % this is from section 15, mouse data

tempList = [diffGenes; ERmarkers'; Bmarkers'];
[a, b] = ismember(tempList, bloodExpGenes);
[a, b] = ismember(diffGenes, bloodExpGenes);
[a, b] = ismember( bloodExpGenes, diffGenes);
[a, b] = ismember(gtexCluster.syms, bloodExpGenes);

% see file corrCluster_function_GTExExpGenes.m in
% cellTypeExpression
geneList = tempList(a);
smallMat = GTExFiveNets.nets(1).net01(b(a), b(a));

% geneList = diffGenes(a);
% smallMat = GTExFiveNets.nets(1).net01(logical(a), logical(a))

c = cluster(z, 'cutoff', 1.1547)
max(c)
objectCount = hist(c, unique(c));
inCs = find(objectCount >= 5)
inMarkers = ismember(c, inCs);

for i = 1:length(inCs)
    inM(i).list = geneList(c == inCs(i));
end
save('~/resultsAndFigures/secondProject/bloodCorrelatedMarkers_BfromMouseIncluded_Distance_07.mat', ...
     'inM')
inMB = inM;
clear inM
load('~/resultsAndFigures/secondProject/bloodCorrelatedMarkers.mat')

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% which clusters they are in 
for i = 1:length(inM)
    thisDiffGenes = inM(i).list;
    [a, b] = ismember(thisDiffGenes, gtexCluster.syms);
    sum(a)
    
    % get the genes
    cids = gtexCluster.cs1547(b(a));
    myBar = hist(cids, unique(cids));
    h = figure
    bar(myBar)        
    set(gca, 'XTick', [1:length(unique(cids))], 'XTickLabel', ...
             unique(cids))
    ct = cellTypes{i};
    title(sprintf('listC %d, geneC %d', i, sum(a)))
end

% how is their expression level (this is from part 9.) 
for i = 1:length(inM)
    cgenes = inM(i).list;
    cgc = length(cgenes)
    [a, b] = ismember(cgenes, [Ads.geneSyms', Bds.geneSyms']);
    sum(a)
    wholeMat = [Ads.mat' Bds.mat']';
    wholeSMat = [Ads.samat' Bds.samat']';
    myboxS = wholeSMat(b(a), :);
    mybox = wholeMat(b(a), :);
    h = figure
    subplot(2, 2, 1)
    boxplot(mybox)
    title(sprintf('listC %d, geneC %d', i, sum(a)))
    
    subplot(2, 2, 2)
    boxplot(zscore(mybox')')
    title(sprintf('listC %d, geneC %d', i, sum(a)))
    
    subplot(2, 2, 3)
    boxplot(myboxS)
    title(sprintf('listC %d, geneC %d', i, sum(a)))
    set(gca, 'XTick', [1:12], 'XTickLabels', groups)
    xtickangle(90)
    
    subplot(2, 2, 4)
    boxplot(zscore(myboxS')')
    title(sprintf('listC %d, geneC %d', i, sum(a)))
    set(gca, 'XTick', [1:12], 'XTickLabels', groups)
    xtickangle(90)
end
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


h = figure
plot(bloodExpMat(b, :))

smallMat = bloodExpMat(b(a), :);
sib = corr(smallMat');
h = figure
heatmap(sib, 'GridVisible', 'off')

myMat = bloodExpMat(b, :);
tic
Y = tsne(myMat, 'Algorithm', 'exact', 'Distance', 'cosine');
toc

h = figure
gscatter(Y(:,1), Y(:,2), gcolor)

gcolor = zeros(1, length(gtexCluster.syms));
gcolor(gtexCluster.cs1547 == 153) = 2;
gcolor(gtexCluster.cs1547 == 232) = 1;
gcolor(gtexCluster.cs1547 == 86) = 3;
gcolor(gtexCluster.cs1547 == 215) = 4;

s1 = Y(:, 1) < -40;
s2 = Y(:, 2) > -20;

s = (s1 + s2) == 2;
tempGenes = find(s);
myGenes = gtexCluster.syms(tempGenes);
h = figure
hist(gtexCluster.cs1547(tempGenes), ...
     unique(gtexCluster.cs1547(tempGenes)))

% 13. Blood TAN work - Cluster bloodtan network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 1
tissue = 'blood'

load('~/data/general/GPL570GemmaMapNEW.mat')
affyGeneSyms = gpl570.uniqueSymbols;
clear gpl570

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'};
affyExpMat = zeros(18494, 5);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    affyExpMat(:, t) = expGenesInd;
end

tissue = 'blood'
load( ['~/networks/tissues/' tissue '/' ...
       'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
thisExpGenes = logical(affyExpMat(:, t));
affyBloodNet = binNet(thisExpGenes, thisExpGenes);
affyFullNet = affyBloodNet + affyBloodNet';
affyBloodSyms = affyGeneSyms(thisExpGenes);

% 13.1 get the clusters and check the marker presence 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
affyBloodNet;
affyBloodSyms;
affyFullNet;
affyBND = sum(affyFullNet);

affySNet = affyFullNet(affyBND > 1, affyBND > 1);
affyFSyms = affyBloodSyms(affyBND > 1);

tic
affyOverlap = affySNet * affySNet;
toc

affyNormMatUp = zeros(size(affyOverlap));
affyMaxOs = diag(affyOverlap);
tic
for i = 1:size(affyOverlap, 1);
    i
    for j = i:size(affyOverlap, 1);
        a = affyMaxOs(i);
        b = affyMaxOs(j);
        tempMin = min(a, b);
        affyNormMatUp(i, j) = tempMin;
        % gtexNormMat(j, i) = tempMin;
    end
end
toc
affyNormMat = triu(affyNormMatUp, 1) + affyNormMatUp' + 1 - affySNet;

kado = diag(affyOverlap);
affyTop = affyOverlap ./ affyNormMat;

affyDist = 1 - affyTop;
affyVDist = squareform(tril(affyDist, -1));
 
z = linkage(affyVDist, 'average');
h = figure
dendrogram(z, 2000)

t = cluster(z, 'cutoff', 1.1547); % go for 1.1547 to get the 300 clusters

%t = cluster(z, 'cutoff', .3, 'criterion', 'distance');
max(t)
h = figure
hist(t, max(t))
book = hist(t, max(t)); % count of members in each cluster
affyCluster.cs1547 = t;
affyCluster.syms = affyFSyms;
save('~/resultsAndFigures/secondProject/affyClusters_blood.mat', ...
     'affyCluster')
sum(book >= 20)
affyInClusterIds = find(book>=20);

% 13.2 get the inter and intra cluster density for affy stuff
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ccount = length(inClusterIds);
affyOverlapDensity = zeros(ccount, ccount);
for i = 1:ccount
    i
    geneSet1 = gtexCluster.syms(gtexCluster.cs1547 == inClusterIds(i)); 
    [a1, b1] = ismember(geneSet1, affyBloodSyms);
    sum(a1)
    affyOverlapDensity(i, i) = sum(sum(affyFullNet(b1(a1), b1(a1))))/(sum(a1) *(sum(a1)-1));
    for j = i+1:ccount
        geneSet2 = gtexCluster.syms(gtexCluster.cs1547 == inClusterIds(j)); 
        [a2, b2] = ismember(geneSet2, affyBloodSyms);
        chunk = affyFullNet(b1(a1), b2(a2));
        affyOverlapDensity(i, j) = sum(sum(chunk))/(sum(a1) *(sum(a2)));
    end
end

fullaod = (affyOverlapDensity + triu(affyOverlapDensity, ...
                                                  1)');

% plotting the interclusterovlerap for TAN
h = figure
h = figure('units', 'centimeters', 'position', [0,0, 25, 25])
fullaod(fullaod >.4) = .4;
heatmap(clusterLabels(outperm), clusterLabels(outperm), fullaod(outperm, ...
                                                  outperm), ...
        'GridVisible', 'off')
% removing clusters with low exps
heatmap(cinfo.clusterLabels(selectedCs), cinfo.clusterLabels(selectedCs), fullaod(selectedCs, ...
                                                  selectedCs), ...
        'GridVisible', 'off')

colormap(1-copper)
figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%sTANoverlap_norms_', figFolder);
file = sprintf('%sTANoverlap_norms_cfiltered', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

heatmap(inClusterIds(outperm(1:54)), inClusterIds(outperm(1:54)), fullaod(outperm(1:54), outperm(1:54)))

t = 1
gtexGeneSyms = GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(t).expGenes);
gtexNet = GTExFiveNets.nets(t).net005;
gtexND = sum(gtexNet + gtexNet');

gtexSNet = gtexNet(gtexND > 2, gtexND > 2);
gtexsfNet = gtexSNet + gtexSNet';
gtexFSyms = gtexGeneSyms(gtexND >2);

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

% save the two clusters
gtexCluster.cs1547 = t;
gtexCluster.cs1546 = t;
gtexCluster.syms = gtexFSyms;
save('~/resultsAndFigures/secondProject/gtexClusters_blood.mat', ...
     'gtexCluster')
load('~/resultsAndFigures/secondProject/gtexClusters_blood.mat')

% 14. expression level of genes in datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meanExps = zeros(1, length(gtexCluster.cs1547));
cgroups = zeros(1, length(gtexCluster.cs1547));
meanClusterExps = zeros(1, 63);
varClusterExps = zeros(1, 63);
for i = 1: length(inClusterIds)
    cid = inClusterIds(outperm(i));
    thisGeneSet = gtexCluster.cs1547 == cid;
    cgenes = gtexCluster.syms(thisGeneSet);
    [a, b] = ismember(cgenes, bloodExpGenes);
    meanExps(thisGeneSet) = mean(bloodExpMat(b, :)');
    cgroups(thisGeneSet) = i;
    meanClusterExps(i) = mean(mean(bloodExpMat(b, :)'));
    varClusterExps(i) = var(mean(bloodExpMat(b, :)'));
end

h = figure('units', 'centimeters', 'position', [0,0, 30, 12])
boxplot(meanExps, cgroups, 'PlotStyle', 'compact')
set(gca, 'XTick', [1:64], 'XTickLabels',[0, clusterLabels(outperm)])
xtickangle(90)
xlim([1 65])
title('expression level of genes in GTExBulkBlood clusters')

figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%sgeneExpLevelsInClusters', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 15. blood mouse data for Erythrocyte markers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 15.1 getting the platform files
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fid = fopen(['/home/mfarahbod/data/affyArray/bloodMouseCellType/' ...
             '2023_GSE6506_expmat.data.txt'])
gpl = textscan(fid, [repmat('%s', 1, 28)], 'Delimiter', '\t', ...
               'Headerlines',7);

%GPL96 probe ID's
gpl1261.probeID = gpl{1};

%GPL96 gene Symbols
gpl1261.geneSymbol = gpl{3};
gpl1261.uniqueGeneSymbols = unique(gpl{3});

%GPL96 gene names
gpl1261.geneNames = gpl{3};

save(['~/data/affyArray/bloodMouseCellType/gpl1261.mat'], 'gpl1261')

% 15.2 load the dataset 
fileName = ['/home/mfarahbod/data/affyArray/' ...
               'bloodMouseCellType/GSE6506_series_matrix.txt'];
fid = fopen(fileName)
linum = 75
header = textscan(fid, '%s', 1, 'delimiter', '\n', 'headerlines', ...
                  linum)

%getting the header:
tempHeader = header{1}{1};
headers = strsplit(tempHeader, '\t');
sampleCount = length(headers) - 1;

%getting the GOEID from a string. 
GEOID = cell(1, sampleCount);
expression = 'GSM[0123456789]+';
for j = 1:sampleCount
    GEOID(j) = regexp(headers{j+1}, expression, 'match');
end

%I need length of header for the sampleCount
fileName = ['/home/mfarahbod/data/affyArray/' ...
               'bloodMouseCellType/GSE6506_series_matrix.txt'];
fid = fopen(fileName)
linum = 76
data = textscan(fid, [repmat('%s', 1,1) repmat('%f', 1, ...
                                               sampleCount)], ...
                'Headerlines', linum, 'Delimiter', '\t');
fclose(fid);

gCount = length(gpl1261.uniqueGeneSymbols)
pCount = length(data{1})
dataExpr = cell2mat(data(2:21));
dataExpr(end, :) = [];

for j = 1:length(dataExpr)
    if(strcmp(data{1}{j}(1), '"'))
        data{1}{j} = data{1}{j}(2:end);
    end
    if(strcmp(data{1}{j}(end), '"'))
        data{1}{j} = data{1}{j}(1:(end-1));
    end
end

[a, b] = ismember(data{1}, gpl1261.probeID);
geneSymbols = gpl1261.geneSymbol(b(a)); %gene ID 22k
myProbes = data{1}(a); %probe names - 34k
myDataExpr = dataExpr(a ,:); %probe values
pCount = length(myProbes)

geneMap = containers.Map(gpl1261.uniqueGeneSymbols, [1:gCount]);
probeMap = containers.Map(myProbes, geneSymbols);%  NCBI geneID

probeDataSet.mat = myDataExpr;
probeDataSet.GEOID = GEOID;

% % on data, but on the platform
divMat = zeros(gCount, sampleCount);

%    adding the maps
tic
for j = 1:pCount
    sib = values(geneMap, values(probeMap, myProbes(j)));
    divMat(sib{1}, :) = divMat(sib{1}, :) + 1;
end
toc

%getting the porbe list for each gene
%max probe count for each gene
probeInd = zeros(gCount, max(max(divMat))); 
for j = 1:pCount
    sib = values(geneMap, values(probeMap, myProbes(j)));
    temp = find(probeInd(sib{1}, :) == 0);
    pInd = temp(1);
    probeInd(sib{1}, pInd) = j;
end

probeMean = mean(probeDataSet.mat, 2);
a = mean(probeMean)
newGeneExpr = zeros(gCount, sampleCount);
affectedGeneCount = 0;

q = quantile(probeDataSet.mat(:), 2)
thr = q(1);

% for each gene
for j = 1: size(probeInd, 1)
    pCount = sum(probeInd(j,:)>0);
    %if there are any with exp > thr, get the mean of bigger ones
    largeProbes = find(probeMean(probeInd(j,1:pCount))> thr);

    if(length(largeProbes) > 0)
        largeProbesInd = probeInd(j, largeProbes);
        newGeneExpr(j,:) = mean(probeDataSet.mat(largeProbesInd, ...
                                                   :), 1);
        
        if(pCount > length(largeProbes))
            affectedGeneCount = affectedGeneCount + 1;
        end
    else
        newGeneExpr(j, :) = mean(probeDataSet.mat(probeInd(j, ...
                                                          1:pCount), ...
                                                    :), 1);
    end
end

finalMat = newGeneExpr;
dataSet.mat = finalMat;
dataSet.GEOID = probeDataSet.GEOID;
dataSet.geneSyms = gpl1261.uniqueGeneSymbols;

save([dataFolder 'GSE6506.mat'], 'dataSet');
load([dataFolder 'GSE6506.mat']);

% >> get the homologue gene lists
file = '~/data/cellTypeVarianceFiles/humanAndMouse.txt'
fid = fopen(file)

linum = 0;
homoData = textscan(fid, [repmat('%d', 1, 3), repmat('%s', 1, 3)], ...
                'Headerlines', linum, 'Delimiter', '\t');

mouseGenes = cell(1, length(bloodExpGenes));
mouseOrt = zeros(1, length(bloodExpGenes));
homo = homoData{2} == 10090;
for i = 1:length(bloodExpGenes)
    i
    [a, b] = ismember(bloodExpGenes(i), homoData{4});
    mouseGenes(i) = {'nan'};
    if(a)
        [ai, bi] = ismember(homoData{1}, homoData{1}(b));
        mouseOrt(i) = sum((ai + homo) == 2);
        if (mouseOrt(i) == 1)
            mouseGenes(i) = homoData{4}((ai + homo) == 2);
        end
    end
end

% extract the "human dataset"
% mm is mouse ortholog for bloodExpGenes 
[mm, b] = ismember(mouseGenes, dataSet.geneSyms); 
dataSet.hgeneSyms = mouseGenes(mm);
dataSet.hhomoGeneSyms = bloodExpGenes(mm);
dataSet.hMat = dataSet.mat(b(mm), :);
save([dataFolder 'GSE6506_homoMapped.mat'], 'dataSet');
load([dataFolder 'GSE6506_homoMapped.mat']);

ctDataSet = load([dataFolder 'GSE6506_homoMapped.mat']);

% >>>> get groups of markers: Erythrocytes, Granulocytes,
% Monocyte, GM, B, BCN, C, CN, N

comTable = ones(9, 18) * -1;
% Erythrocytes
comTable(1, [17, 18]) = 1; 

% Granulocytes
comTable(2, [15, 16]) = 1; 

% Monocyte
comTable(3, [13, 14]) = 1; 
comTable(3, [15, 16]) = 0; 

% GM
comTable(4, [13:16]) = 1; 

% B 
comTable(5, [11, 12]) = 1; 

% BNC
comTable(6, [1:12]) = 1; 

% C
comTable(7, [3:10]) = 1; 

% CN
comTable(8, [1:10]) = 1; 

% N
comTable(9, [1:2]) = 1; 

h = figure
heatmap(comTable)

% getting the markerList
hGeneCount = length(dataSet.hMat);
mds = dataSet.hMat(:, 3:end);
for i = 1:9
    ns = find(comTable(i, :) == -1); % columns to be subtracted 
    nCount = length(ns) % count of them
    rs = find(comTable(i, :) == 1);
    rCount = length(rs) % count of main columns
    rSumm = zeros(hGeneCount, rCount);
    rExp = sum((mds(:, rs)' > 8)) == 2;
    for k = 1:rCount % for each of the main columns
        rmat = zeros(hGeneCount, nCount);
        for j = 1:nCount
            rmat(:, j) = mds(:, rs(k)) - mds(:, ns(j));
        end
        rSumm(:, k) = sum(rmat'>= 2) == nCount;
    end
    homoMarkerList{i} = dataSet.hhomoGeneSyms( sum([rSumm, rExp']') == 3);
end

% get the Erythrocyte markers:
remains1 = zeros(length(dataSet.hMat), 18);
remains2 = zeros(length(dataSet.hMat), 18);
for i = 1:18
    remains1(:, i) = dataSet.hMat(:, 19) - dataSet.hMat(:, i);
    remains2(:, i) = dataSet.hMat(:, 20) - dataSet.hMat(:, i);
end

r1l = remains1 >= 2;
r2l = remains2 >= 2;
e1 = dataSet.hMat(:, 19) > 8;
e2 = dataSet.hMat(:, 20) > 8;

s1 = sum(r1l');
s2 = sum(r2l');
sum((s1 + s2 + e1' + e2') == 38)
ERs = (s1 + s2 + e1' + e2') == 38;

h = figure
boxplot(dataSet.hMat((s1 + s2) == 36, :))

ERmarkers = dataSet.hhomoGeneSyms(ERs);

% Bmarkers
remains1 = zeros(length(dataSet.hMat), 18);
remains2 = zeros(length(dataSet.hMat), 18);
inds = [1:12, 15:20];
for i = 1:18
    remains1(:, i) = dataSet.hMat(:, 13) - dataSet.hMat(:, inds(i));
    remains2(:, i) = dataSet.hMat(:, 14) - dataSet.hMat(:, inds(i));
end

r1l = remains1 >= 1.5;
r2l = remains2 >= 1.5;
e1 = dataSet.hMat(:, 13) > 8;
e2 = dataSet.hMat(:, 14) > 8;

s1 = sum(r1l');
s2 = sum(r2l');
sum((s1 + s2 + e1' + e2') == 38)
Bs = (s1 + s2 + e1' + e2') == 38;

h = figure
boxplot(dataSet.hMat((s1 + s2) == 36, :))

Bmarkers = dataSet.hhomoGeneSyms(Bs);

% extract the homo dataset genes
% which of the mouse genes are in the list
% from these genes, how many of them have human, and if they do,
% get the equivalent human annotation
% humanGenes = cell(1, length(dataSet.geneSyms));
% humanOrt = zeros(1, length(dataSet.geneSyms));
% homo = homoData{2} == 9606;
% for i = 1:length(dataSet.geneSyms)
%     [a, b] = ismember(dataSet.geneSyms(i), homoData{4});
%     if(a)
%         [ai, bi] = ismember(homoData{1}, homoData{1}(b));
%         if (sum((ai + homo) == 2))
%             human
%         end
%     end
% end
    
% 16. first report plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% >>> 16.1. expression of genes in each cluster, the R2 for each cluster 
% >>> 16.2 heatmap of the clusters r2 - GTEx
% >>> 16.3 heatmap of the clusters CTProfiles
% >>> 16.4 plot the CTs for mouse and clusters (from 11.5)
% >>> 16.5 plot for Rs correlation in different datasets
% >>> 16.6 plot for cluster R2s in different datasets
% >>> 16.7 correlation of cluster R2s 
% >>> 16.8 correlation of CT profiles
% >>> 16.9 expression of final marker genes

% load blood dataset, genes
load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat')
bloodNet = GTExFiveNets.nets(1);
bloodExpGenes = bloodDS.genes(bloodNet.expGenes);
clear GTExFiveNets

% get bloodDS 
[a, b] = ismember(bloodExpGenes, bloodDS.genes);
bloodExpMat = log2(bloodDS.mat(b(a), :)+1);

% load the Rs
load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240_Eryth177_Eryth163.mat']) 
regResM = result; 
clear result

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240.mat'])
regRes = result;
clear result

% load the clusters
load('~/resultsAndFigures/secondProject/gtexClusters_blood.mat')

% >>> build the cluster ID-gc labels:
clusterLabels = cell(1, length(inClusterIds));
for i = 1:length(inClusterIds)
   clusterLabels{i} = sprintf('ID%d - %d', inClusterIds(i), incmc(i));
end

% >>> 16.1. expression of genes in each cluster, the R2 for each cluster 
whos meanClusterExps % mean exp level in each cluster
whos varClusterExps % mean var in each cluster
whos avgCTP
whos outperm
whos inClusterIds
whos incr2;
incr2p = incr2(outperm);

book = [avgCTP, meanClusterExps', incr2'];
book(41, :) = [];
h = figure
heatmap(corr(book))
colormap(1- [1-gray' copper']')

% >>> 16.2 heatmap of the clusters r2 - GTEx
h = figure('units', 'centimeters', 'position', [0,0, 25, 12])
heatmap(clusterLabels(outperm), [1], incr2p);
colormap(1- [1-gray' copper']')
figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%sclusterR2Vals', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% >>> 16.3 heatmap of the clusters CTProfiles

grouporder = [1 8 3 4 5 9 10]
h = figure('units', 'centimeters', 'position', [0,0, 25, 12])
book = zscore(avgCTP(selectedCs, grouporder)');
heatmap(cinfo.clusterLabels(selectedCs), groups(grouporder), (book), ...
        'GridVisible', 'off');
colormap(1- [1-gray' copper']')
figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%sclusterCTprofiles_selectedCs', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


h = figure('units', 'centimeters', 'position', [0,0, 25, 12])
book = zscore(avgCTP(outperm, grouporder)');
colormap(1- [1-gray' copper']')
heatmap(clusterLabels(outperm), groups(grouporder), (book), ...
        'GridVisible', 'off');
h = figure('units', 'centimeters', 'position', [0,0, 25, 12])
medPlotMat = medCTP;
medPlotMat(medPlotMat > 9.5) = 9.5;
medPlotMat(medPlotMat < 5.5) = 5.5;

 h = figure
heatmap(clusterLabels(outperm([1:57, 59:end])), groups(grouporder), (medPlotMat(outperm([1:57, 59:end]), grouporder))', ...
        'GridVisible', 'off');

colormap(ocmap/256)
colormap(1- [1-gray' copper']')
colormap(1- [ copper']')

figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%sclusterCTprofiles_doubleZ_7ct', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% >>>>>
i = 1
h = figure('units', 'centimeters', 'position', [0,0, 25, 5])
heatmap(clusterLabels(outperm([1:57, 59:end])), groups(grouporder(i)), (medPlotMat(outperm([1:57, 59:end]), grouporder(i)))', ...
        'GridVisible', 'off');

colormap(bcmap/256)

figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%s_expCluster%d', figFolder, i);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')
% <<<<<<

% >>> 16.4 plot the CTs for mouse and clusters (from 11.5)
grouporder = [8, 10, 1 3 4 5 9 32]
plotMat = [avgCTP, avgCTPmouse];
h = figure('units', 'centimeters', 'position', [0,0, 25, 12])
heatmap([clusterLabels(outperm)], [groups(grouporder(1:7)), 'erythrocyte'] , zscore(plotMat(outperm, grouporder ...
                                                  )'), 'GridVisible', ...
        'off');
myCmap = 1- [1-gray' copper']';
colormap(1- [1-gray' copper']')
figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%sclusterCTprofiles_erythroIncluded_gridOff', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% report the correlation or plot them or whatever
book = zscore(plotMat(outperm, grouporder )')';
book(58, :) = [];
sib = corr(book(:, 1:7)');
h = figure
heatmap(sib)
colormap(1- [1-gray' copper']')
figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%sCTcorrs', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% >>> 16.5 plot for Rs correlation in different datasets
affy16Rs
affy16RsE
affy27ERs
affy27Rs
regResM
regRes
whos bloodExpGenes
whos affyGeneSyms

[a, b] = ismember(affyGeneSyms, bloodExpGenes);
sum(a)

allRvectors = zeros(sum(a), 6);
allRvectors(:, 1) = affy16Rs.result.r2(a);
allRvectors(:, 2) = affy16RsE.result.r2(a);
allRvectors(:, 3) = affy27Rs.result.r2(a);
allRvectors(:, 4) = affy27ERs.result.r2(a);
allRvectors(:, 5) = regRes.r2(b(a));
allRvectors(:, 6) = regResM.r2(b(a));

sib = corr(allRvectors);
h = figure
heatmap(sib)

[as, bs] = ismember(gtexCluster.syms, bloodExpGenes);
clustersR2 = regRes.r2(bs(as));
clustersR2M = regResM.r2(bs(as));
[a, b] = ismember(affyGeneSyms, gtexCluster.syms);
sum(a)

allRvectors = zeros(sum(a), 8);
allRvectors(:, 1) = affy16Rs.result.r2(a);
allRvectors(:, 2) = affy16RsE.result.r2(a);
allRvectors(:, 3) = affy27Rs.result.r2(a);
allRvectors(:, 4) = affy27ERs.result.r2(a);
allRvectors(:, 5) = clustersR2(b(a));
allRvectors(:, 6) = clustersR2M(b(a));
allRvectors(:,  7) = affy13Rs.result.r2(a);
allRvectors(:,  8) = affy77Rs.result.r2(a);
sib = corr(allRvectors(:, [1 3 5]));
sib = corr(allRvectors);
h = figure
heatmap({'affy16', 'affy27', 'GTExblood'}, {'affy16', 'affy27', 'GTExblood'}, sib)

figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%sR2corrs_clusterGenes_affy16_affy27_GTExblood', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

h = figure
scatter(allRvectors(:, 1), allRvectors(:, 3))
scatter(allRvectors(:, 5), allRvectors(:, 3))
scatter(regRes.r2(b(a)), affy77Rs.result.r2(a))

% >>> 16.6 plot for cluster R2s in different datasets

clusterR2s = zeros(3, length(inClusterIds))
geneInCount = zeros(1, length(inClusterIds));
for i = 1:length(inClusterIds)
    % get genes in the cluster
    i
    geneSet = gtexCluster.syms(gtexCluster.cs1547 == inClusterIds(outperm(i)));
              
    [a, b] = ismember(geneSet, bloodExpGenes);
    clusterR2s(3, i) = mean(regRes.r2(b(a)));
    
    [a, b] = ismember(geneSet, affyGeneSyms);
    clusterR2s(1, i) = mean(affy16Rs.result.r2(b(a)));
    clusterR2s(2, i) = mean(affy27Rs.result.r2(b(a)));
end

h = figure('units', 'centimeters', 'position', [0,0, 28, 8])
heatmap(clusterLabels(outperm), {'affy16', 'affy27', 'GTExblood'}, ...
        clusterR2s)
colormap(1- bone)

figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%sclusterR2heatmap_affy16_affy27_GTExblood', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% >>> 16.7 correlation of cluster R2s 
book = clusterR2s;
book(:, 58) = [];
[s, p] = corr(book')

h = figure
heatmap({'affy16', 'affy27', 'GTExblood'},{'affy16', 'affy27', 'GTExblood'}, ...
        s)

figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%sclusterR2Correlation_affy16_affy27_GTExblood', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% >>> 16.8 correlation of CT profiles
whos avgCTPwhole % it is 63 by 114 , which is the avg gene
                 % cell type expression

groups = {'B', 'CD14+', 'CD4+Te', 'CD4+T', 'CD8+T', 'IgGAmemB', ...
          'IgMmemB', 'Monocyte', 'NK', 'Neutrophils', 'pbm', 'ppb'}

labels = [repmat(5, 1, 4), repmat(4, 1, 14), repmat(3, 1, 6), ...
          repmat(9, 1, 15), repmat(1, 1, 7), repmat(6, 1, 4), ...
          repmat(7, 1, 4), repmat(12, 1, 3), repmat(11, 1, 4), ...
          repmat(8, 1, 36), repmat(2, 1, 12), repmat(10, 1, 5)]

[a, b] = ismember(labels, [8]);
myMat = avgCTPwhole(:, a);
myMat = avgCTPwhole;

[s, p] = corr(myMat');

plotMat = -log10(p(outperm, outperm));
plotMat = -log10(p(selectedCs, selectedCs));

h = figure('units', 'centimeters', 'position', [0,0, 25, 25])
heatmap(cinfo.clusterLabels(selectedCs), cinfo.clusterLabels(selectedCs), plotMat,...
        'GridVisible', 'off')
colormap(1-copper)
figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%savgCTPWhole_corr_selectedCs', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

colormap(1-copper)

sib = corr(avgCTP');
h = figure
heatmap(sib(outperm, outperm))
colormap(1- [1-gray' copper']')
colormap

% >>>> 16.9 
load(['~/resultsAndFigures/secondProject/' ...
      'bloodCorrelatedMarkers.mat'])
whos inM
% count of markers in inM
 
mCount = length(inM(1).list);
mCounts = zeros(1,15);
mCounts(1) = mCount;
markerLists = inM(1).list;
for i = 2:13
    mCount = mCount + length(inM(i).list);
    mCounts(i) = length(inM(i).list);
    markerLists = [markerLists; inM(i).list];
end

markerLists = [markerLists; Bmarkers'; ERmarkers'];
markerLists = [markerLists; m177List'; m163List']; % this is what I
                                                   % am plotting
mCounts(14) = length(m177List);
mCounts(15) = length(m163List);
inMlistOrder = [1 2 3 4 12 6 7 9 10 5 8 11 13]
 
mlabels = zeros(1, 240);
s = 1;
e = length(inM(1).list);
markerListsReO = inM(1).list;
mlabels(s:(s + e -1)) = 1;
s = e + 1;
for i = 2:13
    order = inMlistOrder(i);
    markerListsReO = [markerListsReO; inM(order).list];
    e = length(inM(order).list);
    mlabels(s:(s + e -1)) = i;
    s = s + e;
end
%markerListsReO = [markerListsReO; m177List'; m163List']; % this is
%what I
[a, b] = ismember(markerListsReO, [Ads.geneSyms', Bds.geneSyms']);
wholesMat = [Ads.samat' Bds.samat']';
smallExp = wholesMat(b(a), :);
grouporder = [1 8 3 4 5 9 10]

h = figure('units', 'centimeters', 'position', [0,0, 12, 30])
kado = smallExp(:, grouporder);
heatmap(zscore(kado')',...
        'GridVisible', 'off')
colormap(1- [1-gray' copper']')
figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%smarkerListExp_zscoreH', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

[a, b] = ismember(markerListsReO, [Ads.geneSyms', Bds.geneSyms']);
h = figure('units', 'centimeters', 'position', [0,0, 12, 30])
heatmap(mlabels(a)', 'GridVisible', 'off')
colormap(lines)
figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%smarkerListExp_zscoreH_mlabels', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


%17. affy blood dataset human bulk GSE16028
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('~/data/general/GPL570GemmaMapNEW.mat')
affyGeneSyms = gpl570.uniqueSymbols;
clear gpl570

load(['~/data/affyArray/tissues/blood/matFiles/geneExprMatGemmaMapBER/',...
      'GSE16028_BERevomed.mat'])
myNet = load(['~/networks/tissues/blood/singleNet/' ...
              'bloodNetExpThr0.6_GSE16028_0.005.mat'])

myFullNet = full(myNet.sparseSingleNet + myNet.sparseSingleNet');
myNDs = sum(myFullNet);
smallNet = myFullNet(myNDs >= 2, myNDs>=2);
myNetGenes = affyGeneSyms(myNDs >= 2);

affyInClusters = find(objectCount > 20);
affyGSE16028Clusters.c1547 = c;
affyGSE16028Clusters.genes = myNetGenes;
save([dataFolder 'affyGSE16028Clusters.mat'], ...
     'affyGSE16028Clusters');

% get the expression level of genes in thse clusters
myexps = zeros(1, length(gtexCluster.syms));
mygroups = zeros(1, length(gtexCluster.syms));
binds = ismember(gtexCluster.syms, affyGeneSyms);
meanExpsAffy = zeros(1, 63);
for i = 1:length(inClusterIds)
    thisID = inClusterIds(outperm(i));
    [a, b] = ismember(gtexCluster.syms(gtexCluster.cs1547 == ...
                                       thisID), ...
                      affyGeneSyms);
    sum(a)
    tempInds = ((gtexCluster.cs1547 == thisID) + binds') == ...
        2;
    sib = mean(dataSet.mat(b(a), :)');
    meanExpsAffy(i) = mean(sib);
    myexps(tempInds) = sib;
    mygroups(tempInds) = i;
end

h = figure
boxplot(myexps, mygroups)
set(gca, 'XTick', [1:64], 'XTickLabels',[0, clusterLabels(outperm)])
xtickangle(90)
    
% >>>>>>>>>> now get the R2 way 1: 

load('~/resultsAndFigures/secondProject/bloodDiffGenesInfo.mat') 
diff.cellTypes{:}
% cellTypes we are getting: 
tempSelected = [1, 2, 4:13, 22];

inGenes = diff.genes((sum(diff.mat(:, tempSelected)') > 0));
%inGenes = diff.genes;
length(inGenes)

[a, b]  = ismember(inGenes, affyGeneSyms);
markerExpMat = dataSet.mat(b(a), :);

markerLists = inGenes;

% >>>>>>>>>> way 2
markerLists = regResM.markers;
markerLists = cgenes;
[a, b] = ismember(markerLists, affyGeneSyms);
markerExpMat = dataSet.mat(b(a), :);

% getting the prediction
[w, score, latent, tsquared, explained, mu] = pca(markerExpMat', ...
                                                  'NumComponents', 10);

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

result.markers = markerLists;
result.coeff = coeffs;
result.featureMat = featureMat;
result.regOut = regOut;
result.r2 = r2;
save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_diffGenes_all_AffyDSGSE16028.mat'],'result', ...
     '-v7.3')

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240PErythro_AffyDSGSE16028.mat'],'result', ...
     '-v7.3')

affy16Rs = load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240_AffyDSGSE16028.mat'])
affy16RsE = load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240PErythro_AffyDSGSE16028.mat'])


affy77Rs = load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240_AffyDSGSE7753.mat'])
affy13Rs = load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240_AffyDSGSE13849.mat'])


meanRsAffy16 = zeros(1, length(inClusterIds));
for i = 1:length(inClusterIds)
    thisID = inClusterIds(outperm(i));
    [a, b] = ismember(gtexCluster.syms(gtexCluster.cs1547 == ...
                                       thisID), ...
                      affyGeneSyms);
    sum(a)
    meanRsAffy16(i) = mean(affyRs.result.r2(b(a)));
end

h = figure
heatmap(clusterLabels(outperm), [1, 2] ,[meanRs', meanRs2']')

pMat = zeros(length(inClusterIds), 7);
geneInCount = zeros(1, length(inClusterIds));
for i = 1:length(inClusterIds)
    % get genes in the cluster
    i
    geneSet = gtexCluster.syms(gtexCluster.cs1547 == inClusterIds(i));
              
    %>>> doing it with just Ads
    % [a, b] = ismember(geneSet, Ads.geneSyms);
    % % get them in the expression data
    % smallExp = (Ads.samat(b(a), :));
    % whos smallExp
    % %inGenes = sum(smallExp' > 8) > 3;
    
    % % the average 
    % %avgCTP(i, :) = mean(smallExp(inGenes,:));
    % avgCTP(i, :) = mean(zscore(smallExp')');
    
    %>>> doing it with just Ads and Bds
    [a, b] = ismember(geneSet, affyGeneSyms);
    pMat(i, :) = mean(result.coeff(b(a), 2:end));
end
myAffyFMat = pMat(outperm, :)';
h= figure
kado = corr(pMat(outperm, :)');
heatmap(clusterLabels(outperm), clusterLabels(outperm), kado)
colormap(1- [1-gray' copper']')

kado = [myAffyFMat'; myGTExFMat];

sib = corr(kado');
h= figure
heatmap(sib)


% >>> get the cluster overlap
ccount = length(inClusterIds);
affyFullNet = myNet.sparseSingleNet + myNet.sparseSingleNet' + 0;
affyOverlapDensity16 = zeros(ccount, ccount);
for i = 1:ccount
    i
    geneSet1 = gtexCluster.syms(gtexCluster.cs1547 == inClusterIds(i)); 
    [a1, b1] = ismember(geneSet1, affyGeneSyms);
    sum(a1)
    affyOverlapDensity16(i, i) = sum(sum(affyFullNet(b1(a1), b1(a1))))/(sum(a1) *(sum(a1)-1));
    for j = i+1:ccount
        geneSet2 = gtexCluster.syms(gtexCluster.cs1547 == inClusterIds(j)); 
        [a2, b2] = ismember(geneSet2, affyGeneSyms);
        chunk = affyFullNet(b1(a1), b2(a2));
        affyOverlapDensity16(i, j) = sum(sum(chunk))/(sum(a1) *(sum(a2)));
    end
end

fullaod = (affyOverlapDensity16 + triu(affyOverlapDensity16, ...
                                                  1)');

% plotting the interclusterovlerap for TAN
h = figure
h = figure('units', 'centimeters', 'position', [0,0, 25, 25])
heatmap(clusterLabels(outperm), clusterLabels(outperm), fullaod(outperm, ...
                                                  outperm), ...
        'GridVisible', 'off')
colormap(1-copper)
figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%sTANoverlap', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


% 18. hk genes in the clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('~/data/general/GPL570GemmaMapNEW.mat')
affyGeneSyms = gpl570.uniqueSymbols;
clear gpl570

load('~/data/general/hkgInd.mat')
myhkGenes = affyGeneSyms(hkgInd);

hkpor = zeros(1, length(inClusterIds));
for i = 1:length(inClusterIds)
    [a, b] = ismember(gtexCluster.syms(gtexCluster.cs1547 == ...
                                       inClusterIds(outperm(i))), ...
                      myhkGenes);
    hkpor(i) = sum(a)/length(a);
end
h = figure('units', 'centimeters', 'position', [0,0, 12, 35])
heatmap(clusterLabels(outperm), [1], hkpor)
figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%shk_percent', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


% 19. that regression thing : are the Rs from markers more
% correlated than Rs from random groups of genes. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get R2s for marker (this if from section 6.)

% >> for GTEx
[a, b] = ismember(markerLists, bloodExpGenes);
markerR2s = regRes.r2(b);
qsa = quantile(markerR2s, 30);

[a, b] = ismember(gtexCluster.syms, bloodExpGenes);
clusterGeneR2s = regRes.r2(b);
qsb = quantile(clusterGeneR2s, 30)
h = figure
plot(qsa, qsb, '.')
xlabel('markers')
ylabel('cluster')

xlim([.5 1])
ylim([.5 1])
hold on
plot([.5 1], [.5 1])

% >> for affy
[a, b] = ismember(markerLists, affyGeneSyms);
sum(a)
markerR2sAffy = affyRs.result.r2(b(a));
qsa = quantile(markerR2sAffy, 30);

[a, b] = ismember(gtexCluster.syms, affyGeneSyms);
sum(a)
clusterGeneR2sAffy = affyRs.result.r2(b(a));
qsb = quantile(clusterGeneR2sAffy, 30)
h = figure
plot(qsa, qsb, '.')
xlabel('markers')
ylabel('cluster')

xlim([0 1])
ylim([0 1])
hold on
plot([0 1], [0 1])


% now, for selection of 10 sets of random genes, from each dataset,
% with similar expression levels, we get sets of R2s for each
% dataset. Now those R2s, are they as highly correlated or not. 


%20. affy blood dataset human bulk GSE27562
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['~/data/affyArray/tissues/blood/matFiles/geneExprMatGemmaMapBER/',...
      'GSE27562_BERevomed.mat'])
myNet = load(['~/networks/tissues/blood/singleNet/' ...
              'bloodNetExpThr0.6_GSE27562_0.005.mat'])

dataSet.mat = dataSet.mat(:, [1:10, 52:72]);

myFullNet = full(myNet.sparseSingleNet + myNet.sparseSingleNet');
myNDs = sum(myFullNet);
smallNet = myFullNet(myNDs >= 2, myNDs>=2);
myNetGenes = affyGeneSyms(myNDs >= 2);

affyInClusters = find(objectCount > 20);
affyGSE16028Clusters.c1547 = c;
affyGSE16028Clusters.genes = myNetGenes;
save([dataFolder 'affyGSE16028Clusters.mat'], ...
     'affyGSE16028Clusters');

% get the expression level of genes in thse clusters
myexps = zeros(1, length(gtexCluster.syms));
mygroups = zeros(1, length(gtexCluster.syms));
binds = ismember(gtexCluster.syms, affyGeneSyms);
meanExpsAffy = zeros(1, 63);
for i = 1:length(inClusterIds)
    thisID = inClusterIds(outperm(i));
    [a, b] = ismember(gtexCluster.syms(gtexCluster.cs1547 == ...
                                       thisID), ...
                      affyGeneSyms);
    sum(a)
    tempInds = ((gtexCluster.cs1547 == thisID) + binds') == ...
        2;
    sib = mean(dataSet.mat(b(a), :)');
    meanExpsAffy(i) = mean(sib);
    myexps(tempInds) = sib;
    mygroups(tempInds) = i;
end

h = figure
boxplot(myexps, mygroups)
set(gca, 'XTick', [1:64], 'XTickLabels',[0, clusterLabels(outperm)])
xtickangle(90)
    
% now get the R2 

% >>>>>>>>>>>>> way 1 all diffMat marker genes
load('~/resultsAndFigures/secondProject/bloodDiffGenesInfo.mat') 
diff.cellTypes{:}
% cellTypes we are getting: 
tempSelected = [1, 2, 4:13, 22];

inGenes = diff.genes((sum(diff.mat(:, tempSelected)') > 0));
%inGenes = diff.genes;
length(inGenes)

[a, b]  = ismember(inGenes, affyGeneSyms);
markerExpMat = dataSet.mat(b(a), :);

markerLists = inGenes;


% >>>>>>> way 2 selected marker genes

markerLists = regResM.markers;
markerLists = regRes.markers;
[a, b] = ismember(markerLists, affyGeneSyms);
markerExpMat = dataSet.mat(b(a), :);

% getting the prediction
[w, score, latent, tsquared, explained, mu] = pca(markerExpMat', ...
                                                  'NumComponents', 10);

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

result.markers = markerLists;
result.coeff = coeffs;
result.featureMat = featureMat;
result.regOut = regOut;
result.r2 = r2;
save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_diffGenes_all_AffyDSGSE27562.mat'],'result', ...
     '-v7.3')

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240PErythro_AffyDSGSE27562.mat'],'result', ...
     '-v7.3')

save(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240_AffyDSGSE27562.mat'],'result', ...
     '-v7.3')

affy27Rs = load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240_AffyDSGSE27562.mat']);

affy27ERs = load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240PErythro_AffyDSGSE27562.mat']);


meanRsAffy27 = zeros(1, length(inClusterIds));
for i = 1:length(inClusterIds)
    thisID = inClusterIds(outperm(i));
    [a, b] = ismember(gtexCluster.syms(gtexCluster.cs1547 == ...
                                       thisID), ...
                      affyGeneSyms);
    sum(a)
    meanRsAffy27(i) = mean(result.r2(b(a)));
end

h = figure
heatmap(clusterLabels(outperm), [1], meanRsAffy27)


[a, b] = ismember(markerLists, affyGeneSyms);
sum(a)
markerR2sAffy2 = result.r2(b(a));
qsa = quantile(markerR2sAffy2, 30);

[a, b] = ismember(gtexCluster.syms, affyGeneSyms);
sum(a)
clusterGeneR2sAffy2 = result.r2(b(a));
qsb = quantile(clusterGeneR2sAffy2, 30)
h = figure
plot(qsa, qsb, '.')
xlabel('markers')
ylabel('cluster')

xlim([0 1])
ylim([0 1])
hold on
plot([0 1], [0 1])


% 21. correlation of different celltypes for clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
whos avgCTPwhole % it is 63 by 114 , which is the avg gene
                 % cell type expression

groups = {'B', 'CD14+', 'CD4+Te', 'CD4+T', 'CD8+T', 'IgGAmemB', ...
          'IgMmemB', 'Monocyte', 'NK', 'Neutrophils', 'pbm', 'ppb'}

labels = [repmat(5, 1, 4), repmat(4, 1, 14), repmat(3, 1, 6), ...
          repmat(9, 1, 15), repmat(1, 1, 7), repmat(6, 1, 4), ...
          repmat(7, 1, 4), repmat(12, 1, 3), repmat(11, 1, 4), ...
          repmat(8, 1, 36), repmat(2, 1, 12), repmat(10, 1, 5)]

[a, b] = ismember(labels, [8]);
myMat = avgCTPwhole(:, a);
myMat = avgCTPwhole;

[s, p] = corr(myMat');

plotMat = -log10(p(outperm, outperm));
plotMat = -log10(p(outperm, outperm));
h = figure
heatmap(plotMat, 'GridVisible', 'off')
colormap(1-copper)

sib = corr(avgCTP');
h = figure
heatmap(sib(outperm, outperm))
colormap(1- [1-gray' copper']')
colormap

% 22. correlation of R2, CT profiles and variance for affy for
% blood genes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('~/data/general/GPL570GemmaMapNEW.mat')
affyGeneSyms = gpl570.uniqueSymbols;
clear gpl570

affy16Rs = load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240_AffyDSGSE16028.mat'])
load(['~/data/affyArray/tissues/blood/matFiles/geneExprMatGemmaMapBER/',...
      'GSE16028_BERevomed.mat'])
ds16 = dataSet;
clear dataSet

affy27Rs = load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240_AffyDSGSE27562.mat']);
load(['~/data/affyArray/tissues/blood/matFiles/geneExprMatGemmaMapBER/',...
      'GSE27562_BERevomed.mat'])
ds27 = dataSet;
ds27.mat = ds27.mat(:, [1:10, 52:72]);
clear dataSet

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240.mat'])
regRes = result;
load('~/resultsAndFigures/secondProject/bloodDS.mat') % bloodDS
bloodExpMat = log2(bloodDS.mat+1);
load('~/resultsAndFigures/secondProject/gtexClusters_blood.mat') % gtexCluster
load('~/resultsAndFigures/secondProject/bloodClusterInfo.mat')  % cinfo 
selectedCs = cinfo.outperm([1:54, 57, 59, 61]); % removing clusters
                                                % with low exp

[a, b] = ismember(markerListsReO, [Ads.geneSyms', Bds.geneSyms']);
wholesMat = [Ads.samat' Bds.samat']';
wholeMat = [Ads.mat' Bds.mat']';

% getting genes in relevant clusters
[a, b] = ismember(gtexCluster.cs1547, cinfo.inIDs(selectedCs));
cgenes = gtexCluster.syms(a);
[a, b] = ismember(cgenes, [Ads.geneSyms', Bds.geneSyms']);
fcgenes = cgenes(a);
cts = wholesMat(b(a), :);
%cts = wholeMat(b(a), :);
vcts = var(cts');

[a, b] = ismember(fcgenes, bloodExpGenes);
bvar = var(bloodExpMat(b(a), :)');
[a, b] = ismember(fcgenes, bloodExpGenes);
brs = result.r2(b(a));
% brs = rg.r2(b(a)); % this is from 24, for another R with
% different marker set
corrMat = [bvar; vcts; brs'];
kadogtex = corr(corrMat');
h = figure('units', 'centimeters', 'position', [0,0, 35, 12])
subplot(1, 3, 1)
heatmap({'bulkvar', 'ctvar', 'r2s'}, {'bulkvar', 'ctvar', 'r2s'}, ...
        kadogtex)
title('gtex')

% for affy27
[a, b] = ismember(fcgenes, affyGeneSyms);
affyvcts = vcts(a);
av27 = var(ds27.mat(b(a), :)');
kado = corr(av27', affyvcts');
a27rs = affy27Rs.result.r2(b(a));
% a27rs = r27.r2(b(a));
corrMat = [av27; affyvcts; a27rs'];
kadoaffy27 = corr(corrMat');
subplot(1, 3, 2)
heatmap({'bulkvar', 'ctvar', 'r2s'}, {'bulkvar', 'ctvar', 'r2s'}, kadoaffy27)
title('affy27')

% for affy16
[a, b] = ismember(fcgenes, affyGeneSyms);
affyvcts = vcts(a);
av16 = var(ds16.mat(b(a), :)');
kado = corr(av16', affyvcts');
a16rs = affy16Rs.result.r2(b(a));
% a16rs = r16.r2(b(a));
corrMat = [av16; affyvcts; a16rs'];
kadoaffy16 = corr(corrMat');
subplot(1, 3, 3)
heatmap({'bulkvar', 'ctvar', 'r2s'}, {'bulkvar', 'ctvar', 'r2s'}, ...
        kadoaffy16)
title('affy16')

figFolder = ['~/resultsAndFigures/secondProject/bloodReportFigures/']
file = sprintf('%scorrelations_bulkVar_CTVar_R2', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


% correlation of R2 between datasets
r27; % rs are from 24
r16;
rg;
[a, b] = ismember(affyGeneSyms, gtexCluster.syms);

myGenes = affyGeneSyms(a);
[ag, bg] = ismember(myGenes, bloodExpGenes);

sib = [r27.r2(a), r16.r2(a), rg.r2(bg(ag))];
kado = corr(sib)

% 23. plot for distribution of R2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('~/resultsAndFigures/secondProject/bloodDS.mat')

load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat')
bloodNet = GTExFiveNets.nets(1);
bloodExpGenes = bloodDS.genes(bloodNet.expGenes);
clear GTExFiveNets


load('~/data/general/GPL570GemmaMapNEW.mat')
affyGeneSyms = gpl570.uniqueSymbols;
clear gpl570

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240.mat'])
regRes = result;

load('~/resultsAndFigures/secondProject/gtexClusters_blood.mat') % gtexCluster

[a, b] = ismember(bloodExpGenes, gtexCluster.syms);

myvals = regRes.r2(a);
myvals = rg.r2(a);


[f, xi] = ksdensity(myvals);
h = figure;
plot(xi, f)
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sbloodR2dist_ksdensity_r2DiffGenes_SelectedAll', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


h = figure;
histogram(myvals, 100, 'DisplayStyle', 'stairs')
set(gca, 'YTick', [0:50:450], 'XTickLabels',[0:50:450]/length(myvals))

figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sbloodR2dist_histStairsV01', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

affy16Rs = load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240_AffyDSGSE16028.mat'])
load(['~/data/affyArray/tissues/blood/matFiles/geneExprMatGemmaMapBER/',...
      'GSE16028_BERevomed.mat'])
ds27 = dataSet;
clear dataSet

affy27Rs = load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_240_AffyDSGSE27562.mat']);
load(['~/data/affyArray/tissues/blood/matFiles/geneExprMatGemmaMapBER/',...
      'GSE27562_BERevomed.mat'])
ds16 = dataSet;
ds16.mat = ds16.mat(:, [1:10, 52:72]);
clear dataSet
h = figure;
histogram(myvals, 100, 'DisplayStyle', 'stairs')
set(gca, 'YTick', [0:50:450], 'XTickLabels',[0:50:450]/ ...
         length(myvals))
hold on

[a, b] = ismember( affyGeneSyms, gtexCluster.syms);
mv16 = r16.r2(a);
mv16 = affy16Rs.result.r2(a);
histogram(mv16, 100, 'DisplayStyle', 'stairs')

mv27 = r27.r2(a);
mv27 = affy27Rs.result.r2(a);
histogram(mv27, 100, 'DisplayStyle', 'stairs')

set(gca, 'YTick', [0:50:450], 'XTickLabels',[0:50:450]/length(myvals))

figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sbloodR2dist_histStairsAllThree', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


[f, xi] = ksdensity(myvals);
h = figure;
plot(xi, f)
hold on
[f, xi] = ksdensity(mv16);
plot(xi, f)
[f, xi] = ksdensity(mv27);
plot(xi, f)
legend([1, 2, 3])
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sbloodR2dist_ksdensityAllThree_DiffMarkers_selectedAll', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


% 24. R2 for clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('~/data/general/GPL570GemmaMapNEW.mat')
affyGeneSyms = gpl570.uniqueSymbols;
clear gpl570


load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_diffGenes_all_AffyDSGSE27562.mat' ...
      ''])
r27 = result;
clear result

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_diffGenes_all_AffyDSGSE16028.mat'])
r16 = result;
clear result

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodRegression_mSelected_diffGenes_selectedTypes.mat'])
rg = result;
clear result

affyGeneSyms; % affy genes
bloodExpGenes; % gtex genes

[a, b] = ismember(affyGeneSyms, bloodExpGenes);

selectedCs = cinfo.outperm([1:54, 57, 59, 61]); % removing clusters

meanRsAffy27 = zeros(1, length(selectedCs));
meanRsAffy16 = zeros(1, length(selectedCs));
meanRs = zeros(1, length(selectedCs));
for i = 1:length(selectedCs)
    thisID = selectedCs(i);
    [a, b] = ismember(gtexCluster.syms(gtexCluster.cs1547 == ...
                                       thisID), ...
                      affyGeneSyms);
    sum(a)
    meanRsAffy27(i) = mean(r27.r2(b(a)));
    meanRsAffy16(i) = mean(r16.r2(b(a)));

[a, b] = ismember(gtexCluster.syms(gtexCluster.cs1547 == ...
                                       thisID), ...
                      bloodExpGenes);
    meanRs(i) = mean(rg.r2(b(a)));
end

% TODO: the plot here
plotMat = 
h = figure;
heatmap()
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sbloodR2dist_ksdensity_r2DiffGenes_SelectedAll', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


% correlation of new Rs with CT var  % get it from 22. 
[a, b] = ismember(markerListsReO, [Ads.geneSyms', Bds.geneSyms']);
wholesMat = [Ads.samat' Bds.samat']';
wholeMat = [Ads.mat' Bds.mat']';

% correlation of the three of them

% 25. writing supplemental files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 25.1 functional terms for clusters. 
load('~/resultsAndFigures/secondProject/gtexClusters_blood.mat') % gtexCluster
cs = gtexCluster.cs1547;
clusterMemberCount = hist(cs, unique(cs))
inClusterIds = find(clusterMemberCount >= 20);

load('~/resultsAndFigures/secondProject/bloodClusterInfo.mat') %
                                                               % cinfo 
selectedCs = cinfo.outperm([1:54, 57, 59, 61]); % removing clusters

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'bloodCluster_fenrich.mat'])

myFile = sprintf(['~/resultsAndFigures/secondProject/suppFiles/redo/' ...
                  'blood_clusterFunctionTerm_secondFormat_try.txt']);
fid = fopen(myFile, 'w')
fprintf(fid, ['clusterID\tGOID\tGOTerm\tFDR\n'])

for i = 1:length(selectedCs)
    i
    thisCind = selectedCs(i);
    myFInds = find(clusterFE(thisCind).passedFDRIDs);
    
    % writing the file
    for j = 1:length(myFInds)
        kado = cinfo.inIDs(thisCind);
        fprintf(fid, '%d\t%d\t%s\t%e\n', kado, full(clusterFE(thisCind).passedFDRIDs(j)), ...
                clusterFE(thisCind).passedFDRTerms{j},full( clusterFE(thisCind).passedFDRs(j)));
    end
end

% 25.2 gene list in different cell types

load('~/resultsAndFigures/secondProject/bloodDiffGenesInfo.mat') 
diff.cellTypes{:}
% cellTypes we are getting: 
tempSelected = [1, 2, 4:13, 22];
inTypes = diff.cellTypes(tempSelected);

inGenes = diff.genes((sum(diff.mat(:, tempSelected)') > 0));
inGenesMat = diff.mat((sum(diff.mat(:, tempSelected)') > 0), tempSelected);

fileName = ['~/resultsAndFigures/secondProject/suppFiles/redo/' ...
            'GTExBlood_diffGenes327.txt'];
fid = fopen(fileName, 'w')
fprintf(fid, ['geneSymbol\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\' ...
              't%s\t%s\t%s\n'], inTypes{1}, inTypes{2}, inTypes{3}, ...
        inTypes{4},inTypes{5},inTypes{6},inTypes{7},inTypes{8},inTypes{9},inTypes{10},inTypes{11},inTypes{12},inTypes{13})

for i = 1:length(inGenes)
    fprintf(fid, ['%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\' ...
                  't%d\t%d\n'], inGenes{i}, inGenesMat(i, 1),inGenesMat(i, 2),inGenesMat(i, 3),inGenesMat(i, 4),inGenesMat(i, 5),inGenesMat(i, 6),inGenesMat(i, 7),inGenesMat(i, 8),inGenesMat(i, 9),inGenesMat(i, 10),inGenesMat(i, 11),inGenesMat(i, 12),inGenesMat(i, 13));
end

% 25.3 gene cluster labels and R2 that is for all genes in GTEx
% clusters

fileName = ['~/resultsAndFigures/secondProject/suppFiles/redo/' ...
            'GTExBlood_geneClusterAndR2_redo.txt'];
fid = fopen(fileName, 'w')
fprintf(fid, ['geneSymbol\tclusterID\tR2\n'])

[a, b] = ismember(gtexCluster.syms, bloodExpGenes);
grs = rg.r2(b(a));

for i = 1:length(gtexCluster.syms)
    fprintf(fid, '%s\t%d\t%.3f\n', gtexCluster.syms{i}, gtexCluster.cs1547(i), ...
            grs(i));
end


% D1. testing exp level of diff genes
%%%%%%%%%%%%%%% testing exp level of diff genes

thisGeneSet = diffGenes(logical((sum(diffMat(:, 21:22)') - diffMat(:,13)') >0));
thisGeneSet = diffGenes(logical((diffMat(:, 4))));
thisGeneSet = diffGenes(((diffMat(:, 2) - diffMat(:,1))>0));
thisGeneSet = diffGenes(((diffMat(:, 12) - diffMat(:,4))>0));
thisGeneSet = diffGenes(((diffMat(:, 4)' - sum(diffMat(:, [ ...
                    5:12])'))>0));
thisGeneSet = diffGenes((sum(diffMat(:, 4:10)') == 7) - sum(diffMat(:, ...
                                                  [1:3, 11:22])')>0);

[a, b] = ismember(thisGeneSet, Ads.geneSyms);
sum(a)
mybox = Ads.samat(b(a), :);
h = figure
boxplot(mybox)
boxplot(zscore(mybox')')
set(gca, 'XTick', [1:12], 'XTickLabels', groups)
xtickangle(90)

h = figure
hold on
for i = 1:12
    plot(repmat(i, 1, size(mybox, 1)), mybox(:, i), 'o')
    set(gca, 'XTick', [1:12], 'XTickLabels', groups)
    xtickangle(90)
end

h = figure
heatmap([groups], [1:length(mybox)], mybox)

% The T's agree with each other and also share a set with NK. it is
% the same for NK. So we have TCell markers and NKcell markers and
% NKTCell markes. NKactivated, minus tcells, are distinguisehd from
% Tcells, but this is not true for NKresting, their expression is
% also high in other cell types. So for clusters, we should look at
% Tcells, NKactivated filtered, BCells, TCELLSNK. Neutrophils stand
% out, but Monocytes do not, even minus Neutrophils. Plasma minus
% Bcells stands out. So let's look at these sets of markers in the
% clusters and see if we can see them. 
% Bcells agree with each other and distinguish
% NK cells agree 

%% get different gene sets:

% 1. B cells gene set (the overlap)
thisGeneSet = diffGenes(logical(sum(diffMat(:, 1:2)') == 2));
% > also for fold2 diff exp
thisGeneSet = diffGenes(logical(diffMat(:, 1))); % myboxColID = 1
[a, b] = ismember(thisGeneSet, Ads.geneSyms);
sum(a)
mybox = Ads.samat(b(a), :);
minMat = zeros(size(mybox));
myboxColID = 1 % the column in mybox that belongs to this cell type
               % (see the variable <groups>)
for i = 1:12
    minMat(:, i) = mybox(:, myboxColID) - mybox(:, i) ;
end

sib = minMat >= 1;
sum(sib)

% 2. NK activated minus TCells distinguish NK (NKresting doesn't distinguish NK)
thisGeneSet = diffGenes(((diffMat(:, 12)' - sum(diffMat(:, [ ...
                    4:10])'))>0));

% 

% avreage expression level of genes in each cluster: CT of each
% cluster, basically. and how close and far are they, and how
% coherent are they. which clusters are closer to each other, which
% are far away etc etc. you can do the cluster plot for blood, you
% can do CT profiles, enrichment and etc

% gropuing of the cells from the marker mat, and their expression
% in the cell dataset, and their presence in cluster. 

% also, genes that are expressed hihgly in one cell type. 

% do the cluster plot?

% the regression > yes, the regression with these genes. but then
% you have to do the random too. 

% plot the expression level of genes in different clusters - do
% they have markers of different cell types?

% uh, the rest: exp profiles, correlation and regression

%% group the genes in diffMat based on their cell type expression
%% profiles
[a, b] = ismember(diffGenes, Ads.geneSyms);
expProfile = Ads.samat(b(a), :);

sib = corr(expProfile');

% let's cluster them: 


% Draft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


myDS = zeros(size(Ads.mat));
sib = mean(Ads.mat);
for i = 1:size(Ads.mat, 1)
    i
    myDs(i,:) = Ads.mat(i, :) - sib;
end

h = figure
plot(labels, myDs(1, :), 'o')
h = figure
plot(labels, Bds.mat(1, :), 'o')

[a, b] = ismember(data{1}, gpl96.probeID);
geneID = gpl96.geneSymbol(b(a));
geneMap = containers.Map(gpl96.uniqueGeneSymbols, [1:gCount]);
probeMap = containers.Map(data{1}, geneID);
finalMat = zeros(gCount, sampleCount);
divMat = zeros(gCount, sampleCount);

tic
for j = 1:pCount
    sib = values(geneMap, values(probeMap, data{1}(j)));
    finalMat(sib{1}, :) = finalMat(sib{1}, :) + dataExpr(j, :);
    divMat(sib{1}, :) = divMat(sib{1}, :) + 1;
end
toc

finalMat = finalMat./divMat;

dataSet.mat = finalMat;
dataSet.GEOID = GEOID;
dataSet.Gemma = 'null';

clear expr
load(['/home/mfarahbod/data/affyArray/bloodCellTypeDataset/' ...
      'GSE22886_B.mat'])
dmwrite(expr, ['/home/mfarahbod/data/affyArray/' ...
               'bloodCellTypeDataset/GSE22886_B.txt'])



% check out the heatmap with labels and sample correltaions. 

% identification of samples and cell types 

% extract the markers

% load the blood bulk tissue 

% identify the clusters 

% uh, the rest


% 


 