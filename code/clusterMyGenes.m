
% get the network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% getting the genes expressed in SC
n = 1
filDataSet
% myNet = scBinNets.nets(localInd).net01;
% thisNetSyms = ...
%     scBinNets.geneSyms(scBinNets.nets(localInd).expGenes);
load(sprintf(['~/networks/allenBrainSC/allFiveNets/' ...
              'binNets_Exon_V4_allFive_net%d.mat'], n))
myNet = net.net005;
thisNetSyms = filDataSet.geneSyms(net.expGenes);

[a, b] = ismember(geneList, thisNetSyms);
sum(a)
smallNet = myNet(b(a), b(a));

% >>> GTEx 
t = 2
myNet = GTExFiveNets.nets(t).net005;
thisNetSyms = ...
    GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(t).expGenes);

[a, b] = ismember(geneList, thisNetSyms);
[a, b] = ismember(myPlotSyms, thisNetSyms);
smallNet = myNet(b(a), b(a)); 
sum(smallNet(:))

smallMat = smallBrainTSNNDBased;

% get the genes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
localGeneList = geneList(a);

smallNet = book;
smallNet = ndSortedNet;
smallNet = ndSelectedNet;
smallNet = ndSortedSumNet;
smallNet = ndSelectedSumNet;
smallNet = overlapRes.inhSum;
smallNet = overlapRes.excSum;
smallNet = ndSelectNet;

% throw away the ones that are not correlated at all
sfm = smallNet + smallNet';
tic
overlap = sfm * sfm;
toc

normMatUp = zeros(size(overlap));
maxOs = diag(overlap);
tic
for i = 1:size(overlap, 1);
    i
    for j = i:size(overlap, 1);
        a = maxOs(i);
        b = maxOs(j);
        tempMin = min(a, b);
        normMatUp(i, j) = tempMin;
        % gtexNormMat(j, i) = tempMin;
    end
end
toc
normMat = triu(normMatUp, 1) + normMatUp' + 1 - sfm;

myTop = overlap./normMat;

myDist = 1 - myTop;
myVDist = squareform(tril(myDist, -1));

z = linkage(myVDist, 'average');
h = figure
[H, T, outperm] = dendrogram(z, size(normMat, 1));
figFolder = ['~/resultsAndFigures/secondProject/overlappingLinksBetweenSCNets/']
file = sprintf('%sneuronalSumNetClusters_1720genes_GE5RepLinks', figFolder, i);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

c = cluster(z, 'cutoff', .5, 'criterion', 'distance');
c = cluster(z, 'cutoff', 1.1547);
max(c)
% selecting the cluster
m = max(c);
h = figure
objectCount = hist(c, unique(c));
hist(c, unique(c)) 

highCounts = find(objectCount >= 15);
labels = unique(c);

[a, b] = ismember(c, highCounts);

cmod = zeros(size(c));
cmod(a) = c(a);

[a, b] = sort(cmod);
orderedsfm = sfm(b(116:end), b(116:end));
orderedsfm = orderedsfm - eye(size(orderedsfm));
orderedsfm(58, 59) = .4;
orderedsfm(59, 58) = .4;
h = figure
heatmap(orderedsfm)
colormap(jet)
whos s2ExpMat

highCounts
se1 = s2ExpMat(c == 69, :);
se2 = s2ExpMat(c == 71, :);
se3 = s2ExpMat(c == 76, :);

varse1 = zeros(size(se1));
for i = 1:size(se1, 1)
    
end

cgcs = zeros(1, length(highCounts));
clcs = zeros(1, length(highCounts));
localInGenes = zeros(size(localGeneList));
for i = 1:length(highCounts)
    tempList = c == labels(highCounts(i));
    localInGenes(tempList) = 1;
    subNet = smallNet(tempList, tempList);
    cgcs(i) = sum(tempList);
    clcs(i) = sum(subNet(:));
end
localInGenes = logical(localInGenes);
densities = clcs ./ (cgcs .* (cgcs-1) ./2)

% now give me the genes which contribtue to these dense clusters
localClustersGenes = localGeneList(logical(localInGenes));

book = smallNet(localInGenes, localInGenes);
sum(book(:))

% find the count of links between these genes in all the SC
% datasets

subsclc = zeros(1:69);
subscgc = zeros(1:69);
for n = 1:69
    n
    load(sprintf(['~/networks/allenBrainSC/allFiveNets/' ...
                  'binNets_Exon_V4_allFive_net%d.mat'], n))
    myNet = net.net005;
    thisNetSyms = filDataSet.geneSyms(net.expGenes);
    
    [a, b] = ismember(thisNetSyms, localClustersGenes);
    smallNet = myNet(a, a);
    
    subsclc(n) = sum(sum(smallNet));
    subscgc(n) = sum(a);
end

h = figure
heatmap(1, netLabels(7:end), (subsclc ./ sclc)')
colormap(jet)
h = figure
heatmap(1, netLabels(7:end), (sclc)')


% now back to the original GTExbrain. 

t = 5
myNet = GTExFiveNets.nets(t).net005;
thisNetSyms = ...
    GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(t).expGenes);

[a, b] = ismember(localClustersGenes, thisNetSyms);
smallNet = myNet(b(a), b(a)); 
sum(smallNet(:))








