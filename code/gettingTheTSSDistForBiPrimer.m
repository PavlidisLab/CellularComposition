% I first had this code for the links present in more than 10
% networks and after that, genes with more than 2 links. I want to
% revisit it, have it for links present in more than 5 networks out
% of the neuronal sum netwok only. Inh and Exc

% 1. everything. 
% 2. getting the neuronal sum network and cluster it
% 3. just printing the count of links present in multiple networks

% sumNet comes from networksOverlapLR.m
load(['~/resultsAndFigures/secondProject/' ...
      'overlappingLinksBetweenSCNets/sumNets.mat'])
sumNet = sumNets.exc + sumNets.inh;
sumNetNDs = sum(sumNet>0);
sum(sumNetNDs > 0)

upol = sumNet(logical(triu(ones(size(sumNet)), 1)));
book = upol(upol >0);

h = figure
upolHist = hist(book, 1:69);
bar(log10(upolHist+1))
ylabel('log10 count of links', 'FontSize', 16)
xlabel('count of overlaps', 'FontSize', 16)
title(['count of links with overlapping presence between the ' ...
       'networks'], 'FontSize', 16)
figFolder = ['~/resultsAndFigures/secondProject/overlappingLinksBetweenSCNets/']
file = sprintf('%sHistogramOfTheLinkCountsForRepeats', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% load the ensemble file 
% >> getting the NCBI id for my genes 
file = '~/data/brainSingleCell/human_MTG_2018-06-14_genes-rows.csv'
fid = fopen(file)

headers = fgets(fid);

genesInfo = textscan(fid, [repmat('%s', 1, 5)], ...
                'Headerlines', 1, 'Delimiter', ',');

geneSyms = cell(1, length(genesInfo{1}));
for i = 1:length(geneSyms)
    i
    book = genesInfo{1}{i};
    geneSyms{i} = book(2:end-1);
end

load('~/data/brainSingleCell/filDataSet_exon_V4.mat')
scGeneSyms = filDataSet.geneSyms;
[a, b] = ismember(scGeneSyms, geneSyms);
scNCBIIDs = (genesInfo{3}(b));

% >> get the markers for cortex mouse
file = '~/data/genePositionWNCBI.txt'
fid = fopen(file)

headers = fgets(fid);

markerListsRaw = textscan(fid, [repmat('%s', 1, 3) repmat('%s', 1, 5)], ...
                'Headerlines', 1, 'Delimiter', '\t');
markerListsRaw
refGeneSyms = markerListsRaw{3};

[selectedGenes1, tempPosition1] = ismember(scNCBIIDs, ...
                                         markerListsRaw{8});
[selectedGenes2, tempPosition2] = ismember(scGeneSyms, ...
                                         markerListsRaw{3});

% for each gene, we want the postion. 
geneCHRs = zeros(length(selectedGenes1), 22);
dir = zeros(length(selectedGenes1));
tss = zeros(length(selectedGenes1));
for i = 1:length(selectedGenes1)
    if(selectedGenes1(i))
        chstr = markerListsRaw{4}(tempPosition1(i));
        ch = str2num(chstr{1});
        book = markerListsRaw{7}(tempPosition1(i));
        dir(i) = str2num(book{1});
        geneCHRs(i, ch) = 1;
        book = markerListsRaw{5}(tempPosition1(i));
        tss(i) = str2num(book{1});
    else 
        if(selectedGenes2(i))
            chstr = markerListsRaw{4}(tempPosition2(i));
            ch = str2num(chstr{1});
            book = markerListsRaw{7}(tempPosition2(i));
            dir(i) = str2num(book{1});
            geneCHRs(i, ch) = 1;
            book = markerListsRaw{5}(tempPosition2(i));
            tss(i) = str2num(book{1});
        end
    end
end    

halva = geneCHRs * geneCHRs';

[a, b, c] = find(triu(halva, 1));

diffMat = zeros(size(halva));
for i = 1:length(a)
    diffMat(b(i), a(i)) = abs(tss(b(i)) - tss(a(i))) * dir(b(i)) *dir(a(i));
end

fil1 = diffMat > -1000;
fil2 = diffMat < 0;

bid = (fil1 + fil2) == 2;
[a, b, c] = find(bid);

% looking for the genes with two bidirectional partners
nds = sum(bid + bid');
find(nds > 1)
nds(nds>1)

% example: (3 genes) 
kado = bid + bid';
scGeneSyms(logical(kado(466,:)))
scGeneSyms(466)

% example: (4 genes)
scGeneSyms(logical(kado(1946,:)))
scGeneSyms(1946)

% example: (4 genes)
scGeneSyms(logical(kado(15738,:)))
scGeneSyms(15738)

% overlap of bid with binNet05 from the gtex bulk tissue
kado = scGeneSyms(nds >0);
[a, b] = ismember(kado, gtexGeneSyms);
fullBid = bid + bid';
book = fullBid .* binNet05;

i = 1

scGeneSyms(a((i)))
scGeneSyms(b((i)))
i = i + 1

kado = sum(bid);

exp = sumNet .* fullBid;

[a, b, c] = find(exp);
h = figure
plot(c)
[aso, bso, cso] = find(triu(sumNet, 1));

% here, I change the threshold from 10 to 5
inds = find(cso >= 5);

i = 1
as(inds(i))
scGeneSyms(as(inds(i)))
bs(inds(i))
scGeneSyms(bs(inds(i)))
cs(inds(i))
i = i + 1

% repeats in GTEx and other reasons for observed repeated coexpression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% how many of high repeats are also represented in GTEx? it should
% be really really high. 
load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat') 

gtexNet = GTExFiveNets.nets(2).net01;
gtexGeneSyms = ...
    GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(2).expGenes);

whos sumNet
whos scGeneSyms

[a, b] = ismember(scGeneSyms, gtexGeneSyms);
sum(a)

gtexSCMapped = zeros(size(sumNet));
gtexSCMapped(a, a) = gtexNet(b(a), b(a));

% again, change 10 to five here, I want more links
highCorrSumNet = triu(sumNet >= 5);

book = highCorrSumNet .* gtexSCMapped;
sum(book(:))

gtexFlagged = book .* 100;
highCorrSumNetCounts = highCorrSumNet .* sumNet;

sib = gtexFlagged + highCorrSumNetCounts;
[as, bs, cs] = find(sib);

present = (cs > 100) .* 2;
for i = 1:length(as)
    i
    %    [scGeneSyms(as(i)) scGeneSyms(bs(i)) cs(i)]
    if cs(i) < 100
        %        [scGeneSyms(as(i)) scGeneSyms(bs(i)) cs(i)]
        [localA, localB] = ismember([scGeneSyms(as(i)) scGeneSyms(bs(i))], ...
                                    gtexGeneSyms);
        NotBothIn = ~(localA(1) && localA(2));
        present(i) = gtex05(as(i), bs(i)) + -1*(NotBothIn);
    end
end

% present == -1 : gene is absent

i = 1

i = i + 1
[scGeneSyms(as(i)) scGeneSyms(bs(i)) cs(i)]
[localA, localB] = ismember([scGeneSyms(as(i)) scGeneSyms(bs(i))], ...
                            gtexGeneSyms)
present(i)

% find the pair of genes which overlap 
coverFlag = zeros(1, length(present));
for i = 1:length(as)
    book = [scGeneSyms(as(i))];
    temp = strfind(scGeneSyms(bs(i)), [book{1} '-']);
    sf1 = temp{1};
    book = [scGeneSyms(bs(i))];
    temp = strfind(scGeneSyms(as(i)), [book{1} '-']);
    sf2 = temp{1};
    if (length([sf1 sf2]))
        coverFlag(i) = 1;
    end
end

sib = present(logical(coverFlag));

familyFlag = zeros(1, length(present));
for i = 1:length(as)
    
    n = 4;
    if (length(scGeneSyms{as(i)}) < 4) || (length(scGeneSyms{bs(i)}) ...
                                           < 4)
        n = 3;
    end
    
    book = [scGeneSyms(as(i))];
    book = book{1};
    temp = strfind(scGeneSyms(bs(i)), [book(1:n)]);
    sf1 = temp{1};
    book = [scGeneSyms(bs(i))];
    book = book{1};
    temp = strfind(scGeneSyms(as(i)), [book(1:n)]);
    sf2 = temp{1};
    if (length([sf1 sf2]) == 2)
        familyFlag(i) = 1;
    end
end

sib = present(logical(familyFlag));

% if we remove the familyFlag and coverFlag links, what is the
% remanining

whos aso
whos bso

filas = as;
filbs = bs;

filas = aso(inds);
filbs = bso(inds);

sib = familyFlag + coverFlag;

filas = as(sib == 0);
filbs = bs(sib == 0);

length(unique([filas, filbs]))
highConservedNet = sparse(filas, filbs, ones(1,length(filbs)), 16789, ...
                          16789);

selectedSumNet = sparse(aso(inds), bso(inds), cso(inds), 16789, ...
                          16789);

nds = sum(highConservedNet + highConservedNet');

whos highConservedNet
ndSelectedNet = highConservedNet(nds>2, nds> 2);
ndSelectedSumNet = selectedSumNet(nds > 2, nds>2);
ndSelectedGeneSyms = scGeneSyms(nds>2);

ndSelectedNDs = nds(nds>2);

[sorA, sorB] = sort(ndSelectedNDs, 'descend');
ndSelectedGeneSymsSorted = ndSelectedGeneSyms(sorB);
ndSortedNet = ndSelectedNet(sorB, sorB);
ndSortedSumNet = ndSelectedSumNet(sorB, sorB);

h = figure
plotMat = ndSortedNet + ndSortedNet';
plotMat = ndSelectedSumNet + ndSelectedSumNet';
heatmap(ndSelectedGeneSymsSorted(1:200), ndSelectedGeneSymsSorted(1:200), ...
        plotMat(1:200, 1:200))


[cSortA, cSortB] = sort(c); % from clusterMyGenes.m, take
                            % ndSortedNet there and bring back c! I
                            % also got the outperm out from the
                            % dendrogram and going to use that as
                            % the order
h = figure
heatmap(ndSelectedGeneSymsSorted(outperm(plotInds)), ndSelectedGeneSymsSorted(outperm(plotInds)), ...
        plotMat(outperm(plotInds), outperm(plotInds)), 'GridVisible', 'off')
colormap(1-gray)

plotMat = ndSortedSumNet + ndSortedSumNet';
h = figure
plotInds = [1:100, 800:900]
heatmap(ndSelectedGeneSymsSorted(outperm(plotInds)), ndSelectedGeneSymsSorted(outperm(plotInds)), ...
        log2(plotMat(outperm(plotInds), outperm(plotInds))), 'GridVisible', 'off')
colormap(1-(gray*80/100))


figFolder = ['~/resultsAndFigures/secondProject/overlappingLinksBetweenSCNets/']
file = sprintf('%s', figFolder, i);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

h = figure
heatmap(ndSelectedGeneSymsSorted(cSortB), ndSelectedGeneSymsSorted(cSortB), ...
        plotMat(cSortB, cSortB), 'GridVisible', 'off')

result.SCSumBinaryNet = plotMat(cSortB, cSortB);

h = figure
heatmap(log2(sumNet(tampb, tampb)+1), 'GridVisible', 'off')
% clearly two clusters with different hubs 
result.sumNet = sumNet(tampb, tampb);

% now the same plot from GTEx network
[tempa, tampb] = ismember(ndSelectedGeneSymsSorted(cSortB), ...
                          scGeneSyms);
h = figure
fullBinNet05 = binNet05SC + binNet05SC';
heatmap(fullBinNet05(tampb, tampb)+0, 'GridVisible', 'off')

h = figure
fullBinNet01 = binNet01SC + binNet01SC';
heatmap(fullBinNet01(tampb, tampb)+0, 'GridVisible', 'off')
result.gtexBinNet01 = fullBinNet01(tampb, tampb);
result.geneSyms = ndSelectedGeneSymsSorted(cSortB);
result.clusters = cSortA;

myGeneList = ndSelectedGeneSymsSorted(cSortB);
geneList = myGeneList(1:18)
geneList = myGeneList(67:352);

save(['~/resultsAndFigures/secondProject/' ...
      'overlappingLinksBetweenSCNets/highRepLinks.mat'], 'result')
% I later save this file with the datastructure overlapRes instead
% of result

load(['~/resultsAndFigures/secondProject/' ...
      'overlappingLinksBetweenSCNets/highRepLinks.mat'])

[a, b] = ismember(result.geneSyms, scGeneSyms);
h = figure
heatmap(log2(inhSumNet(b, b)+1), 'GridVisible', 'off');

h = figure
heatmap(log2(excSumNet(b, b) + 1), 'GridVisible', 'off');


% I should document and plot the whole group 

% load('~/data/GTEx/Brain_Cortex_expGenes.mat')
[a, b] = ismember(scGeneSyms, gtexGeneSyms);
binNet01 = GTExFiveNets.nets(2).net01;
binNet01SC = zeros(16789);
binNet01SC(a, a) = binNet01(b(a), b(a));

% building the more dens network
smallMat = dataSet.mat(b(a), :);
whos smallMat
sib = corr(smallMat');

upCorr = sib(logical(triu(ones(size(sib)), 1)));
qs = quantile(upCorr, [.9, .95 ,.99, .995]);

binNet05 = triu(sib  > qs(2), 1);

binNet05SC = zeros(16789);

binNet05SC(a, a) = binNet05;

highCountGtexSCMapped(a, a) = binNet05;
gtex05 = highCountGtexSCMapped + highCountGtexSCMapped';
% book = highCorrSumNet .* highCountGtexSCMapped;

% see how they are distributed between clusters in GTEx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a, b] = ismember(result.geneSyms, gtexCluster.syms);

book = gtexCluster.cs1547(b(a));
inClusters = unique(book);
h = figure
myBar = hist(book, inClusters);

[myBar(myBar > 20); inClusters(myBar > 20)'; meanR2(highCountCs)]

h = figure
heatmap(log2(result.sumNet+1), 'GridVisible', 'off')

h = figure
heatmap(result.gtexBinNet01, 'GridVisible', 'off')

% TODO: sort the genes based on the GTEx cluster 
% getting the clusterIDs with >20 genes
highCountCs = (inClusters(myBar>20))

smallSumNet = result.sumNet(a, a);
smallGTExNet = result.gtexBinNet01(a, a);

[a2, b2] = ismember(book, highCountCs);

smallSumNet2 = smallSumNet(a2, a2);
smallGTExNet2 = smallGTExNet(a2, a2);
h = figure
heatmap(smallSumNet2, 'GridVisible', 'off')

h = figure
heatmap(smallGTExNet2, 'GridVisible', 'off')
title('smallGTExNet')

% sorting based on GTExClusters
[as, bs] = sort(book(a2));

h = figure
heatmap(smallSumNet2(bs, bs), 'GridVisible', 'off')

h = figure
heatmap(smallGTExNet2(bs, bs), 'GridVisible', 'off')

% the effect of R2 in the clusters: are some of these genes
% hevaeliy affected by R2 and therefore are part of other cell type
% clusters? 

% clusters 9 and 10 are the most dense, what are the GTEx clusters
% in them
dummy = zeros(1, length(a));
dummy(a) = gtexCluster.cs1547(b(a));

[ad, bd] = ismember(dummy, highCountCs);
dummy(~ad) = 0;
geneClusterLabelsHC = dummy;

sub9 = geneClusterLabelsHC(result.clusters == 9);
sub10 = geneClusterLabelsHC(result.clusters == 10);

bar9 = hist(sub9, unique(geneClusterLabelsHC))
bar10 = hist(sub10, unique(geneClusterLabelsHC))

% >>>>>>>>>>>
load('~/data/GTEx/GTExGeneIDfier') %'gtexGenes'
load('~/data/brainSingleCell/scGeneIDfier.mat') %sc
%<<<<<<<<<<<

% the clusters are a well mix of group of genes from different
% clusters 

[myBar(myBar > 20); inClusters(myBar > 20)'; meanR2(highCountCs)]

% let's do an ensmebl geneID go and back for each of the clusters's
% gene lists

cluster1Genes = result.geneSyms(result.clusters == 1);
cluster9Genes = result.geneSyms(result.clusters == 9);
cluster10Genes = result.geneSyms(result.clusters == 10);

[a, b] = ismember(cluster1Genes, sc.geneSyms);
enID = sc.ensemble(b(a));
for i = 1: length(enID)
    if length(enID{i}) < 3
        enID{i} = 'null'
    end
end

[a, b] = ismember(enID, gtexGenes.ensemble);
geneList = gtexGenes.symbols(b(a));

h = figure
snet = result.sumNet(1:18, 1:18)
heatmap(geneList, geneList, snet(a, a))

% checking the network in string - I also add the gtex symbols to
% the genes here 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
load(['~/resultsAndFigures/secondProject/' ...
      'overlappingLinksBetweenSCNets/highRepLinks.mat'])
load('~/networks/stringExperimental.mat')

% >>>>>>>>>>>
load('~/data/GTEx/GTExGeneIDfier.mat') %'gtexGenes'
load('~/data/brainSingleCell/scGeneIDfier.mat') %sc
% get the ensemble IDs
[a, b] = ismember(result.geneSyms, sc.geneSyms);
enID = sc.ensemble(b(a));
for i = 1: length(enID)
    if length(enID{i}) < 3
        enID{i} = 'null';
    end
end

gtexEns = gtexGenes.ensemble;
for i = 1: length(stringEns)
    if length(gtexEns{i}) < 3
        i
        gtexEns{i} = 'gnull';
    end
end

[a, b] = ismember(enID, gtexEns);
gtexGeneSyms = cell(1, length(a));
gtexGeneSyms(a) = gtexGenes.symbols(b(a));
for i = 1: length(gtexGeneSyms)
    if length(gtexGeneSyms{i}) < 3
        i
        gtexGeneSyms{i} = 'g-null';
    end
end
result.gtexGeneSyms = gtexGeneSyms;
% also refresh the gtex network:
load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat') 
thisNet = GTExFiveNets.nets(2).net005;
thisNetSyms = ...
    GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(2).expGenes);
gtexOverlapNet = zeros(530, 530);
[a, b] = ismember(gtexGeneSyms, thisNetSyms);
gtexOverlapNet(a, a) = thisNet(b(a), b(a));
fullg = gtexOverlapNet + gtexOverlapNet';
result.gtexBinNet005 = fullg;
result.gtexBinNet01 = fullg;
save(['~/resultsAndFigures/secondProject/' ...
      'overlappingLinksBetweenSCNets/highRepLinks.mat'], 'result')
%<<<<<<<<<<<

% get the ensemble IDs
[a, b] = ismember(result.geneSyms, sc.geneSyms);
enID = sc.ensemble(b(a));
for i = 1: length(enID)
    if length(enID{i}) < 3
        enID{i} = 'null';
    end
end

stringEns = string.ENS;
for i = 1: length(stringEns)
    if length(stringEns{i}) < 3
        i
        stringEns{i} = 'snull';
    end
end

[a, b] = ismember(enID, stringEns);

smappedNet = zeros(530, 530);
smappedNet(a, a) = string.net(b(a), b(a));

h = figure
heatmap(log10(smappedNet+1), 'GridVisible', 'off')

result.stringExperimental = smappedNet;
resSC.r2(1:10)
save(['~/resultsAndFigures/secondProject/' ...
      'overlappingLinksBetweenSCNets/highRepLinks.mat'], 'result')

cc = result.clusters == 9;
geneList = result.gtexGeneSyms(cc);

[passedFDRIDs, passedFDRTerms, inGenes, enrichMat] = ...
    fenrich_function(geneList, .1);
% filter for inGenes

whos enrichMat
book = sum(enrichMat');
whos book

plotMat = enrichMat(book > 0, :);
plotGenes = inGenes(book>0);
plotFs = passedFDRTerms;
h = figure
heatmap(plotGenes, plotFs, plotMat)


h = figure
h = figure('units', 'centimeters', 'position', [0,0, 35, 25])
for i = 9:10
    inGenes = find(result.clusters == i)
    
    heatmap(result.geneSyms(inGenes), result.geneSyms(inGenes), ...
            log2(result.sumNet(inGenes, inGenes)+1), 'FontSize', 6, 'GridVisible', ...
            'off')
    figFolder = ['~/resultsAndFigures/secondProject/overlappingLinksBetweenSCNets/']
    file = sprintf('%ssumNet530Genes_%d', figFolder, i);
    set(h, 'PaperOrientation', 'landscape')
    print(h, '-deps', [file '.eps'])
    print(h, '-dpdf', [file '.pdf'])
    saveas(h, [file '.eps'], 'epsc')
end

h = figure
heatmap(result.gtexBinNet005, 'GridVisible', 'off')
figFolder = ['~/resultsAndFigures/secondProject/overlappingLinksBetweenSCNets/']
file = sprintf('%sGTEx005Net530Genes_%d', figFolder, i);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 2. getting the neuronal sum network and cluster it
% get the sumNet from the networksOverlapLR.m

whos sumNet

[a, b, c] = find(triu(sumNet));
inds = c >=5;

myNet = sparse(a(inds), b(inds), c(inds), 16789, 16789);

nds = sum(myNet>0) + sum(myNet'>0);

ndSelect = nds > 2;
ndSelectNet = myNet(ndSelect, ndSelect) >= 5;
ndSelectSumNet = myNet(ndSelect, ndSelect);
gsSelect = scGeneSyms(ndSelect);

% the c and the outperm come from clusterMyGenes, ndSelectNet goes
% to it
h = figure
plotMat = ndSelectSumNet + ndSelectSumNet';
heatmap(log2(plotMat(outperm, outperm)), 'GridVisible', ...
        'off')
colormap(1-gray)
figFolder = ['~/resultsAndFigures/secondProject/overlappingLinksBetweenSCNets/']
file = ...
    sprintf('%sneuronalSumNetClusters_1720genes_GE5RepLinks_net01', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

objectCount = hist(c, unique(c));
clusterIDs = find(objectCount >= 10);

[a, b] = ismember(c, clusterIDs);
finalCs = zeros(size(c));
finalCs(a) = c(a);

[sortC, orderC] = sort(finalCs);
medGeneSyms = gsSelect(orderC(sortC>0));
medGeneSyms(195:234) = [];
finalGeneSyms = medGeneSyms;

medOrder = orderC(sortC>0);
medOrder(195:234) = [];
finalOrder = medOrder;

h = figure
plotMat = ndSelectSumNet + ndSelectSumNet';
heatmap(log2(plotMat(finalOrder, finalOrder)), 'GridVisible', ...
        'off')
colormap(1-gray)
figFolder = ['~/resultsAndFigures/secondProject/overlappingLinksBetweenSCNets/']
file = ...
    sprintf('%sneuronalSumNetClusters_1720genes_GE5RepLinks_net03finalOrder', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% printing the sum network
% sumNet comes from networksOverlapLR.m
load('~/resultsAndFigures/secondProject/overlappingLinksBetweenSCNets/sumNets.mat')
sumNetNDs = sum(sumNet>0);
sum(sumNetNDs > 0)

[inha, inhb, inhc] = find(scSumNets.inh);
[exca, excb, excc] = find(scSumNets.exc);
[botha, bothb, bothc] = find(scSumNets.both);


