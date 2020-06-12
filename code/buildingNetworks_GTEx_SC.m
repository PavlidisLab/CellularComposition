% I build the binary networks for GTEx and allen brain single cell
% data here.

% >> GTEx networks
clear 
networkFolder = '~/networks/GTEx/binNets/'

GTExDSInds = [54, 13, 36, 37, 39]; 

% load datasets
load(['~/data/GTEx/' ...
      'GTExAllDataSetFromGeneLevel_v6p_newWithBlood_RPM.mat'])

uniqueGenes = unique(rpmGTExDS.genes);
myConMap = containers.Map(rpmGTExDS.genes, 1: ...
                          length(rpmGTExDS.genes));

GTExFiveNets.uniqueGeneSyms = uniqueGenes;
GTExFiveNets.GeneSyms = rpmGTExDS.genes;
GTExFiveNets.transCriptIDs = rpmGTExDS.Trans;
qs = cell(1, 5);
for i = 1:5
    i
    myDS = rpmGTExDS.dataSets(GTExDSInds(i));
    net(i).tissue = rpmGTExDS.dataSets(GTExDSInds(i)).tissue;
    
    % collapsing the transcripts with the same symbol
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
    net(i).expGenes = expGenes;
    
    smallMat = myMat(expGenes, :);
    [sib p] = corr(log2(smallMat+1)');
    upCorr = sib(logical(triu(ones(size(sib)), 1)));
        qs = quantile(upCorr, [.99, .995]);
    %qs{i} = quantile(upCorr, 1000);
    
    binNet01 = triu(sib  > qs(1), 1);
    net(i).net01 = sparse(binNet01);
    
    binNet005 = triu(sib  > qs(2), 1);
    net(i).net005 = sparse(binNet005);
end
GTExFiveNets.nets = net;

save('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat', ...
     'GTExFiveNets')
load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat')

% >> brain SC
clear 

% loading the data - exon + intron
dataFolder = ['~/data/brainSingleCell/']
load('~/data/brainSingleCell/filDataSet_intronAndExon_V4.mat')
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exonAndIntron_V4.mat'])
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exonAndIntron_V4_clusterLabels.mat'])


% loading data - exon
load('~/data/brainSingleCell/filDataSet_exon_V4.mat')
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exon_V4.mat'])
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exon_V4_clusterLabels.mat'])

expMat = filDataSet.expMat(:, clusterMeta.inCells);
% normalizing for the million count
% sumExp = sum(expMat);
% milFac = sumExp ./ 1000000;
% normExp = zeros(size(expMat));
% sampleCount = size(expMat, 2);
% for i = 1:sampleCount
%     normExp(:, i) = expMat(:, i)./milFac(i);
% end

scBinNets.geneSyms = filDataSet.geneSyms;
for i = 15:length(clusterMeta.sortedClusterNames) % for each cluster
    i
    net(i).clusterName = clusterMeta.sortedClusterNames(i);
    [a, b] = ismember(clusterMeta.clusters, ...
                      clusterMeta.sortedClusterNames(i));
    %    thisExpMat = log2(expMat(:, a)+1);
    thisExpMat = log2(expMat(:, a)+1);    
    thisSCount = sum(a);
    thisExpGenes = logical(sum(thisExpMat>0, 2) > thisSCount/2);
    net(i).expGenes = thisExpGenes; 
        
    % get the corr and the networks
    smallMat = thisExpMat(thisExpGenes, :);
    sib = corr(smallMat');

    upCorr = sib(logical(triu(ones(size(sib)), 1)));
    qs = quantile(upCorr, [.99, .995]);
    
    binNet01 = triu(sib  > qs(1), 1);
    net(i).net01 = sparse(binNet01);
    
    binNet005 = triu(sib  > qs(2), 1);
    net(i).net005 = sparse(binNet005);        
end

scBinNets.nets = net;
save('~/networks/allenBrainSC/binNets_exon_V4.mat', ...
     'scBinNets')

save('~/networks/allenBrainSC/binNets_exonAndIntron_V4_cpmNorm.mat', ...
     'scBinNets')

whos expMat

expMat = filDataSet.expMat(:, clusterMeta.inCells);
% normalizing for the million count
sumExp = sum(expMat);
milFac = sumExp ./ 1000000;
normExp = zeros(size(expMat));
sampleCount = size(expMat, 2);
for i = 1:sampleCount
    normExp(:, i) = expMat(:, i)./milFac(i);
end

expMat = normExp;

scBinNets.geneSyms = filDataSet.geneSyms;
expGeneCounts = zeros(1, 76);
fdrFacRatio = zeros(1, 76);
for i = 1:length(clusterMeta.sortedClusterNames) % for each cluster
    i
    net.clusterName = clusterMeta.sortedClusterNames(i);
    [a, b] = ismember(clusterMeta.clusters, ...
                      clusterMeta.sortedClusterNames(i));
    %    thisExpMat = log2(expMat(:, a)+1);
    thisExpMat = log2(expMat(:, a)+1);    
    thisSCount = sum(a);
    thisExpGenes = logical(sum(thisExpMat>0, 2) > thisSCount/5);
    net.expGenes = thisExpGenes; 
    expGeneCounts(i) = sum(thisExpGenes);
    
    % get the corr and the networks
    smallMat = thisExpMat(thisExpGenes, :);
    
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
    [a, b] = sort(upP1, 'ascend');
    qs = quantile(a(1:fdrFac), 1000);

    count = fdrFac;
    myvalue = 0;
    k = 0
    while myvalue < .1
        k;
        myvalue = a(k + 1) * (count - k);
        k = k + 1;
    end
    
    BHBinNet01 = triu(p1 <= (a(k)), 1);

    myvalue = 0;
    k = 0
    while myvalue < .2
        k;
        myvalue = a(k + 1) * (count - k);
        k = k + 1;
    end
    
    BHBinNet02 = triu(p1 <= (a(k)), 1);

    % getting the network with certain density
    q = quantile(a(1:fdrFac), [.005, 0.01]);
    binNet005 = p1 < q(1);
    binNet01 = p1 < q(2);
    
    net.expGenes = thisExpGenes;
    net.net005 = binNet005;
    net.net01 = binNet01;
    net.BH01 = BHBinNet01;
    net.BH02 = BHBinNet02;
    net.Bon = BonBinNet;
    save(sprintf('~/networks/allenBrainSC/rpm_binNets_Exon_V4_allFive_net%d.mat', ...
                 i), 'net') 
end

h = figure('units', 'centimeters', 'position', [0,0, 20, 10])
plot(fdrFacRatio(1:75), 'o')
hold on
plot(fdrFacRatio(1:75))
ylim([-.1, 1.1])
title(['ratio of the gene-pairs jointly expressed in >= 20% of ' ...
       'samples'])
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sratioOfLinks', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


h = figure('units', 'centimeters', 'position', [0,0, 20, 10])
plot(expGeneCounts(1:75), 'o')
hold on
plot(expGeneCounts(1:75))
title('Count of genes with exp-value > 0 in >=20% of samples')
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sexpressedGeneCounts', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


scBinNets.nets = net;
save('~/networks/allenBrainSC/binNets_exon_V4_fdr.mat', ...
     'scBinNets')

save('~/networks/allenBrainSC/binNets_exonAndIntron_V4_cpmNorm.mat', ...
     'scBinNets')

% checking the reproducibility of the networks by shuffling the
% data : Given a dataset, given a count of genes, give me the count
% of links between a randomly selected genes in the network for
% their shuffled expression level. Since the expression level
% itself is a coexpression partner, do the shuffle for all. 

%%% This is the shuffling test
% getting the genes in the cluster

expMat = filDataSet.expMat(:, clusterMeta.inCells);
for gc = [200 210 220]
    for i = 1:69 % the SC cluster 
        i

        load(sprintf('~/networks/allenBrainSC/pQuantiles/binNets_Exon_V4_quantiles_net%d.mat', ...
                     i))
        thr01 = qs(10);
        thr005 = qs(5);

        net.clusterName = clusterMeta.sortedClusterNames(i);
        [a, b] = ismember(clusterMeta.clusters, ...
                          clusterMeta.sortedClusterNames(i));
        %    thisExpMat = log2(expMat(:, a)+1);
        thisExpMat = log2(expMat(:, a)+1);    
        thisSCount = sum(a);
        thisExpGenes = logical(sum(thisExpMat>0, 2) > thisSCount/5);

        selectedExpMat = thisExpMat(thisExpGenes, :);
        linkCounts01 = zeros(1000, 1);
        linkCounts005 = zeros(1000, 1);
        
        % I should get the same count of genes, shuffle their
        % expression each time and do the thing
        for j = 1:1000
            j
            rs = datasample(1:sum(thisExpGenes), gcd);
            smallMat = selectedExpMat(rs, :);
            smallMatNAN = smallMat;
            smallMatNAN(smallMatNAN == 0) = nan;
            tic
            [sib1, p1] = corr(smallMatNAN', 'rows', ...
                              'pairwise', 'Tail', 'right');
            toc

            kado = (smallMat > 0) + 0;
            book = kado * kado';

            throwOff = book < 20;
            p1(throwOff) = 1;
            upP1 = p1(logical(triu(ones(size(p1)), 1)));
            linkCounts01(j) = sum(upP1 < thr01);
            linkCounts005(j) = sum(upP1 < thr005);
        end
        linkCounts.thr01 = linkCounts01;
        linkCounts.thr005 = linkCounts005;
        save(sprintf(['~/networks/allenBrainSC/randomNetLC_perCluster/GTExCluster_1574/c%d/' ...
                      'clusterLinkCount_cellType%d_gc%d.mat'], gc, i, gc), 'linkCounts')
    end
end

%% This is the network test (in each given network, what are the odds)

load(['~/resultsAndFigures/secondProject/gtexClusterRep/' ...
      'scClusterRep.mat'])

for i = 1:69 % the SC cluster 
    i
    
    % load the net
    load(sprintf('~/networks/allenBrainSC/allFiveNets/binNets_Exon_V4_allFive_net%d.mat', ...
                 i))

    myNet01 = net.net01;
    myNet005 = net.net005;

    % for each count of test genes, get the sample from the
    % correlation matrix 
    geneCounts = unique(scGTExClusterRep.gc(i, :));
    linkCounts01 = zeros(length(geneCounts), 10000);
    linkCounts005 = zeros(length(geneCounts), 10000);
    for k = 1:length(geneCounts)
        k
        gc = geneCounts(k);
        if gc >= 5
            tic
            for j = 1:10000
                rs = datasample(1:size(net.net01, 1), gc);

                smallMat = myNet01(rs, rs);
                linkCounts01(k, j) = sum(smallMat(:));

                smallMat = myNet005(rs, rs);
                linkCounts005(k, j) = sum(smallMat(:));
            end
            toc
        else 
            linkCounts01(k, :) = -1;
            linkCounts005(k, :) = -1;
        end
    end
    linkCounts.thr01 = linkCounts01;
    linkCounts.thr005 = linkCounts005;
    linkCounts.geneCounts = geneCounts;
    save(sprintf(['~/networks/allenBrainSC/randomNetLC_perCluster/gtexCluster_1547/' ...
                  'rpm_clusterLinkCount_cellType%d_10k.mat'], i), 'linkCounts')
end

% over the SC datasets

% over the clusters 

% getting the count of links and network densityes for raw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear net
lcs = zeros(1, 69);
dns = zeros(1, 69);
for i = 39:69
    i
    load(sprintf('~/networks/allenBrainSC/allFiveNets/binNets_Exon_V4_allFive_net%d.mat', ...
                 i)) 
    lcs(i) = sum(sum(net.BH01));
    gc = size(net.BH01, 1);
    dns(i) = lcs(i)/ (gc * (gc-1)/2);
end

min(dns)
h = figure 
plot(log10(dns))
hold on
plot(log10(dns), 'o')
xlabel('networks')
ylabel('log10 Density')
set(gca, 'FontSize', 15)
figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%sSCnetDensities', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')




