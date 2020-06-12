% the networks I have are the SC, the GTEx brain, liver, blood, the
% TAN blood, TAN liver, TAN brain, the CTC and the 10 simulated. I have to write a code to run
%  over them and load them pairwise. 69, 3 3 1 10: 86 networks

% 1. just getting the overlap between the networks 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 

% SC 
load('~/data/brainSingleCell/filDataSet_exon_V4.mat')
scGeneSyms = filDataSet.geneSyms;
clear filDataSet;

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

netLikeRatiooverlap = zeros(86, 86);
netLinkCountoverlap = zeros(86, 86);
netGeneOverlap = zeros(86, 86);
netLinkCounts = zeros(86, 86);

for i = 1:86
    'First'
    i
    % >>>>>>> first net SC nets
    if i <= 69
        n = i % to 69
        load(sprintf(['~/networks/allenBrainSC/rpmAllFiveNets/' ...
                      'rpm_binNets_Exon_V4_allFive_net%d.mat'], n))
        % scNet = net.net005;
        % thisNetSyms = scGeneSyms(net.expGenes);
        firstNet = triu(net.net005, 1);
        firstSyms = scGeneSyms(net.expGenes);
    end

    % >>>>>>>> first net SC bulk Sim
    if (i > 69) && (i <= 74)
        n = i - 69 % to 10
        load(sprintf('~/resultsAndFigures/secondProject/SimBulkNetworksFromSC/withCounts/bulkFromSC_%d_newComb3_withCounts.mat', ...
                     n)) % bulkFromSC
                         % thisNetSyms = scGeneSyms;
                         % myNet = bulkFromSC.binNet005;
        firstNet = triu(bulkFromSC.binNet005, 1);
        firstSyms = scGeneSyms;
    end
    
        % >>>>>>>> first net SC bulk Sim
    if (i > 74) && (i <= 79)
        n = i - 69 % to 10
        load(sprintf('~/resultsAndFigures/secondProject/SimBulkNetworksFromSC/bulkFromSC_%d_newComb7_noVar.mat', ...
                     n)) % bulkFromSC
                         % thisNetSyms = scGeneSyms;
                         % myNet = bulkFromSC.binNet005;
        firstNet = triu(bulkFromSC.binNet005, 1);
        firstSyms = scGeneSyms;
    end

    
    % >>>>>>>>>> TAN
    % I get the symbols for mapping, however the density comes
    % directly from the density of the mapped network, for affy it
    % is different
    if (i > 79) && (i <= 82) 
        t = i - 79  % 1 2 3 

        tissue = tissues{t}
        load( ['~/networks/tissues/' tissue '/' ...
               'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
        thisExpGenes = logical(affyExpMat(:, t));
        firstNet = binNet(thisExpGenes, thisExpGenes);
        firstSyms = affyGeneSyms(thisExpGenes);
        clear thisExpGenes
    end

    % >>>>>>>> GTEx 
    if (i > 82) && (i <= 85) % I should fix affy density
        t = i - 82 % 1 2 3
        
        load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat') 
        firstNet = GTExFiveNets.nets(t).net005;
        firstSyms = ...
            GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(t).expGenes);
    end

    % CTC
    if i == 86
        % load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
        %       'correctedBinNets_logCorrected_jusResiduals_scBasedMarkers' ...
        %       '.mat'])

        load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
              'correctedBinNets_logCorrected_jusResiduals_scBasedMarkers_redo.mat'])

        firstNet = ctc.net005;
        firstSyms = ctc.geneSyms;
    end
    
    % >>>>>>>>>>>>>>>>>>>>>>>>> inner loop starts
    % >>>>>>>>>>>>>>>>>>>>>>>>>
    for j = (i+1):86
        j
        % >>>>>>> first net SC nets
        if j <= 69
            n = j % to 69
            load(sprintf(['~/networks/allenBrainSC/rpmAllFiveNets/' ...
                          'rpm_binNets_Exon_V4_allFive_net%d.mat'], n))
            % scNet = net.net005;
            % thisNetSyms = scGeneSyms(net.expGenes);
            secondNet = triu(net.net005, 1);
            secondSyms = scGeneSyms(net.expGenes);
        end

        % >>>>>>>> first net SC bulk Sim
            if (i > 69) && (i <= 74)
        n = i - 69 % to 10
        load(sprintf('~/resultsAndFigures/secondProject/SimBulkNetworksFromSC/withCounts/bulkFromSC_%d_newComb3_withCounts.mat', ...
                     n)) % bulkFromSC
                         % thisNetSyms = scGeneSyms;
                         % myNet = bulkFromSC.binNet005;
        secondNet = triu(bulkFromSC.binNet005, 1);
        secondSyms = scGeneSyms;
    end
    
        % >>>>>>>> first net SC bulk Sim
    if (i > 74) && (i <= 79)
        n = i - 69 % to 10
        load(sprintf('~/resultsAndFigures/secondProject/SimBulkNetworksFromSC/bulkFromSC_%d_newComb7_noVar.mat', ...
                     n)) % bulkFromSC
                         % thisNetSyms = scGeneSyms;
                         % myNet = bulkFromSC.binNet005;
        secondNet = triu(bulkFromSC.binNet005, 1);
        secondSyms = scGeneSyms;
    end

        
        % >>>>>>>>>> TAN
        % I get the symbols for mapping, however the density comes
        % directly from the density of the mapped network, for affy it
        % is different
        if (j > 79) && (j <= 82) 
            t = j - 79  % 1 2 3 

            tissue = tissues{t};
            load( ['~/networks/tissues/' tissue '/' ...
                   'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
            thisExpGenes = logical(affyExpMat(:, t));
            secondNet = binNet(thisExpGenes, thisExpGenes);
            secondSyms = affyGeneSyms(thisExpGenes);
            clear thisExpGenes
        end

        % >>>>>>>> GTEx 
        if (j > 82) && (j <= 85) % I should fix affy density
            t = j - 82 % 1 2 3
            load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat') 
            secondNet = GTExFiveNets.nets(t).net005;
            secondSyms = ...
                GTExFiveNets.uniqueGeneSyms(GTExFiveNets.nets(t).expGenes);
        end

        % CTC
        if j == 86
            % load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
            %       'correctedBinNets_logCorrected_jusResiduals_scBasedMarkers' ...
            %       '.mat'])

         load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
              'correctedBinNets_logCorrected_jusResiduals_scBasedMarkers_redo.mat'])
            secondNet = ctc.net005;
            secondSyms = ctc.geneSyms;
        end
        
        % now I have the firstNet, firstSyms, secondNet, secondSyms
        %% The actual network comparisons ...
        [a, b] = ismember(firstSyms, secondSyms);
        subNet1 = firstNet(a, a);
        subNet2 = secondNet(b(a), b(a));
        ogc = sum(a);
        
        d1 = sum(subNet1(:)) / (ogc*(ogc-1)/2);
        d2 = sum(subNet2(:)) / (ogc*(ogc-1)/2);
        nullRatio = d1 * d2;
        
        netGeneOverlap(i, j) = ogc;
        
        netLinkCounts(i, j) = sum(subNet1(:));
        netLinkCounts(j, i) = sum(subNet2(:));
        netLinkCounts(i, i) = sum(firstNet(:));

        book = (subNet1 + subNet2) == 2;
        netLinkCountoverlap(i, j) = sum(book(:));
        d3 = sum(book(:)) / (ogc*(ogc-1)/2);
        netLikeRatiooverlap(i, j) = d3 / nullRatio;
    end
end


noResult.netLinkCounts = netLinkCounts;
noResult.netGeneOverlap = netGeneOverlap;
noResult.netLinkCountOverlap = netLinkCountoverlap;
noResult.netLikeRatioOverlap = netLikeRatiooverlap;
save('~/resultsAndFigures/secondProject/allNetworkOverlaps_noVarAdded_redo.mat', ...
     'noResult')
load('~/resultsAndFigures/secondProject/allNetworkOverlaps_noVarAdded.mat')

load('~/resultsAndFigures/secondProject/allNetworkOverlaps.mat')
h = figure
heatmap(noResult.netGeneOverlap + noResult.netGeneOverlap' + eye(86, ...
                                                  86).* 4000)
colormap(jet)

h = figure
heatmap(netLikeRatiooverlap(80:end, 80:end))
plotMat = noResult.netLikeRatioOverlap + noResult.netLikeRatioOverlap';
heatmap(log2(plotMat))

h = figure
netLCs = diag(noResult.netLinkCounts);
inNets = find(netLCs > 10000);
labels = 1:86;
plotMat = netLikeRatiooverlap(inNets, inNets) + netLikeRatiooverlap(inNets, inNets)';
heatmap(labels(inNets), labels(inNets), log2(plotMat+1));
colormap(jet)

% building labels
inNets(59:60)
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exon_V4_clusterLabels.mat'])
scLabels = clusterMeta.sortedClusterNames(1:69);
[scSortedLabels, scInds] = sort(scLabels(inNets(1:59)));

plotMat(1:59, 1:59) = plotMat(scInds, scInds);

ml = cell(76, 1);
ml(1:59) = scSortedLabels;
for i = 1:10
    ml{59+i} = sprintf('scSim_%d', i);
end

ml{70} = 'TAN_blood';
ml{71} = 'TAN_brain';
ml{72} = 'TAN_liver';
ml{73} = 'GTEx_blood';
ml{74} = 'GTEx_brain';
ml{75} = 'GTEx_liver';
ml{76} = 'GTEx_CTC';

h = figure
h = figure('units', 'centimeters', 'position', [0,0, 30, 27])
heatmap(ml, ml, log2(plotMat+1))
colormap(jet)

title(['log2 Likelihood ratio of the overlap of the network edges to the ' ...
       'null density of the network'])

figFolder = ['~/resultsAndFigures/secondProject/']
file = sprintf('%snetworkOverlaps_LRAll_boneColormap', figFolder);
set(h, 'PaperOrientation', 'landscape')
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


% 2. getting the gene pairs with disproportionally high overlap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('~/data/brainSingleCell/filDataSet_exon_V4.mat')
scGeneSyms = filDataSet.geneSyms;
clear filDataSet;
whos scGeneSyms

sumNet = zeros(length(scGeneSyms));
sumNetExc = zeros(length(scGeneSyms));
sumNetInh = zeros(length(scGeneSyms));
for n = 1:69
    n
    load(sprintf(['~/networks/allenBrainSC/rpmAllFiveNets/' ...
                  'rpm_binNets_Exon_V4_allFive_net%d.mat'], n))
    carrier = zeros(length(scGeneSyms));
    carrier(net.expGenes, net.expGenes) = net.net005;
    temp = net.clusterName{1};
    % just getting the neuronal networks
    if (length(strfind(temp, 'Inh')) || length(strfind(temp, 'Exc')))
        sumNet = sumNet + carrier; 
    end
    
    if (length(strfind(temp, 'Inh')))
        sumNetInh = sumNetInh + carrier; 
    end

    if (length(strfind(temp, 'Exc')))
        sumNetExc = sumNetExc + carrier; 
    end
end
load(['~/resultsAndFigures/secondProject/' ...
      'overlappingLinksBetweenSCNets/sumNets.mat'])

upCount = sumNet(logical(triu(ones(size(sumNet)), 1)));

book = hist(upCount, [0:69]);
h = figure
plot(log10(book+1))
hold on
plot(log10(book+1), 'o')

sib = sumNet > 15;
sibnds = sum(sib);
sum(sibnds > 0)

[a, b, c] =  find(triu(sumNet));

sa = a(c >= 15);
sb = b(c >= 15);
sc = c(c >= 15);

i = 1
scGeneSyms(sa(i))
scGeneSyms(sb(i))
i = i + 1

% how does it overlap with TSN
[a, b] = ismember(gpl570.uniqueSymbols, scGeneSyms);
subTSN = tsn(a, a);
subSib = sib(b(a), b(a));
combSyms = gpl570.uniqueSymbols(a);
combNet = ((subTSN + subSib) == 2);
combNDs = sum(combNet + combNet');
inGenes = combNDs > 0;

plotNet = combNet(inGenes, inGenes);
plotSyms = combSyms(inGenes);
h = figure
heatmap(plotSyms, plotSyms, plotNet+plotNet')

[e f g] = find(combNet);

halva = sum(sum((subTSN + subSib) == 2))
% the count of links is low (~200), but it is 18.27 times more
% likely than by chance (LR = 18.27)

% getting sorted clusters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sortedSCNames, sortedSCInds] = sort(clusterMeta.sortedClusterNames);
whos scGeneSyms
inhInds = sortedSCInds(3:27);
excInds = sortedSCInds(28:72);

inhSumNet = zeros(length(scGeneSyms));
excSumNet = zeros(length(scGeneSyms));
for i = 1:69
    n = sortedSCInds(i);
    if n < 69
        if ismember(n, inhInds)
            
            load(sprintf(['~/networks/allenBrainSC/rpmAllFiveNets/' ...
                          'rpm_binNets_Exon_V4_allFive_net%d.mat'], n))
            carrier = zeros(length(scGeneSyms));
            carrier(net.expGenes, net.expGenes) = net.net005;
            inhSumNet = inhSumNet + carrier; 
            
        end
        
        if ismember(n, excInds)
            load(sprintf(['~/networks/allenBrainSC/rpmAllFiveNets/' ...
                          'rpm_binNets_Exon_V4_allFive_net%d.mat'], n))
            carrier = zeros(length(scGeneSyms));
            carrier(net.expGenes, net.expGenes) = net.net005;
            excSumNet = excSumNet + carrier; 
        end
    end
end

sumNets.inh = sparse(inhSumNet);
sumNets.exc = sparse(excSumNet);
save(['~/resultsAndFigures/secondProject/' ...
      'overlappingLinksBetweenSCNets/sumNets.mat'], 'sumNets')

[a, b] = ismember(myGenes, scGeneSyms);
subInh = inhSumNet(b(a), b(a));
subExc = excSumNet(b(a), b(a));
result.excSum = subExc;
result.inhSum = subInh;

upCount = sumNet(logical(triu(ones(size(sumNet)), 1)));


