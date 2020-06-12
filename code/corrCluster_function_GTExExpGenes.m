% function [geneList2, p] = corrCluster_function(geneList, cthr,
% c clusterCutOff
% cthr corrlation cutoff

function [corrMarkers] = corrCluster_function_GTExExpGenes(geneList, ct)
    load(['~/data/GTEx/Brain_Cortex_expGenes.mat'])
    
    % sib = corr(log2(dataSet.mat + 1)');
    % upperSingle = sib(logical(triu(ones(size(sib)), 1)));
    % q = quantile(upperSingle, cthr)
    
    load('~/networks/GTEx/fiveTissues_rpmFromGeneLevel_binNets.mat')

    %    geneList = markerSet(:, 18);
    sum(geneList)
        %c = .9
    smallMat = GTExFiveNets.nets(2).net01(logical(geneList), ...
                                          logical(geneList));
    
    % throw away the ones that are not correlated at all
    sfm = smallMat + smallMat';
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
    %    dendrogram(z, size(normMat, 1));
    
    c = cluster(z, 'cutoff', ct, 'criterion', 'distance');
    
    % selecting the cluster
    m = max(c);
    objectCount = hist(c, unique(c));
    
    [a, b] = max(objectCount);
    finalList = c == b;

    geneList1 = find(geneList);
    h = figure
    heatmap(smallMat(finalList, finalList)+0);
    corrMarkers = geneList1(finalList);

    % hh = figure
    %  heatmap(book(finalList, finalList), [], [], [], 'TickAngle', 45, 'ShowAllTicks', ...
    %       true, 'Colorbar', true, 'Colormap', 'bone')
    
    % What are the odds of having this count of genes from the
    % given list connected to each other
end
