% In this file I do the simulations 
% 1. for a given gene, give me its average expression level in 6 cell types, with varying variance. 
% 2. code for the supplement
% 3. simulation code for the supplement
% 4. correlation of R2 and all variance stuff 

%here, I want to demonstrate that the higher the variance of the
% average expression level of a gene in a cell type, the more
% likely that the observed coexpression in the bulk tissue
% represent the correlaion of average exprssion level of the genes
% in different cell types 

% 1. for a given gene, give me its average expression level
% in 6 cell types, with varying variance. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A =[   13.0000
   10.0000
    7.0465
    7.5970
    3.7662
    3.3596
    5.6350
    4.9803
    5.5305
    5.5566]'

B =[1.0394
    5.0807
    3.8169
    2.4795
    2.6128
    6.5629
    2.3805
    1.5264
    9.0000
   11.0000]'

A = [13, 10 (randn(1,8)*2 + 4)]
B = [(randn(1,8)*2 + 4), 9, 11]
centerc = corr(A', B')
var(A)
var(B)

m = 100
err = zeros(1, 1000);
centerc = zeros(1, 1000);
sampleCorr = zeros(1, 1000);
for k = 1:1000
    k
    A = (randn(1,10)*3 + 4);
    B = (randn(1,10)*3 + 4);
    A = A - mean(A);
    B = B - mean(B);
    var(A);
    var(B);
    Gs = [A; B];

    centerc(k) = corr(A', B');

    cs = zeros(1, 100);
    for j = 1:100
        temp = rand(10,  m);
        w = temp;
        % normalizing the weight matrix
        for i = 1:m
            w(:, i) = w(:, i) ./(sum(w(:,i)));
        end

        % getting the final vectors 
        exps = Gs * w ;
        sib = corr(exps');
        cs(j) = sib(1,2);
    end

    sampleCorr(k) = sib(1, 2);

    err(k) = sqrt(sum((cs - centerc(k)).^2));
end

corr(err', abs(centerc'))

h = figure
scatter(sampleCorr, centerc)
xlim([-1, 1])
ylim([-1, 1])
figFolder = ['~/resultsAndFigures/secondProject/simCorr/']
file = sprintf('%scorrAcBc_corrEAcEBc', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 2. code for the supplement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
m = 100
CTCorr = zeros(1, 1000);
sampleCorr = zeros(1, 1000);
for k = 1:1000
    k
    A = (randn(1,10)*3 + 4);
    B = (randn(1,10)*3 + 4);

    Gs = [A; B];
    CTCorr(k) = corr(A', B');

    w = rand(10,  m);
    % normalizing the weight matrix
    for i = 1:m
        w(:, i) = w(:, i) ./(sum(w(:,i)));
    end

    % getting the final vectors 
    exps = Gs * w ;

    sib = corr(exps');
    sampleCorr(k) = sib(1, 2);
end

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

A = rand(1, 10)
B = rand(1, 10)
corr(A', B')

Ac = A - mean(A);
Bc = B - mean(B);

corr(Ac', Bc')

% getting the cosine
sum(Ac.*Bc) / (norm(Ac)*norm(Bc))

alpha = rand(1,10)
alphac = alpha - mean(alpha);

sum(alphac .*Ac)/(norm(alphac)*norm(Ac))
sum(alphac .*Bc)/(norm(alphac)*norm(Bc))

% see what is happening between cos(a, b) and cos(ac, bc)
% my guess is that, if a and b are highly correlated, cosine of the
% (a,b) and cos(ac, bc) are very close to each other, close ot 1

cs = zeros(1000, 2);
cr = zeros(1000, 1);
for i = 1:10000
    % cos(a,b)
    a = (randn(1,10)*3 + 4);
    b = (randn(1,10)*3 + 4);
    a = a + abs(min(a));
    b = b + abs(min(b));
    ac = a - mean(a);
    bc = b - mean(b);

    % cos(a, b)
    cr(i) = corr(a', b');
    cs(i, 1) = sum(a.*b) /(norm(a)*norm(b));
    cs(i, 2) = sum(ac.*bc) / (norm(ac)*norm(bc));
end

corr(cs)

% in this simulation we show that the higher the correlation of A
% and B, the smaller the angle between them, and therefore, the
% smaller will be the difference between the angle alphaA and alphaB

h = figure
scatter(cs(:,1), cs(:,2))
xlim([-1, 1])
ylim([-1, 1])
figFolder = ['~/resultsAndFigures/secondProject/simCorr/']
file = sprintf('%scosineOfABandAcBc', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% now for the variance
orVar = zeros(1,1000);
resVar = zeros(1, 1000);
norms = zeros(1, 1000);
% when n = m
v = 3
n = 10 % count of tissues
m = 100 % count of bulk samples
        %A = (sqrt(v).*randn(1,n)+4);

temp = rand(n,  m);
w = temp;
% normalizing the weight matrix
for i = 1:m
    w(:, i) = w(:, i) ./(sum(w(:,i)));
end

book = w.^2;

meanEU = zeros(1, 1000);
% a gene with the mean expression level 10
for k = 1:1000
    k
    v = rand/(rand+.1)
    A = (sqrt(v).*randn(1,n)+40);
    orVar(k) = var(A);
    
    norms(k) = norm(A);

    % temp = rand(n,  m);
    % w = temp;
    % % normalizing the weight matrix
    % for i = 1:m
    %     w(:, i) = w(:, i) ./(sum(w(:,i)));
    % end

    % EUC = zeros(m);
    % for i = 1:m-1
    %     for ii = i+1:m
    %         %            EUC(i,ii) = sqrt(sum((w(:, i) - w(:,
    %         %            ii)).^2));
    %         EUC(i, ii) = dot(w(:,i), w(:,ii))/sqrt((sum((w(:,i).^2))*sum((w(:,ii).^2))));
    %     end
    % end

    % meanEU(k) = sum(EUC(:))/(100*50);
    % meanEU(k) = sqrt(sum(var(w)'.^2));

    % getting the final vectors 
    exps = A * w;
    resVar(k) = var(exps);
end


x = [ones(1, length(orVar)) ; log2((orVar))];
[b, bint, r, rint, stats] = regress(log2(resVar'), x');

x = [ones(1, length(orVar)) ; log2(orVar)];
[b, bint, r, rint, stats] = regress(log2(resVar'), x');

[s, p] = corr(resVar', norms')

corr(log2(orVar'), log2(resVar'))
corr((orVar'), (resVar'))
h = figure
scatter(log2(orVar), log2(resVar))
figFolder = ['~/resultsAndFigures/secondProject/simCorr/']
file = sprintf('%sorVarAVarRes_varRV_m%d_n%d_sqrt', figFolder, m, n);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 3. simulation code for the supplement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
t = 10 % count of tissues
m = 100 % count of bulk samples

% normalizing the weight matrix
w = rand(t,  m); 
for i = 1:m
    w(:, i) = w(:, i) ./(sum(w(:,i)));
end

orVar = zeros(1,1000); % variance of the CT profile
obVar = zeros(1,1000); % the observed variance in bulk tissue
for k = 1:1000 % for 1000 genes
    ind = datasample(1:length(myVar), 1);
    v = myVar(ind);
    orVar(k) = v;
    A = normrnd(myMeans(ind), sqrt(v), 1, t);
    exps = A * w;
    obVar(k) = var(exps);
end
corr(obVar', orVar', 'type', 'Spearman')
h = figure
scatter(orVar, obVar, 'filled', 'MarkerFaceAlpha', 0.2)
xlim([-.25, 12])
ylim([-.01, .4])
figFolder = ['~/resultsAndFigures/secondProject/simCorr/']
file = sprintf('%sorVarAVarRes_varRV_m%d_n%d_sqrt_V01', figFolder, m, t);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% 4. correlation of R2 and all variance stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% the induced variance in the bulk tissue is correlated highly with
% the among the tissues. 

% this must show that also in my data, genes with higher variance
% should have higher R2 values 
% loading the GTEx cluster

clear
% loading data - exon
load('~/data/brainSingleCell/filDataSet_exon_V4.mat')
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exon_V4.mat'])
load(['~/data/brainSingleCell/' ...
      'dataSet_meta_filtered_exon_V4_clusterLabels.mat'])

load('~/data/GTEx/Brain_Cortex_expGenes.mat')

load(['~/resultsAndFigures/secondProject/GTExRegression/' ...
      'regressionResult_5PCA_AllenSCBasedMarker.mat'])

% for a giveen sc cluster, get the variation
scsVar = var(thisExpMat');
halva = sum(thisExpMat' > 0 );
varExist = halva > 100;
varEInds = find(varExist);
genesWithVar = filDataSet.geneSyms(varExist);
oneSCVar = zeros(sum(varExist), 1);
for i = 1:sum(varExist)
    temps = find(thisExpMat(varEInds(i), :) > 0);
    mySamples = datasample(temps, 100);
    oneSCVar(i) = var(thisExpMat(varEInds(i), mySamples));
end

[a, b] = ismember(filDataSet.geneSyms, dataSet.genes);
scsVar = scsVar(a);
myInds = (b(a));
myGenes = dataSet.genes(b(a));
[a, b] = ismember(geneList, myGenes);
markerInds = b(a);
highVar = myInds(kado);
h = figure
heatmap(dataSet.mat(highVar, :), 'GridVisible', 'off')
highVarMat = dataSet.mat(highVar, :);
sib = corr(highVarMat);
h = figure
heatmap(sib)
highVarGenes = dataSet.genes(highVar);

bulkVar = var(log2(dataSet.mat(b(a), :)+1)');
bulkMean =  mean(log2(dataSet.mat(b(a), :)+1)');

bulkVar = var((dataSet.mat(b(a), :)+1)');
bulkMean =  mean((dataSet.mat(b(a), :)+1)');

sib = var(geneListExpMat');

scVar = sib(a);
myRs = redoRes.r2(b(a));
corr(scVar', myRs)

mycorrs = zeros(1, 100);
for i = 1:100
    tempRs = result(i).r2((b(a)));
    mycorrs(i) = corr(scVar', tempRs);
end

k1 = scVar < 5;
k2 = bulkVar < .5;
k = (k1 + k2) ==2;
corr(scVar(k)', bulkVar(k)')

h = figure
scatter(log2(scVar), log2(bulkVar), 'filled', 'MarkerFaceAlpha', ...
        .05)

h = figure
scatter((bulkVar), (scVar), 'filled', 'MarkerFaceAlpha', ...
        .02)
xlim([-.25, 2])
ylim([-.25, 5])
title('coexpression of CTvar and bulkVar')
xlabel(' variance of the genes in the bulk tissue')
ylabel(' variance of the CT profiles from AllenSNC')
figFolder = ['~/resultsAndFigures/secondProject/simCorr/']
file = sprintf('%sbulkVsOneSCVar_clipped', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

[s, p] = corr(log2(scVar)',(myRs))
h = figure
scatter(log2(myRs), log2(scVar), 'filled', 'MarkerFaceAlpha', .05)
P = polyfit(myRs, log2(scVar)', 2)
yfit = P(1) * myRs.^2 + P(2)*myRs + P(3);
hold on
plot(myRs, yfit, 'r-')

[s, p] = corr(myRs, log2(bulkVar)')
h = figure
scatter((myRs), log2(bulkVar), 'filled', 'MarkerFaceAlpha', .05)

[s, p] = corr((bulkMean'), (bulkVar)')
[s, p] = corr((bulkVar'), (scVar)')
[s, p] = corr(log2(scVar'), log2(bulkVar)')
[s, p] = corr((scVar(markerInds)'), (bulkVar(markerInds))','type' ...
              ,'Spearman')
[s, p] = corr((scVar'), (bulkVar)','type' ,'Spearman')
[s, p] = corr((scsVar(a)'), (bulkVar)','type' ,'Spearman')

% Q1: Are the bulk tissue variation and the oneSCVar correlated? 
[a, b] = ismember(myGenes, genesWithVar);
[s, p] = corr(log2(bulkVar(a)'), log2(oneSCVar(b(a))))
% A: no, they are significantly negatively correlated. 
h = figure
scatter(log2(bulkVar(a)), log2(oneSCVar(b(a))), 'filled', ...
        'MarkerFaceAlpha', .05)
xlabel('log2 bulk variance')
ylabel('log2 one single cell cluster variance')
figFolder = ['~/resultsAndFigures/secondProject/simCorr/']
file = sprintf('%sbulkVsOneSCVar', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% Q2: are the bulk tissue variation and the scVar correlated?
[s, p] = corr(log2(bulkVar'+1), log2(scVar'+1))
[s, p] = corr(log2(bulkMean'+1), log2(scVar'+1))
[s, p] = corr(log2(bulkMean'), log2(bulkVar'))
[s, p] = corr(log2(bulkMean'), log2(scVar'+1))

% A: yes, they are highly correlated 
h = figure
scatter(log2(bulkVar), log2(scVar), 'filled', 'MarkerFaceAlpha', ...
        .05)
xlabel('log2 bulk variance')
ylabel('log2 variance between CTV from AllenSC')
figFolder = ['~/resultsAndFigures/secondProject/simCorr/']
file = sprintf('%sbulkVsCTVfromAllenSC', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% Q3: are bulk tissue variation and the scsVar (not filtered)
% correlated?
[a, b] = ismember(myGenes, filDataSet.geneSyms);
[s, p] = corr(log2(bulkVar'), log10(scsVar(b(a))'))
% A: yes. They are.
h = figure
scatter(log2(bulkVar), log10(scsVar(b(a))), 'filled', 'MarkerFaceAlpha', ...
        .05)
xlabel('log2 bulk variance')
ylabel('log2 variance from one SC cluster - no filter')
figFolder = ['~/resultsAndFigures/secondProject/simCorr/']
file = sprintf('%sbulkVsVarianceFromOneSC_noFilter', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

h = figure
scatter(log2(bulkVar(a)), log2(oneSCVar(b(a))), 'filled', 'MarkerFaceAlpha', .05)

h = figure
scatter(log2(bulkVar), log2(scVar), 'filled', 'MarkerFaceAlpha', .05)

h = figure
scatter(log2(myVar), log2(book))
title(sprintf('corr %.2f', s ))
xlabel('log var(genes) from SC')
ylabel('log var(genes) from GTEx')
figFolder = ['~/resultsAndFigures/secondProject/simCorr/']
file = sprintf('%slog2VarTobulk_gtex', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


h = figure
scatter(log2(book), result.r2(b(a)))
title(sprintf('corr %.2f', s ))
xlabel('log var(genes)')
ylabel('Rsqrd of the prediction')
figFolder = ['~/resultsAndFigures/secondProject/simCorr/']
file = sprintf('%slog2VarToR2_gtex', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

% so yes. genes with varying average expression among cell types,
% DO have high R2 values for the prediction model. 

% The variance of the gene between the cell type is highly
% correlated with the variance of the gene among the samples in the
% bulk tissue. Simulations show that. 

% Then we go ahead and explain the coexpression thing, I should
% model the noise thing and etc. 

% Simulation 1
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% >> fetch the group of genes with similar mean in the bulk
% tissue. Show that they have higher variance if they have higher
% CT variance

% I want the difference of the mean to be so small that we do not
% detect any correlaion between the mean and variance 

% just for figure file 21
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha1 = [13 29 73 53 17]
alpha2 = 1 - alpha1./100;

mybars = [alpha1./100; alpha2];

h = figure
bar(mybars', 'stacked')
figFolder = ['~/resultsAndFigures/secondProject/simCorr/']
file = sprintf('%salphasFigure21', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')


myCTVs = [8 10 0 3 7; 0 0 9 8 7]
book = mybars' * myCTVs
h = figure
heatmap(book')
colormap(1-bone)
figFolder = ['~/resultsAndFigures/secondProject/simCorr/']
file = sprintf('%sheatmapFigure21_03', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')

book = mybars' * myCTVs
h = figure
hold on
for i = 1:5
    plot([1:5], book(:, i))
end
xlim([0 6])
figFolder = ['~/resultsAndFigures/secondProject/simCorr/']
file = sprintf('%slinesFigure21_02', figFolder);
set(h, 'PaperOrientation', 'landscape')
%print(h, '-dpng', [file '.png'])
print(h, '-deps', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])
saveas(h, [file '.eps'], 'epsc')




