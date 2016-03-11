
%% Analysis of the different treatments
function [] = whMetaAnalysis(mainDirname,metaDataFname,timePerFrame)

% % mainDirname = '/project/cellbiology/gdanuser/collab/hall/MetaAnalysis/20140829/';
% % mainDirname = '/project/cellbiology/gdanuser/collab/hall/MetaAnalysis/20141013/';
% % metaDataFname = 'GEFsScreenMetaData20140909_all.mat';
% % metaDataFname = 'GEFsScreenMetaData20140911_all_kd50.mat';
% % metaDataFname = 'GEFsScreenMetaData20140911_all_kd40.mat';
% % metaDataFname = 'GEFsScreenMetaData20141013_kd50.mat';

% mainDirname = '/project/cellbiology/gdanuser/collab/hall/MetaAnalysis/20150504/';
% metaDataFname = 'GEFsScreenMetaData20150424_kd50.mat';
% timePerFrame = 5;

tic;

close all;

addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/utils/export_fig'));
addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/utils/notBoxPlot'));
addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/Hall'));

load([mainDirname metaDataFname]); % metaData

metaData.timePerFrame = timePerFrame;
metaData.timeToAnalyze = 200;%450;%200; % minutes
metaData.timePartition = 4;
metaData.spatialPartition = 3;
metaData.initialTimePartition = 0; %2

kymographDir = [mainDirname 'kymographs/'];
pcaDir = [mainDirname 'pca/'];
pcaDayDir = [pcaDir 'day/'];
pcaGeneDir = [pcaDir 'gene/'];
pcaDayGeneDir = [pcaDir 'dayGene/'];
associationDir = [mainDirname 'association/'];
healingRateDir = [mainDirname 'healingRate/analysis/'];
varianceDir = [mainDirname 'variance/'];
dayDir = [mainDirname 'day/'];
daySVMDir = [mainDirname 'svm/'];
dayClusterDir = [mainDirname 'clusters/'];
dayGeneDir = [mainDirname 'day/gene/'];
shStatsDir = [mainDirname 'shSeqEffStats/'];
micrscopeRepeatErrorStatsDir = [mainDirname 'micrscopeError/'];
targetsDir = [mainDirname 'targets/'];
dayGeneControlKymographDir = [mainDirname 'dayGeneControlKymograph/'];
dayGenePlotsDir = [mainDirname 'dayGenePlots/'];
dayWellReplicatesKymographDir = [mainDirname 'dayWellReplicatesKymographs/'];
combinedPropertiesDir = [mainDirname 'combinedProperties/'];
combinedPropertiesAnglesPcsDir = [mainDirname 'combinedProperties/anglesPC/'];
controlInterdayAssessmentDir = [mainDirname 'controlInterdayAssessment/'];

plithotaxisDir = [mainDirname 'plithotaxis/'];
plithotaxisOutDir = [mainDirname 'plithotaxisOut/'];

correctMotionDir = [mainDirname 'correctMotion/'];

RhoGTPasesDir = [mainDirname 'RhoGTPases/'];

if ~exist(kymographDir,'dir')
    mkdir(kymographDir);
end

if ~exist(pcaDir,'dir')
    mkdir(pcaDir);
end

if ~exist(pcaDayDir,'dir')
    mkdir(pcaDayDir);
end

if ~exist(pcaDayGeneDir,'dir')
    mkdir(pcaDayGeneDir);
end

if ~exist(pcaGeneDir,'dir')
    mkdir(pcaGeneDir);
end

if ~exist(associationDir,'dir')
    mkdir(associationDir);
end

if ~exist(healingRateDir,'dir')
    mkdir(healingRateDir);
end

if ~exist(varianceDir,'dir')
    mkdir(varianceDir);
end

if ~exist(dayDir,'dir')
    mkdir(dayDir);
end

if ~exist(daySVMDir,'dir')
    mkdir(daySVMDir);
end

if ~exist(dayClusterDir,'dir')
    mkdir(dayClusterDir);
end

if ~exist(dayGeneDir,'dir')
    mkdir(dayGeneDir);
end

if ~exist(shStatsDir,'dir')
    mkdir(shStatsDir);
end

if ~exist(micrscopeRepeatErrorStatsDir,'dir')
    mkdir(micrscopeRepeatErrorStatsDir);
end

if ~exist(targetsDir,'dir')
    mkdir(targetsDir);
end

if ~exist(dayGeneControlKymographDir,'dir')
    mkdir(dayGeneControlKymographDir);
end

if ~exist(dayGenePlotsDir,'dir')
    mkdir(dayGenePlotsDir);
end

if ~exist(dayWellReplicatesKymographDir,'dir')
    mkdir(dayWellReplicatesKymographDir);
end

if ~exist(combinedPropertiesDir,'dir')
    mkdir(combinedPropertiesDir);
end

if ~exist(combinedPropertiesAnglesPcsDir,'dir')
    mkdir(combinedPropertiesAnglesPcsDir);
end

if ~exist(controlInterdayAssessmentDir,'dir')
    mkdir(controlInterdayAssessmentDir);
end

if ~exist(plithotaxisOutDir,'dir')
    mkdir(plithotaxisOutDir);
end

if ~exist(RhoGTPasesDir,'dir')
    mkdir(RhoGTPasesDir);
end

[labels, strLabels] = getLabels(metaData); % labels are combination of gene & shSeq

%% META ANALYSIS
%   1. KD efficiency
%   2. Micrscopy error correction
%   3. Extract high-dimensional features from the data
%       a. Spatiotemporal associations (time delay) between measures (commented out)
%       b. PCA
%   4. Visualize PCs
%   5. Gene abalysis (not day as atomic unit)
%   6. Correlate PCs of each measure to wound healing rate (healingRate directory)
%   7. MetaDayAnalysis
%   8. Shear strain analysis: correlation shear strain events and motion
%       a. CDC42 vs. Control
%       b. All experiments
%   9. RhoGTPAses analysis

flags.kdEfficiency = 0;
flags.correctMicroscopeError = 0;
flags.visualizePCs = 0;
flags.geneAnalysis = 0;
flags.corrPCsHealingRate = 0;
flags.metaDayAnalysis = 1;
flags.shearStrainAnalysis = 0;
flags.RhoGTPAsesAnalysis = 0;



%% 1. KD efficiency
% Statistics of KD efficiency of all sh-sequences
if flags.kdEfficiency
    whMetaKnockdownEfficiency(metaData,strLabels,shStatsDir);
end

%% 2. Statistics (printouts) of correction microscope repeat errors
if flags.correctMicroscopeError
    whMetaMicrscopeRepeatCorrectionStats(metaData,correctMotionDir);
end

%% 3. Extract high-dimensional features from the data
%       a. Spatiotemporal associations (time delay) between measures (commented out)
%       b. PCA
allFeatsFname = [mainDirname '/allFeatures.mat'];
reuseFeats = true; 

if exist(allFeatsFname,'file') && reuseFeats
    load(allFeatsFname);
else
    [allFeatures] = getAllFeatureVectors(mainDirname,metaData);
    save(allFeatsFname,'allFeatures');
end

% patch for fix directionality > 8
allFeatures.directionalityFeats.features(allFeatures.directionalityFeats.features > 10) = 10;
save(allFeatsFname,'allFeatures');

% associations = spatioTemporalAssociations(allFeatures,strLabels,mainDirname,timePerFrame);

% TODO: fill empty features? Use machine learning to predict empty spots
% one-by-one based on previous and future spatiotemporal data.

% Also plots gene (no sh-seq data);
plotBasicPCAFname = [mainDirname '/pcaFeatures.mat'];
reuseFeats = true; 

if exist(plotBasicPCAFname,'file') && reuseFeats
    load(plotBasicPCAFname);
else
    out = pcaFeatures(allFeatures,metaData,strLabels,mainDirname); % includes subjective clustering
    save(plotBasicPCAFname,'out');
end

%%   4. Visualize PCs
if flags.visualizePCs
    close all;
    visualizePCs(out,metaData,mainDirname);
end



%%  5. Gene abalysis (not day as atomic unit)
if flags.geneAnalysis
    close all;
    % % DOES NOT USE DAILY CONTROLS!
    % geneVsControl(out,strLabels,metaData,mainDirname);
    
    % For each gene & day plot all sh-seqences + control
    geneAnalysis(out,strLabels,metaData,mainDirname);% gene x 3 shRNAs, Control, pSup of the same gene (defined by day), BetaPIX & ARHGEF12 as reference points
    
    
    % similar to geneAnalysis just not using the sh-sequences
    dayAnalysisControl(out,strLabels,metaData,mainDirname);
    
    % % same as dayAnalysisControl? Why is it commented out?
    % dayAnalysisGene(out,strLabels,metaData,mainDirname);
end

%%  6. Correlate PCs of each measure to wound healing rate
if flags.corrPCsHealingRate
    % Correlates PC1 to monolayer speed (same for PC2)
    close all;
    healingRateAnalysis(out,metaData,mainDirname);
end

%%  7. Meta day analysis
if flags.metaDayAnalysis
    whMetaDayAnalysis(allFeatures,out.healingRate,strLabels,metaData,mainDirname);
end

%%  8. Shear strain analysis: correlation shear strain events and motion
%       a. CDC42 vs. Control
%       b. All experiments
if flags.shearStrainAnalysis
    [p,x0,x1] = whShearStrainCntrlCdc42(metaData,mainDirname); % motion-shear strain for all controls + CDC42
    whMotionShearStrainGeneDaySeq(metaData,mainDirname,p,x0,x1);
end

%%  9. RhoGTPAses analysis
if flags.RhoGTPAsesAnalysis
    close all;
    % TODO: implement this
    RhoGTPasesAnalysis(out,metaData,RhoGTPasesDir);
end

% TODO? take pooled correlations (for each well) between the control and
% the KD phenotype
end






%%

%%
function [labels, strLabels] = getLabels(metaData)
nLabels = length(metaData.groupsByTreatments);
nExps = metaData.N;

labels = zeros(1,nExps);
strLabels = cell(1,nExps);

for curLabel = 1 : nLabels
    labels(metaData.groupsByTreatments{curLabel}.inds) = curLabel;
    strLabels(metaData.groupsByTreatments{curLabel}.inds) = metaData.groupsByTreatments{curLabel}.treatment;
end

end

%% 
function [allFeatures] = getAllFeatureVectors(mainDirname,metaData)
kymographDir = [mainDirname 'kymographs/'];
healingRateDir = [mainDirname 'healingRate/'];
measures = {'speed','directionality','strainRate','coordination'};

allFeatures.speedFeats = extractFeatures(kymographDir,metaData,measures{1});
allFeatures.directionalityFeats = extractFeatures(kymographDir,metaData,measures{2});
% allFeatures.strainRateFeats = extractFeatures(kymographDir,metaData,measures{3});
allFeatures.coordinationFeats = extractFeatures(kymographDir,metaData,measures{4});

allFeatures.healingRate = getHealingRate(healingRateDir,metaData);
% TODO: add dfeatures
end



%%
function [feats] = extractFeatures(kymographDir,metaData,measureStr)
feats.kymographs = cell(1,metaData.N);
feats.features = cell(1,metaData.N);
for i = 1 : metaData.N
    load([kymographDir measureStr '/' metaData.fnames{i} '_' measureStr 'Kymograph.mat']);
    if exist('speedKymograph','var')
        kymograph = speedKymograph;
    end
    if exist('directionalityKymograph','var')
        kymograph = directionalityKymograph;
    end
    if exist('strainRateKymograph','var')
        kymograph = strainRateKymograph;
    end
    if exist('coordinationKymograph','var')
        kymograph = coordinationKymograph;
    end
    
    assert(exist('kymograph','var') > 0);
    
    % TODO: deal with different kymographs
    feats.kymographs{i} = kymograph;    
end

[minDistFromWound,allDistancesFromWound] = getDistanceFromWound(feats.kymographs);

[feats.features, feats.kymographs]= getFeatures(feats.kymographs,metaData,minDistFromWound);
end

%% Find maximal distance for which all values are calculated
function [minDistFromWound,allDistancesFromWound] = getDistanceFromWound(kymographs)
minDistFromWound = inf;
n = length(kymographs);
allDistancesFromWound = zeros(1,n);
for i = 1 : n
    curKymograph = kymographs{i};
    ind = find(isnan(curKymograph(:,1)),1,'first')-1;
    if isempty(ind)
        ind = size(curKymograph,1);
    end
    allDistancesFromWound(i) = ind;
    if ~isempty(ind)
        minDistFromWound = min(minDistFromWound,ind);
    end
end
end

%%
function [features kymographs] = getFeatures(kymographs,metaData,nDist)
n = length(kymographs);
nTime = floor(metaData.timeToAnalyze/metaData.timePerFrame);
nFeats = (metaData.timePartition-metaData.initialTimePartition) * metaData.spatialPartition; %metaData.timePartition * metaData.spatialPartition;
features = zeros(nFeats,n);

for i = 1 : n
    curKymograph = kymographs{i};    
    ys = 1 : floor(nDist/(metaData.spatialPartition)) : (nDist+1);
    xs = 1 : floor(nTime/(metaData.timePartition)) : (nTime+1);
    curFeatI = 0;
    for y = 1 : metaData.spatialPartition
        for x = (1+metaData.initialTimePartition) : metaData.timePartition%x = 1 : metaData.timePartition
            curFeatI = curFeatI + 1;
            values = curKymograph(ys(y):(ys(y+1)-1),xs(x):(xs(x+1)-1));
            values = values(~isinf(values));
            values = values(~isnan(values));
            features(curFeatI,i) = mean(values(:));
        end
    end
    kymographs{i} = curKymograph(1:nDist,1:nTime);
end

end

%% 
function [allHealingRates] = getHealingRate(healingRateDir,metaData)
allHealingRates = zeros(1,metaData.N);
nTime = floor(metaData.timeToAnalyze/metaData.timePerFrame);

for i = 1 : metaData.N
    clear averageHealingRate;
    load([healingRateDir '/' metaData.fnames{i} '_healingRate.mat'])
    allHealingRates(i) = averageHealingRate(nTime);    
end
end

%%
function [out] = pcaFeatures(allFeatures,metaData,strLabels,mainDirname)
out.healingRate = allFeatures.healingRate;

[out.speed.coeff,out.speed.score,out.speed.latent] = pca(allFeatures.speedFeats.features');
out.speed.accVariance = cumsum(out.speed.latent)./sum(out.speed.latent);
[out.directional.coeff,out.directional.score,out.directional.latent] = pca(allFeatures.directionalityFeats.features');
out.directional.accVariance = cumsum(out.directional.latent)./sum(out.directional.latent);
% [out.strainRate.coeff,out.strainRate.score,out.strainRate.latent] = pca(allFeatures.strainRateFeats.features');
% out.strainRate.accVariance = cumsum(out.strainRate.latent)./sum(out.strainRate.latent);
[out.coordination.coeff,out.coordination.score,out.coordination.latent] = pca(allFeatures.coordinationFeats.features');
out.coordination.accVariance = cumsum(out.coordination.latent)./sum(out.coordination.latent);
% TODO: add more measures

[out.axesSpeed, out.dataRepSpeed] = pcaPlot(out.speed,metaData,strLabels,'Speed',mainDirname);
[out.axesDirectionality, out.dataRepDirectionality] = pcaPlot(out.directional,metaData,strLabels,'Directionality',mainDirname);
% [out.axesStrainRate, out.dataRepStrainRate] = pcaPlot(out.strainRate,metaData,strLabels,'StrainRate',mainDirname);
[out.axesCoordination, out.dataRepCoordination] = pcaPlot(out.coordination,metaData,strLabels,'Coordination',mainDirname);
close all;
% TODO: add more measures
end

%% 
function [newAxes, extremeData] = pcaPlot(data,metaData,strLabels,measureStr,mainDirname)
close all;
[pc1MaxVal,pc1MaxInd] = max(data.score(:,1));
[pc2MaxVal,pc2MaxInd] = max(data.score(:,2));
[pc1MinVal,pc1MinInd] = min(data.score(:,1));
[pc2MinVal,pc2MinInd] = min(data.score(:,2));
[pc1MedVal, pc1MedInd] = min(abs(median(data.score(:,1)-data.score(:,1))));
[pc2MedVal, pc2MedInd] = min(abs(median(data.score(:,2)) - data.score(:,2)));
extremeData.pc1Max.val = [pc1MaxVal, data.score(pc1MaxInd,2)];
extremeData.pc2Max.val = [data.score(pc2MaxInd,1), pc2MaxVal];
extremeData.pc1Min.val = [pc1MinVal, data.score(pc1MinInd,2)];
extremeData.pc2Min.val = [data.score(pc2MinInd,1),pc2MinVal];
extremeData.pc1Med.val = [pc1MedVal, data.score(pc1MedInd,2)];
extremeData.pc2Med.val = [data.score(pc2MedInd,1),pc2MedVal];

extremeData.pc1Max.ind = pc1MaxInd;
extremeData.pc2Max.ind = pc2MaxInd;
extremeData.pc1Min.ind = pc1MinInd;
extremeData.pc2Min.ind = pc2MinInd;
extremeData.pc1Med.ind = pc1MedInd;
extremeData.pc2Med.ind = pc2MedInd;

extremeData.pc1Max.label = strLabels{pc1MaxInd};
extremeData.pc2Max.label = strLabels{pc2MaxInd};
extremeData.pc1Min.label = strLabels{pc1MinInd};
extremeData.pc2Min.label = strLabels{pc2MinInd};
extremeData.pc1Med.label = strLabels{pc1MedInd};
extremeData.pc2Med.label = strLabels{pc2MedInd};

colorsPerm = 'ymcrgbk'; nColors = length(colorsPerm);
markersPerm ='os^Vph><+*X';
nTreats = length(metaData.groupsByTreatments);

h = figure;
xlabel('1st PC','FontSize',15);
ylabel('2nd PC','FontSize',15);
hold all;
for i = 1 : nTreats
    inds = metaData.groupsByTreatments{i}.inds;
    plot(data.score(inds,1)',data.score(inds,2)',sprintf('%s',markersPerm(mod(i,nColors) + 1)),'MarkerFaceColor',sprintf('%s',colorsPerm(mod(i,nColors)+1)),'MarkerSize',10,...
        'DisplayName',metaData.groupsByTreatments{i}.treatment);    
end
legend show;
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',32);

xlim = get(haxes,'XLim');
ylim = get(haxes,'YLim');
newAxes = [xlim, ylim];

pcaFname = [mainDirname 'pca/' measureStr 'PCA_legend.jpg'];
eval(sprintf('print -djpeg %s', pcaFname));

legend off;

hold off;

pcaFname = [mainDirname 'pca/' measureStr 'PCA.jpg'];
eval(sprintf('print -djpeg %s', pcaFname));

%% Ignore different sh sequences, just label per gene
colorsPerm = 'ymcrgbk'; nColors = length(colorsPerm);
markersPerm ='os^Vph><+*X';

uniqueGenes = whGetUniqueGeneLabels(strLabels);%new
nTreats = length(uniqueGenes);%new

h = figure;
xlabel('1st PC','FontSize',32);
ylabel('2nd PC','FontSize',32);
hold all;
for i = 1 : nTreats
    inds = strcmp(whToPrefix(strLabels),uniqueGenes{i});    
    plot(data.score(inds,1)',data.score(inds,2)',sprintf('%s',markersPerm(floor(i/nColors) + 1)),'MarkerFaceColor',sprintf('%s',colorsPerm(mod(i,nColors)+1)),'MarkerSize',10,...
        'DisplayName',uniqueGenes{i});    
end

axis(newAxes);
legend('Location','EastOutside');

pcaFname = [mainDirname 'pca/' measureStr 'PCA_Gene_legend.jpg'];
eval(sprintf('print -djpeg %s', pcaFname));

legend off;

haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',32);

hold off;

pcaFname = [mainDirname 'pca/' measureStr 'PCA_Gene_.jpg'];
eval(sprintf('print -djpeg %s', pcaFname));

close all;

%% Each gene versus the rest of the data
doGeneVsRest = false;
if doGeneVsRest
    for i = 1 : nTreats
        h = figure;
        xlabel('1st PC','FontSize',32);
        ylabel('2nd PC','FontSize',32);
        hold all;
        
        plot(data.score(:,1)',data.score(:,2)','ok','MarkerFaceColor','g','MarkerSize',8,'DisplayName','Other GEFs');
        indsControl = strcmp(whToPrefix(strLabels),'Control');
        plot(data.score(indsControl,1)',data.score(indsControl,2)','ok','MarkerFaceColor','b','MarkerSize',10,'DisplayName','Control');
        indNT = strcmp(whToPrefix(strLabels),'NT');
        plot(data.score(indNT,1)',data.score(indNT,2)','ok','MarkerFaceColor','k','MarkerSize',10,'DisplayName','NT');
        indsGene = strcmp(whToPrefix(strLabels),uniqueGenes{i});
        plot(data.score(indsGene,1)',data.score(indsGene,2)','ok','MarkerFaceColor','r','MarkerSize',12,'DisplayName',uniqueGenes{i});
        
        
        axis(newAxes);
        
        %     pcaFname = [mainDirname 'pca/' measureStr 'PCA_Gene_' uniqueGenes{i} '_legend.jpg'];
        %     eval(sprintf('print -djpeg %s', pcaFname));
        %
        %     legend off;
        
        haxes = get(h,'CurrentAxes');
        set(haxes,'FontSize',32);
        
        legend('show');
        
        hold off;
        
        pcaFname = [mainDirname 'pca/' measureStr 'PCA_Gene_' uniqueGenes{i} '.jpg'];
        eval(sprintf('print -djpeg %s', pcaFname));
        
        close all;
    end
end




%% by the subjective clustering:
doSubjectiveClustering = false;
if doSubjectiveClustering
    %     cluster1 = {'ARHGEF18', 'Fgd6','Tuba'};
    %     cluster2 = {'ARHGEF3', 'ARHGEF4', 'ARHGEF12',...
    %         'ARHGEF15', 'Fgd3', 'MCF2L2', 'SGEF', 'TIAM1'};
    %     cluster3 = {'beta-PIX', 'Cdc42', 'Rap1', 'Sos1'};
    
    
    cluster1 = {'ARHGEF18', 'Fgd6','Tuba','Vav2','Sos2','BCR','ABR'};
    cluster2 = {'ARHGEF12','MCF2L2', 'SGEF', 'TIAM1','alpha-PIX ','FGD2','FARP2','NGEF','ARHGEF3','ITSN1','FGD1'};
    cluster3 = {'Cdc42','beta-PIX', 'Sos1','Trio','Rac1'}; 
    cluster4 = {'FARBIN'};
    cluster5 = {'NT'};
    
    h = figure;
    xlabel('1st PC','FontSize',32);
    ylabel('2nd PC','FontSize',32);
    hold all;
    
    indsControl = strcmp(whToPrefix(strLabels),'Control');
    plot(data.score(indsControl,1)',data.score(indsControl,2)','ob','MarkerFaceColor','b','MarkerSize',7,'DisplayName','Control');
    
    inds1 = false(1,length(strLabels));
    for i1 = 1 : length(cluster1)
        inds1 = (inds1 | strcmp(whToPrefix(strLabels),cluster1{i1}));
    end
    plot(data.score(inds1,1)',data.score(inds1,2)','ok','MarkerFaceColor','g','MarkerSize',7,'DisplayName',clusterToString(cluster1));
    
    inds2 = false(1,length(strLabels));
    for i2 = 1 : length(cluster2)
        inds2 = (inds2 | strcmp(whToPrefix(strLabels),cluster2{i2}));
    end
    plot(data.score(inds2,1)',data.score(inds2,2)','ok','MarkerFaceColor','y','MarkerSize',7,'DisplayName',clusterToString(cluster2));
    
    inds3 = false(1,length(strLabels));
    for i3 = 1 : length(cluster3)
        inds3 = (inds3 | strcmp(whToPrefix(strLabels),cluster3{i3}));
    end
    plot(data.score(inds3,1)',data.score(inds3,2)','ok','MarkerFaceColor','r','MarkerSize',7,'DisplayName',clusterToString(cluster3));
    
    %     inds4 = false(1,length(strLabels));
    %     for i4 = 1 : length(cluster4)
    %         inds4 = (inds4 | strcmp(whToPrefix(strLabels),cluster4{i4}));
    %     end
    %     plot(data.score(inds4,1)',data.score(inds4,2)','ok','MarkerFaceColor','m','MarkerSize',7,'DisplayName',clusterToString(cluster4));
    %
    %     inds5 = false(1,length(strLabels));
    %     for i5 = 1 : length(cluster5)
    %         inds5 = (inds5 | strcmp(whToPrefix(strLabels),cluster5{i5}));
    %     end
    %     plot(data.score(inds5,1)',data.score(inds5,2)','ok','MarkerFaceColor','c','MarkerSize',5,'DisplayName',clusterToString(cluster5));
    
    
    axis(newAxes);
    legend('Location','SouthOutside');
    
    pcaFname = [mainDirname 'pca/' measureStr 'PCA_Cluster_legend.jpg'];
    eval(sprintf('print -djpeg %s', pcaFname));
    
    legend off;
    
    haxes = get(h,'CurrentAxes');
    set(haxes,'FontSize',32);
    
    hold off;
    
    pcaFname = [mainDirname 'pca/' measureStr 'PCA_Cluster_.jpg'];
    eval(sprintf('print -djpeg %s', pcaFname));
end

end

%%
function visualizePCs(out,metaData,mainDirname)
pcPlot(out.speed,'Speed',metaData,mainDirname); close all;
pcPlot(out.directional,'Directionality',metaData,mainDirname); close all;
% pcPlot(out.strainRate,'StrainRate',metaData,mainDirname); close all;
pcPlot(out.coordination,'Coordination',metaData,mainDirname); close all;
end

function pcPlot(data,measureStr,metaData,mainDirname)
nfeats = metaData.timePartition * metaData.spatialPartition;
iPc = 1;
fprintf(sprintf('%s PC2 = %.2f %.2f %.2f %.2f\n',measureStr,data.coeff(1,2),data.coeff(2,2),data.coeff(3,2),data.coeff(4,2)))
while(data.accVariance(iPc) < 0.96)
    h = figure;
    xlabel('Coefficients','FontSize',32);
    ylabel('Weight','FontSize',32);
    title(sprintf('PC #%d: %.2f variance',iPc,data.accVariance(iPc)),'FontSize',32)
    hold all;
    plot(1:nfeats,data.coeff(:,iPc),'sr','MarkerFaceColor','r','MarkerSize',10);
    
    haxes = get(h,'CurrentAxes');
    set(haxes,'XLim',[0,nfeats+1]);
    set(haxes,'XTick',1:3:nfeats);
    set(haxes,'XTickLabel',1:3:nfeats);
    set(haxes,'YLim',[-0.5,0.6]);
    set(haxes,'YTick',-0.5:0.5:0.5);
    set(haxes,'YTickLabel',-0.5:0.5:0.5);        
    set(haxes,'FontSize',32);
    
    hold off;
    
    pcaFname = [mainDirname 'pcs/' measureStr '_PC_' num2str(iPc) '.bmp'];
    
    if ~exist([mainDirname 'pcs/'],'dir')
        unix(sprintf('mkdir %s',[mainDirname 'pcs/']));
    end
    
    eval(sprintf('print -dbmp16m  %s', pcaFname));            
    
    %     % ----------------- TIME -----------------
    %     h = figure;
    %     xlabel('Space','FontSize',32);
    %     ylabel('Weight','FontSize',32);
    %     title(sprintf('PC #%d: %.2f variance',iPc,data.accVariance(iPc)),'FontSize',32)
    %     hold all;
    %     plot(1:4,data.coeff(1:4,iPc),'s','Color',[0,0.4,0],'MarkerFaceColor',[0,0.4,0],'MarkerSize',10);
    %     plot(1:4,data.coeff(5:8,iPc),'s','Color',[0,0.7,0],'MarkerFaceColor',[0,0.7,0],'MarkerSize',10);
    %     plot(1:4,data.coeff(9:12,iPc),'s','Color',[0,1,0],'MarkerFaceColor',[0,1,0],'MarkerSize',10);
    %     legend('time1','time2','time3','FontSize',32,'Location','NorthEastOutside');
    %
    %
    %     haxes = get(h,'CurrentAxes');
    %     set(haxes,'XLim',[0,5]);
    %     set(haxes,'XTick',1:4);
    %     set(haxes,'XTickLabel',1:4);
    %     set(haxes,'YLim',[-0.5,0.6]);
    %     set(haxes,'YTick',-0.5:0.5:0.5);
    %     set(haxes,'YTickLabel',-0.5:0.5:0.5);
    %     set(haxes,'FontSize',32);
    %
    %     hold off;
    %
    %     pcaFname = [mainDirname 'pcs/' measureStr '_PC_' num2str(iPc) '_Time.bmp'];
    %
    %     eval(sprintf('print -dbmp16m  %s', pcaFname));
    %
    %     % ----------------------------------
    %     h = figure;
    %     xlabel('Time','FontSize',32);
    %     ylabel('Weight','FontSize',32);
    %     title(sprintf('PC #%d: %.2f variance',iPc,data.accVariance(iPc)),'FontSize',32)
    %     hold all;
    %     plot(1:3,data.coeff(1:4:12,iPc),'s','Color',[0,0.4,0],'MarkerFaceColor',[0,0.4,0],'MarkerSize',10);
    %     plot(1:3,data.coeff(2:4:12,iPc),'s','Color',[0,0.6,0],'MarkerFaceColor',[0,0.6,0],'MarkerSize',10);
    %     plot(1:3,data.coeff(3:4:12,iPc),'s','Color',[0,0.8,0],'MarkerFaceColor',[0,0.8,0],'MarkerSize',10);
    %     plot(1:3,data.coeff(4:4:12,iPc),'s','Color',[0,1,0],'MarkerFaceColor',[0,1,0],'MarkerSize',10);
    %     legend('space1','space2','space3','space4','FontSize',32,'Location','NorthEastOutside');
    %
    %     haxes = get(h,'CurrentAxes');
    %     set(haxes,'XLim',[0,4]);
    %     set(haxes,'XTick',1:3);
    %     set(haxes,'XTickLabel',1:3);
    %     set(haxes,'YLim',[-0.5,0.6]);
    %     set(haxes,'YTick',-0.5:0.5:0.5);
    %     set(haxes,'YTickLabel',-0.5:0.5:0.5);
    %     set(haxes,'FontSize',32);
    %
    %     hold off;
    %
    %     pcaFname = [mainDirname 'pcs/' measureStr '_PC_' num2str(iPc) '_Space.bmp'];
    %
    %     eval(sprintf('print -dbmp16m  %s', pcaFname));
    %
    
    iPc = iPc + 1;
end
end


%% Compare each gene seperately to the Control
% function [] = geneVsControl(out,strLabels,mainDirname)
% geneVsControlProperty(out.speed,strLabels,mainDirname,'Speed',out.axesSpeed);
% geneVsControlProperty(out.directional,strLabels,mainDirname,'Directionality',out.axesDirectionality);
% geneVsControlProperty(out.strainRate,strLabels,mainDirname,'StrainRate',out.axesStrainRate);
% geneVsControlProperty(out.coordination,strLabels,mainDirname,'Coordination',out.axesCoordination);
% close all;
% end

% function [] = geneVsControlProperty(data,strLabels,mainDirname,propertyStr,newAxes)
% close all;
% uniqueGenes = whGetUniqueGeneLabels(strLabels);
% N = length(uniqueGenes);
% 
% strControl = 'Control';
% indsControl = strcmp(whToPrefix(strLabels),strControl);
% 
% if isempty(find(indsControl))
%     strControl = 'DMSO';
%     indsControl = strcmp(whToPrefix(strLabels),strControl);
% end
% 
% markersPerm ='s^h*X+';
% for gene = 1 : N
%     [seqStr, seqInds] = whGetSequencesStr(strLabels,uniqueGenes{gene});
%     
%     h = figure;
%     xlabel('1st PC','FontSize',32);
%     ylabel('2nd PC','FontSize',32);
%     hold all;
%     
%     plot(data.score(indsControl,1)',data.score(indsControl,2)','ob','MarkerFaceColor','b','MarkerSize',10,'DisplayName','Control');
%     
%     for seq = 1 : length(seqStr)
%         legendStr = [uniqueGenes{gene} '_{' seqStr{seq} '}'];
%         plot(data.score(seqInds{seq},1)',data.score(seqInds{seq},2)',sprintf('%sr',markersPerm(seq)),'MarkerFaceColor','r','MarkerSize',10,'DisplayName',legendStr);
%     end
%     
%     legend show;
%     haxes = get(h,'CurrentAxes');
%     axis(newAxes);
%     set(haxes,'FontSize',32);
%     
%     pcaFname = [mainDirname 'pca/' propertyStr 'PCA_Control_' uniqueGenes{gene} '.jpg'];
%     eval(sprintf('print -djpeg %s', pcaFname));
% end
% end



%% Spatiotemporal associations between different properties (time delay)
function [associations] = spatioTemporalAssociations(allFeatures,strLabels,mainDirname,timePerFrame)
uniqueGenes = whGetUniqueGeneLabels(strLabels);
nGenes = length(uniqueGenes);

timeShifts = -15:15;
nFrames = 200/5;

associations.gene = cell(1,nGenes);
associations.strainRateDirectionality = cell(1,nGenes);
associations.directionalityCoordination = cell(1,nGenes);
associations.strainRateCoordination = cell(1,nGenes);
for gene = 1 : nGenes
    indsGene = strcmp(whToPrefix(strLabels),uniqueGenes{gene});
    associations.gene = uniqueGenes{gene};
    associations.strainRateDirectionality{gene} = kymographAssociation(allFeatures.strainRateFeats.kymographs(indsGene),allFeatures.directionalityFeats.kymographs(indsGene),timeShifts,nFrames);
    associations.directionalityCoordination{gene} = kymographAssociation(allFeatures.directionalityFeats.kymographs(indsGene),allFeatures.coordinationFeats.kymographs(indsGene),timeShifts,nFrames);
    associations.strainRateCoordination{gene} = kymographAssociation(allFeatures.strainRateFeats.kymographs(indsGene),allFeatures.coordinationFeats.kymographs(indsGene),timeShifts,nFrames);        
    
    srdir = associations.strainRateDirectionality{gene};
    [mSrDir,iSrDir] = max(abs(srdir));
    SrDirCorr = max(srdir(iSrDir),0);
    
    dircoord = associations.directionalityCoordination{gene};
    [mDirCoord,iDirCoord] = max(abs(dircoord));
    DirCoordCorr = max(dircoord(iDirCoord),0);
    
    srcoord = associations.strainRateCoordination{gene};
    [mSrCoord,iSrCoord] = max(abs(srcoord));
    SrCoordCorr = max(srdir(iSrCoord),0);
    
    fprintf('***\n\nGene %s:\n strain rate - directionality %.3f (%d)\n directionality - coordination %.3f (%d)\n strain rate - coordination %.3f (%d)\n',...
        uniqueGenes{gene},SrDirCorr,iSrDir-16,DirCoordCorr,iDirCoord-16,SrCoordCorr,iSrCoord-16);
    
    h = figure;
    xlabel('Time lag (minutes)','FontSize',32);
    ylabel('Cross correlation','FontSize',32);
    hold all;
    plot(timeShifts.*timePerFrame,srdir,'ob','MarkerFaceColor','b','MarkerSize',10,'DisplayName','Control');            
    haxes = get(h,'CurrentAxes');
    set(haxes,'XLim',[timeShifts(1),timeShifts(end)].*timePerFrame);
    set(haxes,'XTick',(timeShifts(1):5:timeShifts(end))*timePerFrame);
    set(haxes,'XTickLabel',(timeShifts(1):5:timeShifts(end))*timePerFrame);
    set(haxes,'YLim',[-1,1]);
    set(haxes,'YTick',-1:0.5:1);
    set(haxes,'YTickLabel',-1:0.5:1);        
    set(haxes,'FontSize',32);                
    associationFname = [mainDirname 'association/StrainRateDirectionality_' uniqueGenes{gene} '.bmp'];
    eval(sprintf('print -dbmp16m %s', associationFname));    
    hold off;
    
    h = figure;
    xlabel('Time lag (minutes)','FontSize',32);
    ylabel('Cross correlation','FontSize',32);
    hold all;
    plot(timeShifts.*timePerFrame,dircoord,'ob','MarkerFaceColor','b','MarkerSize',10,'DisplayName','Control');            
    haxes = get(h,'CurrentAxes');
    set(haxes,'XLim',[timeShifts(1),timeShifts(end)].*timePerFrame);
    set(haxes,'XTick',(timeShifts(1):5:timeShifts(end))*timePerFrame);
    set(haxes,'XTickLabel',(timeShifts(1):5:timeShifts(end))*timePerFrame);
    set(haxes,'YLim',[-1,1]);
    set(haxes,'YTick',-1:0.5:1);
    set(haxes,'YTickLabel',-1:0.5:1);        
    set(haxes,'FontSize',32);                
    associationFname = [mainDirname 'association/DirectionalityCoordination_' uniqueGenes{gene} '.bmp'];
    eval(sprintf('print -dbmp16m %s', associationFname));    
    hold off;
    
    h = figure;
    xlabel('Time lag (minutes)','FontSize',32);
    ylabel('Cross correlation','FontSize',32);
    hold all;
    plot(timeShifts.*timePerFrame,srcoord,'ob','MarkerFaceColor','b','MarkerSize',10,'DisplayName','Control');            
    haxes = get(h,'CurrentAxes');
    set(haxes,'XLim',[timeShifts(1),timeShifts(end)].*timePerFrame);
    set(haxes,'XTick',(timeShifts(1):5:timeShifts(end))*timePerFrame);
    set(haxes,'XTickLabel',(timeShifts(1):5:timeShifts(end))*timePerFrame);
    set(haxes,'YLim',[-1,1]);
    set(haxes,'YTick',-1:0.5:1);
    set(haxes,'YTickLabel',-1:0.5:1);        
    set(haxes,'FontSize',32);                
    associationFname = [mainDirname 'association/StrainRateCoordination_' uniqueGenes{gene} '.bmp'];
    eval(sprintf('print -dbmp16m %s', associationFname));    
    hold off;
    
    close all;
end


end

%%
function [genesStr] = clusterToString(genesCellArray)
genesStr = genesCellArray{1};
for i = 2 : length(genesCellArray)
    genesStr = [genesStr ', ' genesCellArray{i}];
end
end


%%
function [] = geneAnalysis(out,strLabels,metaData,mainDirname)
geneAnalysisProperty(out.speed,strLabels,metaData,mainDirname,'Speed',out.axesSpeed);
geneAnalysisProperty(out.directional,strLabels,metaData,mainDirname,'Directionality',out.axesDirectionality);
% geneAnalysisProperty(out.strainRate,strLabels,metaData,mainDirname,'StrainRate',out.axesStrainRate);
geneAnalysisProperty(out.coordination,strLabels,metaData,mainDirname,'Coordination',out.axesCoordination);
close all;
end


% TODO: upadte to show NT distinctly from pSup (now they are combined)
function [] = geneAnalysisProperty(data,strLabels,metaData,mainDirname,propertyStr,newAxes)
close all;
fontsize = 24;

uniqueGenes = whGetUniqueGeneLabels(strLabels);
N = length(uniqueGenes);

strPSup = 'pSuper';
strNT = 'NT';
indsPSup = strcmp(whToPrefix(strLabels),strPSup) | strcmp(whToPrefix(strLabels),strNT); % then filter by days for the gene

markersPerm ='sdh';

nDays = length(metaData.groupsByDays);

for day = 1 : nDays
    indsDay = false(1,length(indsPSup));
    indsDay(metaData.groupsByDays{day}.inds) = true;
    strDay = metaData.groupsByDays{day}.dates;
    for gene = 1 : N
        
        if strcmp(uniqueGenes{gene},strPSup) || strcmp(uniqueGenes{gene},strNT)
            continue;
        end
        
        indsGene = strcmp(whToPrefix(strLabels),uniqueGenes{gene});
        %         indsDates = getDatesInds(metaData,indsGene);
        
        indsPSup1 = indsPSup & indsDay;
        indsGene = indsGene & indsDay;
        
        if sum(indsPSup1) < 2
            continue;
        end
        
        if sum(indsGene) < 3
            continue;
        end
        
        [seqStr, seqInds] = whGetSequencesStr(strLabels,uniqueGenes{gene});
        
        h = figure;
        xlabel('1st PC','FontSize',fontsize);
        ylabel('2nd PC','FontSize',fontsize);
        hold on;
        
        plot(data.score(indsPSup1,1)',data.score(indsPSup1,2)','o','MarkerEdgeColor','k','LineWidth',2,'MarkerSize',8,'DisplayName',strPSup);
        
        shSerialNumber = 0;
        for seq = 1 : length(seqStr)
            seqIndsDay = seqInds{seq};
            seqIndsDay = seqIndsDay & indsGene;
            if sum(seqIndsDay) == 0
                continue;
            end
            shSerialNumber = shSerialNumber + 1;
            depletionRateStr = num2str(metaData.KD{find(seqIndsDay,1)}); % to debug
            seqStrStr = seqStr{seq}; seqStrStr = seqStrStr(2:end-1);
            legendStr = [uniqueGenes{gene} '_{' seqStrStr ' (' depletionRateStr ')}'];            
            plot(data.score(seqIndsDay,1)',data.score(seqIndsDay,2)',sprintf('%s',markersPerm(shSerialNumber)),'MarkerEdgeColor','b','LineWidth',2,'MarkerSize',8,'DisplayName',legendStr);
        end
        
        haxes = get(h,'CurrentAxes');
        axis(newAxes);
        set(haxes,'FontSize',fontsize);
        legend('Location','EastOutside');
        set(h,'Color','none');
        %         position = get(h,'position');
        %         set(h,'position',[position(1:2) round(1.2*position(3:4))]);
        %
        pcaFname = [mainDirname 'pca/gene/' propertyStr 'PCA_Gene_' uniqueGenes{gene} '_' strDay '_shSeq.eps'];
        export_fig(pcaFname);
        
        hold off;
        close all;
    end
end
end

function [indsDates] = getDatesInds(metaData,indsGene)
N = length(indsGene);
inds = find(indsGene);
n = length(inds);

indsDates = false(1,N);
for i = 1 : n
    date = num2str(metaData.dates{inds(i)});
    tmp = false(1,N);
    for j = 1 : N
        if strcmp(date,num2str(metaData.dates{j}))
            tmp(j) = true;
        end
    end
    indsDates = indsDates | tmp;    
end
end

%% Day analysis control

function [] = dayAnalysisControl(out,strLabels,metaData,mainDirname)
dayAnalysisControlProperty(out.speed,strLabels,metaData,mainDirname,'Speed',out.axesSpeed);
dayAnalysisControlProperty(out.directional,strLabels,metaData,mainDirname,'Directionality',out.axesDirectionality);
dayAnalysisControlProperty(out.strainRate,strLabels,metaData,mainDirname,'StrainRate',out.axesStrainRate);
dayAnalysisControlProperty(out.coordination,strLabels,metaData,mainDirname,'Coordination',out.axesCoordination);
close all;
end


% TODO: deal with NT independent of pSup
function [] = dayAnalysisControlProperty(data,strLabels,metaData,mainDirname,propertyStr,newAxes)
close all;

fontsize = 24;

N = length(metaData.groupsByDays);
n = length(strLabels);

strPSup = 'pSuper';
strNT = 'NT';
indsPSup = strcmp(whToPrefix(strLabels),strPSup) | strcmp(whToPrefix(strLabels),strNT); % then filter by days for the gene

for day = 1 : N
    dayStr = metaData.groupsByDays{day}.dates;
    daysInds = metaData.groupsByDays{day}.inds;
    indsDay = false(1,n);
    indsDay(daysInds) = true;
    indsPSup1 = indsPSup & indsDay;
    
    if sum(indsPSup1) < 3
        continue;
    end        
    
    
    % Now add the gene for that day (make sure there is only 1 of those)
    indsGenes = indsDay & ~indsPSup;
    
    if sum(indsGenes) < 4
        continue;
    end
    
    genesStr = metaData.treatment(indsGenes);
    
    geneStr = unique(whToPrefix(genesStr));
    
    for igene = 1 : length(geneStr)
        curGeneStr = geneStr{igene};
        indsGene = strcmp(whToPrefix(strLabels),curGeneStr);
        indsGene = indsGene & indsDay;
        
        h = figure;
        xlabel('PC 1','FontSize',fontsize);
        ylabel('PC 2','FontSize',fontsize);
        hold on;
        plot(data.score(indsPSup1,1)',data.score(indsPSup1,2)','ok','MarkerEdgeColor','k','LineWidth',3,'MarkerSize',10,'DisplayName',strPSup);        
        plot(data.score(indsGene,1)',data.score(indsGene,2)','ob','MarkerEdgeColor','b','LineWidth',3,'MarkerSize',10,'DisplayName',curGeneStr);
        legend show;        
    
        
        haxes = get(h,'CurrentAxes');
        axis(newAxes);
        set(haxes,'FontSize',fontsize);
        set(h,'Color','none');                        
               
        pcaFname = [mainDirname 'pca/dayGene/' propertyStr 'PCA_Gene_' curGeneStr '_Day_' dayStr '.eps'];
        export_fig(pcaFname);
        
        hold off;
    end
              
    close all;        
end
end


%% Day analysis gene - compares gene to that day's control

% function [] = dayAnalysisGene(out,strLabels,metaData,mainDirname)
% dayAnalysisGeneProperty(out.speed,strLabels,metaData,mainDirname,'Speed',out.axesSpeed);
% dayAnalysisGeneProperty(out.directional,strLabels,metaData,mainDirname,'Directionality',out.axesDirectionality);
% dayAnalysisGeneProperty(out.strainRate,strLabels,metaData,mainDirname,'StrainRate',out.axesStrainRate);
% dayAnalysisGeneProperty(out.coordination,strLabels,metaData,mainDirname,'Coordination',out.axesCoordination);
% close all;
% end


% TODO: update by day
% function [] = dayAnalysisGeneProperty(data,strLabels,metaData,mainDirname,propertyStr,newAxes)
% close all;
% 
% N = length(metaData.groupsByDays);
% n = length(strLabels);
% 
% strPSup = 'pSuper';
% indsPSup = strcmp(whToPrefix(strLabels),strPSup); % then filter by days for the gene
% 
% for day = 1 : N
%     dayStr = metaData.groupsByDays{day}.dates;
%     daysInds = metaData.groupsByDays{day}.inds;
%     indsDay = false(1,n);
%     indsDay(daysInds) = true;
%     indsPSup1 = indsPSup & indsDay;
%     
%     if sum(indsPSup1) < 2
%         continue;
%     end
%     
%     indsGenes = indsDay & ~indsPSup;
%     
%     if sum(indsGenes) < 4
%         continue;
%     end
%     
%     genesStr = metaData.treatment(indsGenes);
%     
%     geneStr = unique(whToPrefix(genesStr));        
%     
%     nGene = length(geneStr);
%     
%     for curGene = 1 : nGene
%         
%     end
%     
%     h = figure;
%     xlabel('1st PC','FontSize',16);
%     ylabel('2nd PC','FontSize',16);
%     hold all;
%     
%     
%     plot(data.score(indsPSup1,1)',data.score(indsPSup1,2)','ob','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',8,'DisplayName',strPSup);
%     
%         
%     haxes = get(h,'CurrentAxes');
%     axis(newAxes);
%     set(haxes,'FontSize',16);
%     
%     legend show;
%     set(h,'Color','none');
%     position = get(h,'position');
%     set(h,'position',[position(1:2) round(1.2*position(3:4))]);    
%     hold off;
%     
%     pcaFname = [mainDirname 'pca/day/' propertyStr 'PCA_Day_' dayStr '.jpg'];
%     eval(sprintf('print -djpeg %s', pcaFname));
%     
%     % Now add the gene for that day (make sure there is only 1 of those)
%     indsGenes = indsDay & ~indsPSup;
%     
%     if sum(indsGenes) == 0
%         continue;
%     end
%     
%     genesStr = metaData.treatment(indsGenes);
%     
%     geneStr = unique(whToPrefix(genesStr));
%     
%     colors = 'mcky';
%     for igene = 1 : length(geneStr)
%         curGeneStr = geneStr{igene};
%         indsGene = strcmp(whToPrefix(strLabels),curGeneStr);
%         indsGene = indsGene & indsDay;
%         
%         hold all;
%         legend off;
%         % find the gene and plot it
%         hg = plot(data.score(indsGene,1)',data.score(indsGene,2)',sprintf('o%s',colors(igene)),'MarkerFaceColor',sprintf('%s',colors(igene)),'MarkerEdgeColor','k','MarkerSize',12,'DisplayName',curGeneStr);
%         legend show;
%         hold off;
%     
%         pcaFname = [mainDirname 'pca/dayGene/' propertyStr 'PCA_Day_' dayStr '_Gene_' curGeneStr '.jpg'];
%         eval(sprintf('print -djpeg %s', pcaFname));
%         set(hg,'Visible','off');
%     end
%     
%     close all;        
% end
% end

%% Healing rate analysis
function [] = healingRateAnalysis(data,metaData,mainDirname) % data.speed, data.healingRate
h = figure;
xlabel('Healing rate (\mum hour{-1})','FontSize',28);
ylabel('PC1','FontSize',28);
hold all;
for i = 1 : metaData.N
    plot(data.healingRate(i),data.speed.score(i,1),'ok','MarkerFaceColor','r','MarkerSize',8);
end
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',28);
set(h,'Color','none');
% position = get(h,'position');
% set(h,'position',[position(1:2) round(1.2*position(3:4))]);
hold off;
healingRatePCA1Fname = [mainDirname 'healingRate/analysis/healingRateVsSpeedPC1.eps'];
export_fig(healingRatePCA1Fname);

h = figure;
xlabel('Healing rate (\mum hour{-1})','FontSize',28);
ylabel('PC2','FontSize',28);
hold all;
for i = 1 : metaData.N
    plot(data.healingRate(i),data.speed.score(i,2),'ok','MarkerFaceColor','r','MarkerSize',8);
end
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',28);
set(h,'Color','none');
% position = get(h,'position');
% set(h,'position',[position(1:2) round(1.2*position(3:4))]);
hold off;
healingRatePCA1Fname = [mainDirname 'healingRate/analysis/healingRateVsSpeedPC2.eps'];
export_fig(healingRatePCA1Fname);

h = figure;
xlabel('Healing rate (\mum hour{-1})','FontSize',28);
ylabel('PC3','FontSize',28);
hold all;
for i = 1 : metaData.N
    plot(data.healingRate(i),data.speed.score(i,3),'ok','MarkerFaceColor','r','MarkerSize',8);
end
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',28);
set(h,'Color','none');
% position = get(h,'position');
% set(h,'position',[position(1:2) round(1.2*position(3:4))]);
hold off;
healingRatePCA3Fname = [mainDirname 'healingRate/analysis/healingRateVsSpeedPC3.eps'];
export_fig(healingRatePCA3Fname);


% Directionality
h = figure;
xlabel('Healing rate (\mum hour{-1})','FontSize',28);
ylabel('PC1','FontSize',28);
hold all;
for i = 1 : metaData.N
    plot(data.healingRate(i),data.directional.score(i,1),'ok','MarkerFaceColor','r','MarkerSize',8);
end
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',28);
set(h,'Color','none');
% position = get(h,'position');
% set(h,'position',[position(1:2) round(1.2*position(3:4))]);
hold off;
healingRatePCA1Fname = [mainDirname 'healingRate/analysis/healingRateVsDirectionalityPC1.eps'];
export_fig(healingRatePCA1Fname);

h = figure;
xlabel('Healing rate (\mum hour{-1})','FontSize',28);
ylabel('PC2','FontSize',28);
hold all;
for i = 1 : metaData.N
    plot(data.healingRate(i),data.directional.score(i,2),'ok','MarkerFaceColor','r','MarkerSize',8);
end
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',28);
set(h,'Color','none');
% position = get(h,'position');
% set(h,'position',[position(1:2) round(1.2*position(3:4))]);
hold off;
healingRatePCA2Fname = [mainDirname 'healingRate/analysis/healingRateVsDirectionalityPC2.eps'];
export_fig(healingRatePCA2Fname);

h = figure;
xlabel('Healing rate (\mum hour{-1})','FontSize',28);
ylabel('PC3','FontSize',28);
hold all;
for i = 1 : metaData.N
    plot(data.healingRate(i),data.directional.score(i,3),'ok','MarkerFaceColor','r','MarkerSize',8);
end
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',28);
set(h,'Color','none');
% position = get(h,'position');
% set(h,'position',[position(1:2) round(1.2*position(3:4))]);
hold off;
healingRatePCA3Fname = [mainDirname 'healingRate/analysis/healingRateVsDirectionalityPC3.eps'];
export_fig(healingRatePCA3Fname);

% Coordination
h = figure;
xlabel('Healing rate (\mum hour{-1})','FontSize',28);
ylabel('PC1','FontSize',28);
hold all;
for i = 1 : metaData.N
    plot(data.healingRate(i),data.coordination.score(i,1),'ok','MarkerFaceColor','r','MarkerSize',8);
end
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',28);
set(h,'Color','none');
% position = get(h,'position');
% set(h,'position',[position(1:2) round(1.2*position(3:4))]);
hold off;
healingRatePCA1Fname = [mainDirname 'healingRate/analysis/healingRateVsCoordinationPC1.eps'];
export_fig(healingRatePCA1Fname);

h = figure;
xlabel('Healing rate (\mum hour{-1})','FontSize',28);
ylabel('PC2','FontSize',28);
hold all;
for i = 1 : metaData.N
    plot(data.healingRate(i),data.coordination.score(i,2),'ok','MarkerFaceColor','r','MarkerSize',8);
end
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',28);
set(h,'Color','none');
% position = get(h,'position');
% set(h,'position',[position(1:2) round(1.2*position(3:4))]);
hold off;
healingRatePCA1Fname = [mainDirname 'healingRate/analysis/healingRateVsCoordinationPC2.eps'];
export_fig(healingRatePCA1Fname);

h = figure;
xlabel('Healing rate (\mum hour{-1})','FontSize',28);
ylabel('PC3','FontSize',28);
hold all;
for i = 1 : metaData.N
    plot(data.healingRate(i),data.coordination.score(i,3),'ok','MarkerFaceColor','r','MarkerSize',8);
end
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',28);
set(h,'Color','none');
% position = get(h,'position');
% set(h,'position',[position(1:2) round(1.2*position(3:4))]);
hold off;
healingRatePCA3Fname = [mainDirname 'healingRate/analysis/healingRateVsCoordinationPC3.eps'];
export_fig(healingRatePCA3Fname);

% Distribution
[nelements,centers] = hist(data.healingRate,20); 
healingRateDistribution = nelements./sum(nelements);
plotAlignmentAngleDistribution(healingRateDistribution,centers,[mainDirname 'healingRate/analysis/healingRateDistribution.eps']);
end

function [] = plotAlignmentAngleDistribution(distribution,centers,outFname)
h = figure;
hold on;
bar(centers,distribution,'r');
xlabel('Healing rate (\mum hour{-1})','FontSize',22);
ylabel('Percent','FontSize',22);
haxes = get(h,'CurrentAxes');
% set(haxes,'XLim',[0,90]);
% set(haxes,'XTick',0:45:90);
% set(haxes,'XTickLabel',0:45:90);
% set(haxes,'YLim',[0,0.4]);
% set(haxes,'YTick',0:0.1:0.4);
% set(haxes,'YTickLabel',0:0.1:0.4);
set(haxes,'FontSize',22);
set(h,'Color','none');
hold off;
export_fig(outFname);
end



%%
function [] = RhoGTPasesAnalysis(out,metaData,mainDirname)
nDays = length(metaData.groupsByDays);
inds = [];
for i = 1 : nDays 
    if strcmp(metaData.groupsByDays{i}.dates,'20150401') || strcmp(metaData.groupsByDays{i}.dates,'20150402')
        inds = [inds metaData.groupsByDays{i}.inds];
    end
end
treatments = metaData.treatment(inds);
uniqueTreatments = unique(treatments);
nTreats = length(uniqueTreatments);
indsTreats = cell(1,nTreats);
for i = 1 : nTreats
    indsTreats{i} = inds(strcmp(uniqueTreatments{i},treatments));
end

plotRhoGTPasesPCA(out.speed,uniqueTreatments,indsTreats,mainDirname,'Speed');
plotRhoGTPasesPCA(out.directional,uniqueTreatments,indsTreats,mainDirname,'Directionality');
plotRhoGTPasesPCA(out.coordination,uniqueTreatments,indsTreats,mainDirname,'Coordination');
end

function [] = plotRhoGTPasesPCA(data,treatsStr,treatsInds,dirname,measureStr)

close all;

nTreats = length(treatsStr);
colorsPerm = 'mkbcr';

%% PC1 vs. PC2
h = figure;
xlabel('PC1','FontSize',15);
ylabel('PC2','FontSize',15);
hold all;
for i = 1 : nTreats
    inds = treatsInds{i};    
    plot(data.score(inds,1)',data.score(inds,2)','ok','MarkerFaceColor',sprintf('%s',colorsPerm(i)),'MarkerSize',10,...
        'DisplayName',treatsStr{i});    
end
legend('Location','southwest');
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',32);

xlim = get(haxes,'XLim');
ylim = get(haxes,'YLim');
newAxes = [xlim, ylim];

pcaFname = [dirname measureStr 'PCA_legend.jpg'];
eval(sprintf('print -djpeg %s', pcaFname));

legend off;

hold off;

pcaFname = [dirname measureStr 'PCA12.jpg'];
eval(sprintf('print -djpeg %s', pcaFname));

%% PC1 vs. PC3
h = figure;
xlabel('PC1','FontSize',15);
ylabel('PC3','FontSize',15);
hold all;
for i = 1 : nTreats
    inds = treatsInds{i};    
    plot(data.score(inds,1)',data.score(inds,3)','ok','MarkerFaceColor',sprintf('%s',colorsPerm(i)),'MarkerSize',10,...
        'DisplayName',treatsStr{i});    
end

haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',32);

xlim = get(haxes,'XLim');
ylim = get(haxes,'YLim');
newAxes = [xlim, ylim];


hold off;

pcaFname = [dirname measureStr 'PCA13.jpg'];
eval(sprintf('print -djpeg %s', pcaFname));
%%

cdc42DiffScore =  data.score(treatsInds{1},:) - repmat(mean(data.score(treatsInds{5},:)),size(data.score(treatsInds{1},:),1),1);
rac1DiffScore = data.score(treatsInds{3},:) - repmat(mean(data.score(treatsInds{2},:)),size(data.score(treatsInds{3},:),1),1);
rhoaDiffScore = data.score(treatsInds{4},:) - repmat(mean(data.score(treatsInds{2},:)),size(data.score(treatsInds{4},:),1),1);

%% PC1 vs. PC2
h = figure;
xlabel('PC1 (KD - control)','FontSize',15);
ylabel('PC2 (KD - control)','FontSize',15);
hold all;

plot(cdc42DiffScore(:,1)',cdc42DiffScore(:,2)','ok','MarkerFaceColor',sprintf('%s',colorsPerm(1)),'MarkerSize',10,'DisplayName',treatsStr{1});
plot(rac1DiffScore(:,1)',rac1DiffScore(:,2)','ok','MarkerFaceColor',sprintf('%s',colorsPerm(3)),'MarkerSize',10,'DisplayName',treatsStr{3});
plot(rhoaDiffScore(:,1)',rhoaDiffScore(:,2)','ok','MarkerFaceColor',sprintf('%s',colorsPerm(4)),'MarkerSize',10,'DisplayName',treatsStr{4});

legend('Location','southwest');
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',32);

xlim = get(haxes,'XLim');
ylim = get(haxes,'YLim');
newAxes = [xlim, ylim];

pcaFname = [dirname measureStr 'PCADiff_legend.jpg'];
eval(sprintf('print -djpeg %s', pcaFname));

legend off;

hold off;

pcaFname = [dirname measureStr 'PCADiff12.jpg'];
eval(sprintf('print -djpeg %s', pcaFname));

%% PC1 vs. PC3
h = figure;
xlabel('PC1 (KD - control)','FontSize',15);
ylabel('PC3 (KD - control)','FontSize',15);
hold all;
plot(cdc42DiffScore(:,1)',cdc42DiffScore(:,3)','ok','MarkerFaceColor',sprintf('%s',colorsPerm(1)),'MarkerSize',10,'DisplayName',treatsStr{1});
plot(rac1DiffScore(:,1)',rac1DiffScore(:,3)','ok','MarkerFaceColor',sprintf('%s',colorsPerm(3)),'MarkerSize',10,'DisplayName',treatsStr{3});
plot(rhoaDiffScore(:,1)',rhoaDiffScore(:,3)','ok','MarkerFaceColor',sprintf('%s',colorsPerm(4)),'MarkerSize',10,'DisplayName',treatsStr{4});

haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',32);

xlim = get(haxes,'XLim');
ylim = get(haxes,'YLim');
newAxes = [xlim, ylim];


hold off;

pcaFname = [dirname measureStr 'PCA13Diff.jpg'];
eval(sprintf('print -djpeg %s', pcaFname));

end