
%% Analysis of the different treatments
function [] = whMetaAnalysis2016(mainDirname,metaDataFname,timePerFrame)

addpath(genpath('/home2/azaritsky/code/applications/monolayer/HallScreen/'));

% mainDirname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ScreenFinal/';
% metaDataFname = 'GEFScreenFinal_kd0.mat';
% timePerFrame = 5;

if nargin == 0
    %% Screen    
    %     mainDirname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ScreenFinal/kd50time200_RHOA/';
    %     metaDataFname = 'GEFScreenFinal20160526_RHOA_kd50.mat';
    
    %     mainDirname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ScreenFinal/kd50time200/';
    %     metaDataFname = 'GEFScreenFinal20160526_kd50.mat';
    
    %     mainDirname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ScreenFinal/kd0time200/';
    %     metaDataFname = 'GEFScreenFinal20160513_kd0.mat';
    
    %% Followups, screen + followups    
    %     mainDirname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20161101/';
    %     metaDataFname = 'GEFProjectAll20160523_kd0.mat'; % this is the real "age" of the data
    
    %     mainDirname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160523/';
    %     metaDataFname = 'GEFProjectAll20160523_kd0.mat';
    %     mainDirname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/20160523/';
    %     metaDataFname = 'YunYuFollowup20160523_AZ_kd0.mat';
    %     mainDirname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160516/';
    %     metaDataFname = 'GEFProjectAll20160516_kd0.mat';
    %     mainDirname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll_Lcf/';
    %     metaDataFname = 'GEFProjectAll20160330_kd0.mat';
    
    %% Revision
    % RhoA day 3
    mainDirname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/Revision201612/';
    metaDataFname = 'ShefaliRevision201612_RhoA_Day3_kd0.mat'; % this is the real "age" of the data
    % RhoC hairpin #3
    %     mainDirname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/Revision201612_RhoC3/';
    %     metaDataFname = 'ShefaliRevision201612_RhoC_Hairpin3_kd0.mat';
    % all
    %     mainDirname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/Revision201612/';
    %     metaDataFname = 'ShefaliRevision201612All_kd0.mat'; % this is the real "age" of the data
    % first experiment
    %     mainDirname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/Revision20161215/';
    %     metaDataFname = 'ShefaliRevision20161215_kd0.mat'; % this is the real "age" of the data
    
    timePerFrame = 5;
end

%% Init
close all;
addpath(genpath('/home2/azaritsky/code/extern/export_fig'));
% addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/utils/notBoxPlot'));/project/bioinformatics/Danuser_lab/shared/assaf/oldUTSWCodeAndData/code/utils/notBoxPlot
addpath(genpath('/home2/azaritsky/code/applications/monolayer'));

load([mainDirname metaDataFname]); % metaData

metaData.timePerFrame = timePerFrame;
metaData.timeToAnalyze = 200;
metaData.spaceToAnalyze = 12; % in patches, 12 x 15um = 180um
metaData.timePartition = 4;
metaData.spatialPartition = 3;
% ??
metaData.initialTimePartition = 0; %2

% metaData.inputDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/allData/'; % Revision20161215 was not moved
metaData.inputDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Data201612_shefali/';
% metaData.inputDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Data20161215_shefali/';
% metaData.inputDname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Data20160523/';

%% RETURN THIS ASSERTION BACK!
% assert(strcmp(metaData.inputDname,'/project/bioinformatics/Danuser_lab/GEFscreen/analysis/allData/'));

[dirs] = whMetaInitDirectories2016(mainDirname);

[labels, strLabels] = getMetaLabels2016(metaData); % labels are combination of gene & shSeq

% always = 0;

%% META ANALYSIS - always extract high-dimensional features from the data and perform PCA (based on the screen's controls)
%   1. KD efficiency
%   2. Printouts of micrscopy error corrections (no output)
%   3. Extract high-dimensional features from the data
%       a. Spatiotemporal associations (time delay) between measures (commented out)
%       b. PCA
%   4. Visualize PCs
% Discarded   5. Gene abalysis (not day as atomic unit)
%   6. Correlate PCs of each measure to wound healing rate (healingRate directory)
%   7. MetaDayAnalysis
% Commented out  8. Shear strain analysis: correlation shear strain events and motion
%       a. CDC42 vs. Control
%       b. All experiments
% Discarded   9. RhoGTPAses analysis

flags.kdEfficiency = 0;
flags.correctMicroscopeError = 0;
flags.visualizePCs = 0;
flags.corrPCsHealingRate = 0;
flags.metaDayAnalysis = 1;
% flags.shearStrainAnalysis = 0;

%% 1. KD efficiency
% Statistics of KD efficiency of all sh-sequences
if flags.kdEfficiency
    whMetaKnockdownEfficiency(metaData,strLabels,dirs.hairpinEfficiencyDir);
end

%% 2. Statistics (printouts) of correction microscope repeat errors
if flags.correctMicroscopeError    
    whMetaMicrscopeRepeatCorrectionStats(metaData,[metaData.inputDname 'correctMotion/']);
end

%% 3. Extract high-dimensional features from the data
%       a. Spatiotemporal associations (time delay) between measures (commented out)
%       b. PCA
allFeatsFname = [mainDirname '/allFeatures.mat'];
reuseFeats = false; 

if exist(allFeatsFname,'file') && reuseFeats
    load(allFeatsFname);
else
    [allFeatures] = getAllFeatureVectors(metaData);
    save(allFeatsFname,'allFeatures');
end

% patch for fix directionality > 8
allFeatures.directionalityFeats.features(allFeatures.directionalityFeats.features > 10) = 10;
save(allFeatsFname,'allFeatures');

% associations = spatioTemporalAssociations(allFeatures,strLabels,mainDirname,timePerFrame);

% TODO: fill empty features? Use machine learning to predict empty spots
% one-by-one based on previous and future spatiotemporal data.

% Also plots gene (no sh-seq data);
plotPCAFname = [mainDirname '/pcaFeatures.mat'];
plotBasicPCAFname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ScreenFinal/pcaParamsScreen.mat';
reuseFeats = true; 

if exist(plotPCAFname,'file') && reuseFeats
    load(plotPCAFname); % out
    load(plotBasicPCAFname); % outControls
else
    [indsControl] = whControlInds(strLabels);
    [out,outControls] = pcaFeatures(allFeatures,indsControl); % this includes the healing rate
            
    save(plotPCAFname,'out');
    if strfind(mainDirname,'ScreenFinal') % strcmp('/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ScreenFinal/kd0time200/',mainDirname)        
        save(plotBasicPCAFname,'outControls');
    else
        warning('not the main screen, using other PCA baseline');
    end
end

%%   4. Visualize PCs
if flags.visualizePCs
    close all;
    visualizePCs(out,metaData,[mainDirname 'pcs/all/']);
    visualizePCs(outControls,metaData,[mainDirname 'pcs/control/']);
end

%%  6. Correlate PCs of each measure to wound healing rate
if flags.corrPCsHealingRate
    % Correlates PC1 to monolayer speed (same for PC2)
    close all;
    %     would not work because maxAbsPC is nto here
    %     healingRateAnalysis(out,metaData,[dirs.healingRateDir 'all/']); %    
    maxAbsPC = 10;
    healingRatePCsAssociation(outControls,[dirs.healingRateDir 'control/'],maxAbsPC);
end

%%  7. Meta day analysis
if flags.metaDayAnalysis
    whMetaDayAnalysis2016(allFeatures,out.healingRate,strLabels,metaData,mainDirname,plotBasicPCAFname);
end

% %%  8. Shear strain analysis: correlation shear strain events and motion
% %       a. CDC42 vs. Control
% %       b. All experiments
% if flags.shearStrainAnalysis
%     try
%         [p,x0,x1] = whShearStrainCntrlCdc42(metaData,mainDirname); % motion-shear strain for all controls + CDC42
%         whMotionShearStrainGeneDaySeq(metaData,mainDirname,p,x0,x1);
%     catch exception
%     end
% end

end




%% 
function [allFeatures] = getAllFeatureVectors(metaData)

mainDirname = metaData.inputDname;

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

if ~strcmp(measureStr,'coordination')
    feats.kymographsStd = cell(1,metaData.N);
else
    feats.kymographsStd = nan;
end

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
    
    % Std Kymographs 
    if ~strcmp(measureStr,'coordination')
        load([kymographDir measureStr 'Std/' metaData.fnames{i} '_' measureStr 'KymographStd.mat']);
        if exist('speedKymographStd','var')
            kymograph = speedKymographStd;
        end
        if exist('directionalityKymographStd','var')
            kymograph = directionalityKymographStd;
        end
        
        feats.kymographsStd{i} = kymograph;
    end
end

[minDistFromWound,allDistancesFromWound] = getDistanceFromWound(feats.kymographs);

% assert(minDistFromWound == 12);
if minDistFromWound < metaData.spaceToAnalyze % 180 um 
    warning('Too few cells in %d movies',sum(allDistancesFromWound<metaData.spaceToAnalyze));    
    minDistFromWound = metaData.spaceToAnalyze;
else
    minDistFromWound = metaData.spaceToAnalyze;
end


[feats.features, feats.kymographs, feats.kymographsStd]= getFeatures(feats.kymographs,feats.kymographsStd,metaData,minDistFromWound);
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
function [features, kymographs,kymographsStd] = getFeatures(kymographs,kymographsStd,metaData,nDist)
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
    
    % HERE THE KYMOGRAPH IS TRIMMED!
    kymographs{i} = curKymograph(1:nDist,1:nTime);
    if iscell(kymographsStd)
        curKymographStd = kymographsStd{i};
        kymographsStd{i} = curKymographStd(1:nDist,1:nTime);
    end
end

end

%% 
function [allHealingRates] = getHealingRate(healingRateDir,metaData)
allHealingRates = zeros(1,metaData.N);
nTime = floor(metaData.timeToAnalyze/metaData.timePerFrame);

for i = 1 : metaData.N
    clear averageHealingRate;
    healingRateFname = [healingRateDir '/' metaData.fnames{i} '_healingRate.mat'];
    % TODO: THIS IS AN UGLY PATCH - CALCULATE HEALING RATES EVEN FOR THESE
    % ONES!
    if ~exist(healingRateFname,'file')
        allHealingRates(i) = nan;
        continue;
    end
    load(healingRateFname);
    allHealingRates(i) = averageHealingRate(nTime);    
end
end

%%
function [out,outControls] = pcaFeatures(allFeatures,indsControl)
out.healingRate = allFeatures.healingRate;
outControls.healingRate = allFeatures.healingRate(:,indsControl);

% todo: switch all coeff and score for data.coeff(1,iPc) to make the first weight positive
out.speed = whGetPCA(allFeatures.speedFeats.features);
out.directional = whGetPCA(allFeatures.directionalityFeats.features);
out.coordination = whGetPCA(allFeatures.coordinationFeats.features);

outControls.speed = whGetPCA(allFeatures.speedFeats.features(:,indsControl));
outControls.directional = whGetPCA(allFeatures.directionalityFeats.features(:,indsControl));
outControls.coordination = whGetPCA(allFeatures.coordinationFeats.features(:,indsControl));

% curFeats = allFeatures.speedFeats.features';
% [nObs mFeats] = size(curFeats);
% meanFeat = mean(curFeats);
% stdFeat = std(curFeats);
% curFeatsNorm = (curFeats - repmat(meanFeat,[nObs 1])) ./ repmat(stdFeat,[nObs 1]); % zscore(A)
% 
% % [coeff D] = eig(cov(curFeatsNorm));
% [coeff1,score1,latent1] = pca(curFeatsNorm);
% score = curFeatsNorm * coeff1;
end

%%
function visualizePCs(out,metaData,outDname)
pcPlot(out.speed,'Speed',metaData,outDname); close all;
pcPlot(out.directional,'Directionality',metaData,outDname); close all;
pcPlot(out.coordination,'Coordination',metaData,outDname); close all;
end

function pcPlot(data,measureStr,metaData,outDname)

if ~exist(outDname,'dir')
    unix(sprintf('mkdir %s',outDname));
end

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
    
    pcaFname = [outDname measureStr '_PC_' num2str(iPc) '.bmp'];       
    
    eval(sprintf('print -dbmp16m  %s', pcaFname));                  
    
    iPc = iPc + 1;
end
end



%% Healing rate analysis
function [] = healingRatePCsAssociation(data,dirname,maxAbsPC) % data.speed, data.healingRate

if ~exist(dirname,'dir')
    mkdir(dirname);
end

plotHealingRatePC(data.healingRate,data.speed.score,'PC1',maxAbsPC,[dirname 'healingRateVsSpeedPC1.eps']);
plotHealingRatePC(data.healingRate,data.speed.score,'PC2',maxAbsPC,[dirname 'healingRateVsSpeedPC2.eps']);
plotHealingRatePC(data.healingRate,data.speed.score,'PC3',maxAbsPC,[dirname 'healingRateVsSpeedPC3.eps']);
plotHealingRatePC(data.healingRate,data.directional.score,'PC1',maxAbsPC,[dirname 'healingRateVsDirectionalityPC1.eps']);
plotHealingRatePC(data.healingRate,data.directional.score,'PC2',maxAbsPC,[dirname 'healingRateVsDirectionalityPC2.eps']);
plotHealingRatePC(data.healingRate,data.directional.score,'PC3',maxAbsPC,[dirname 'healingRateVsDirectionalityPC3.eps']);
plotHealingRatePC(data.healingRate,data.coordination.score,'PC1',maxAbsPC,[dirname 'healingRateVsCoordinationPC1.eps']);
plotHealingRatePC(data.healingRate,data.coordination.score,'PC2',maxAbsPC,[dirname 'healingRateVsCoordinationPC2.eps']);
plotHealingRatePC(data.healingRate,data.coordination.score,'PC3',maxAbsPC,[dirname 'healingRateVsCoordinationPC3.eps']);

% Distribution
[nelements,centers] = hist(data.healingRate,20); 
healingRateDistribution = nelements./sum(nelements);
plotAlignmentAngleDistribution(healingRateDistribution,centers,data.healingRate,[dirname 'healingRateDistribution.eps']);
end

function [] = plotAlignmentAngleDistribution(distribution,centers,healingRates,outFname)
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

save([outFname(1:end-4) '.mat'],'centers','distribution','healingRates');
end

function [] = plotHealingRatePC(healingRate,score,pcLabel,maxAbsPC,outFname)
i = str2double(pcLabel(3));
h = figure;
xlabel('Healing rate (\mum hour{-1})','FontSize',28);
ylabel(pcLabel,'FontSize',28);
hold all;
plot(healingRate,score(:,i),'ok','MarkerSize',3,'MarkerFaceColor','k');%'MarkerFaceColor','r','LineWidth',2
% for i = 1 : length(data.healingRate)
%     plot(data.healingRate(i),data.speed.score(i,1),'ok','MarkerFaceColor','r','MarkerSize',8);
% end
ylim([-maxAbsPC,maxAbsPC]);
haxes = get(h,'CurrentAxes');
set(haxes,'FontSize',28);
set(h,'Color','none');
% position = get(h,'position');
% set(h,'position',[position(1:2) round(1.2*position(3:4))]);
hold off;
export_fig(outFname);

pcData = score(:,i);
save([outFname(1:end-4) '.mat'],'healingRate','pcData','pcLabel','maxAbsPC');

end