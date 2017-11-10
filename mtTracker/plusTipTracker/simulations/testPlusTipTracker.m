% Script which simulates EB1 movies of varying densities, generate some
% detection misses, runs the tracking from the plusTipTracker, determines 
% the number of fals posotives/negatives and create contour maps.
%
% Sebastien Besson, June 2011

clear
clc
close all
warning('off','MATLAB:DELETE:DirectoryDeletion');
warning('off','MATLAB:singularMatrix');

% Set up main analysis folder
dataFolder = '/home/sb286/Desktop/EB1Simulations/Tracking';
analysisFolder = [dataFolder filesep 'analysis'];

% Script flads
generateEB1Tracks = 0;
trackComets =0;
analyzeData =1;
plotData=1;


% Constant parameters
nMT=150; % Fix the total number of microtubules
pixelSize=.1;
SNR=10; 
sigmay=1.5;
ecc=2;
ampAboveBG =5000;
samplingRate = 1;
nFrames = 100;

% MT dynamics parameters
meanGrowthSpeed=20;
deltaGrowthSpeed =2;
meanShrinkageSpeed=30;
deltaShrinkageSpeed=3;
mtSimParam = struct('growthSpeed',[meanGrowthSpeed deltaGrowthSpeed],...
    'shrinkageSpeed',[meanShrinkageSpeed deltaShrinkageSpeed],...
    'growthTime',[20 8],'shrinkageTime',[10 4]);

% Varying parameters (handled by anonymous functions)
nParameters = 4;
densityMin=1e-5; % Typical density 1e-4
densityMax=1e-2;
detMissMin=1;
detMissMax=17;
mtDensity = @(x) 10.^(log10(densityMin)+...
    (x-1).*(log10(densityMax)-log10(densityMin))/(nParameters-1));
detMiss = @(x) detMissMin+...
    (x-1).*(detMissMax-detMissMin)/(nParameters-1);
imSize= @(x) round(sqrt(nMT/mtDensity(x)));
roi = @(x)[1 1;1 imSize(x);imSize(x) imSize(x); imSize(x) 1];

% Structure for saving image
saveDir =@(x,y) [dataFolder filesep 'MT-density-' num2str(mtDensity(x)) ...
    '-detectionMisses-' num2str(detMiss(y))];

%% Create project structure

% Initialize image and analysid directories
projList(nParameters,nParameters)=struct('imDir','','anDir','');
for i=1:nParameters
    for j=1:nParameters
        projList(i,j).imDir = [saveDir(i,j) filesep 'images'];
        projList(i,j).anDir = [saveDir(i,j) filesep 'roi_1'];
    end
end

% Create simulation anonymous function
saveStruct =@(i,j)struct('filenameBase','EB1_simulations',...
    'dir2save',[saveDir(i,j) filesep 'images']);
paramSimulation = @(x,y) simEB1Images(imSize(x),pixelSize,mtDensity(x),SNR,...
    ampAboveBG,samplingRate,nFrames,mtSimParam,[sigmay*ecc sigmay],...
    saveStruct(x,y));

%% Initialize simulations
if generateEB1Tracks
    
    
    % Initialize ground trouth detection and tracking results
    % movieInfoGT = cell(1,nParameters);
    tracksGT = cell(1,nParameters);
    tracksSimMiss = cell(1,nParameters);
    detMissesGT = zeros(nParameters,nParameters);
    movieInfoFields = {'xCoord','yCoord','orient','amp'};
    
    hFig=figure;
    progressText(0,'Generating simulations');
    for i=1:nParameters
        [~,tracksGT{1,i}]=paramSimulation(i,1);
        
        %     nCometsMean=mean(arrayfun(@(x) size(x.amp,1),movieInfoGT{1,i}));
        for j=1:nParameters
            featDir = [saveDir(i,j) filesep 'roi_1' filesep 'feat'];
            mkClrDir(featDir);
            % Get the associated movieInfo and randomly remove featuers
            [movieInfo,tracksSimMiss{i,j} detMissesGT(i,j)] = ...
                genMovieInfoFromMtTracks(tracksGT{1,i}, detMiss(j));
            
            gapInfoSimMiss=findTrackGaps(tracksSimMiss{i,j});
            save([featDir filesep 'movieInfo.mat'],'movieInfo');
            
            
            % Show the average number of comets per frame
            figure(hFig)
            plot(arrayfun(@(x) size(x.amp,1),movieInfo));
            xlabel('Frames');
            ylabel('Number of comets');
            print(hFig,'-dpng',[projList(i,j).anDir filesep 'nComets.png'])
            progressText(sub2ind([nParameters,nParameters],j,i)/nParameters^2);
        end
    end
    [mtDensities,detMisses] = meshgrid(mtDensity(1:nParameters),detMiss(1:nParameters));
    save([dataFolder filesep 'simData.mat']);
end
%% Detect comets
% Use default parameters of plusTipTracker
if trackComets
    minTrackLen =2;
    minRadius = 5;
    maxRadius = 15;
    maxFAngle = 30;
    maxBAngle =10;
    maxShrinkFactor = s.meanShrinkageSpeed/s.meanGrowthSpeed;
    fluctRad=1.5;
    
    projList(nParameters,nParameters)=struct('imDir','','anDir','');

    progressText(0,'Tracking comets');
    for index=1:nParameters^2
        [i,j]=ind2sub([nParameters,nParameters],index);
        
        % Create projList struct
        projList.imDir = [saveDir(i,j) filesep 'images'];
        projList.anDir = [saveDir(i,j) filesep 'roi_' detectionMethod(k).name];
        mkClrDir(projList.anDir);
        save([projList.anDir filesep 'roiYX.mat'],'roiYX');
        
    
        gapInfoGT=findTrackGaps(tracksSimMiss{i,j});
        timeWindow = max(gapInfoGT(:,4))+1;
        plusTipCometTracker(projList,timeWindow,minTrackLen,minRadius,...
            maxRadius,maxFAngle,maxBAngle,maxShrinkFactor,fluctRad,[1 nFrames]);
        progressText(index/nParameters^2);
    end
end

%% Post-detection analysis
if analyzeData
    s=load([dataFolder filesep 'simData.mat']);
    
    % Initialize results array
    tracksPT = cell(nParameters,nParameters);
    linkFN=zeros(nParameters,nParameters);
    linkFP=zeros(nParameters,nParameters);
    linkRatio=zeros(nParameters,nParameters);
    nfgapGT=zeros(nParameters,nParameters);
    nbgapGT=zeros(nParameters,nParameters);
    fgapFP=zeros(nParameters,nParameters);
    fgapFN=zeros(nParameters,nParameters);
    bgapFP=zeros(nParameters,nParameters);
    bgapFN=zeros(nParameters,nParameters);
    
    mkClrDir(analysisFolder);
    progressText(0,'Analyzing track links and gaps');
    for index=1:nParameters^2
        % Convert the linear index into (i,j) indices
        [i,j]=ind2sub([nParameters,nParameters],index);
        
        % Load the tracked results
        res=load([saveDir(i,j) filesep 'roi_1' filesep 'track' filesep 'trackResults.mat']);
        tracksPT{i,j} = res.tracksFinal;
        
        % Generate the link and gap stats using the simulated and tracked
        % tracks
        [linkStats,fgapStats,bgapStats,gapType]=mtScoreLinksGaps(tracksPT{i,j},...
            s.tracksSimMiss{i,j});
        
        [linkStats2,gapStats]=scoreLinksGapsMS(tracksPT{i,j},s.tracksSimMiss{i,j});
        
        %Check result consistency)
        assert(isequal(linkStats,linkStats2))
        assert(length(gapStats)==length(fgapStats)+length(bgapStats));
        
        % Calculate the links false negatives, positives and ratios
        linkFN(i,j)=mean(1-linkStats(linkStats(:,1)~=0,3)./linkStats(linkStats(:,1)~=0,1));
        linkFP(i,j)=mean(linkStats(linkStats(:,1)~=0,4)./linkStats(linkStats(:,1)~=0,1));
        linkRatio(i,j) = mean(linkStats(linkStats(:,1)~=0,2)./linkStats(linkStats(:,1)~=0,1));
        
        nbgapGT(i,j) = sum(strcmp(gapType,'b'));
        nfgapGT(i,j) = numel(gapType)-nbgapGT(i,j);
        fgapFP(i,j)=sum(fgapStats(:,2)~=0)/nfgapGT(i,j);
        fgapFN(i,j)=1-sum(fgapStats(:,2)==0)/nfgapGT(i,j);
        
        bgapFP(i,j)=sum(bgapStats(:,2)~=0)/nbgapGT(i,j);
        bgapFN(i,j)=1-sum(bgapStats(:,2)==0)/nbgapGT(i,j);
        
        gapInfoSimMiss=findTrackGaps(s.tracksSimMiss{i,j});
        fgapInfoSimMiss =gapInfoSimMiss(strcmp(gapType,'f'),:);
        gapInfoPT=findTrackGaps(tracksPT{i,j});
        
        [~,ia]=setdiff(fgapInfoSimMiss(:,3:4),gapInfoPT(:,3:4),'rows');
        
        progressText(index/nParameters^2);
    end
    save([analysisFolder filesep 'analysisData.mat'],'linkFP','linkFN',...
        'fgapFP','fgapFN','bgapFP','bgapFN','nfgapGT','nbgapGT','linkRatio');
end

%% Plot results
if plotData
    load([dataFolder filesep 'simData.mat']);
    load([analysisFolder filesep 'analysisData.mat']);

    close all
    % define small and large fonts
    tfont = {'FontName', 'Helvetica', 'FontSize', 14, 'FontAngle', 'italic'};
    sfont = {'FontName', 'Helvetica', 'FontSize', 18};
    lfont = {'FontName', 'Helvetica', 'FontSize', 22};
    
    
    % bgapFN
    gapOutput.var = {linkFN,linkFP,linkRatio,bgapFP,bgapFN,fgapFP,fgapFN};
    gapOutput.name = {'Link false negatives','Link false positives',...
        'Link ratio','Backward gaps false positives',...
        'Backward gaps false negatives','Forward gaps false positives',...
        'Forward gaps false negatives'};
    gapOutput.range={.01:.01:.03,.01:.01:.03,.995:.001:1.01,.2:.2:.8,.2:.2:.8,...
        .02:.02:.1,.1:.1:.6};
    for i=1:numel(gapOutput.var);
        figure('PaperPositionMode', 'auto','Position',[50 50 500 500],...
            'Name',gapOutput.name{i}); % enable resizing
        hold on;
        [C,h] = contour(interp2(s.mtDensities,8),interp2(s.detMisses,8),...
            interp2(gapOutput.var{i},8),gapOutput.range{i},'LineWidth',2);
        hText = clabel(C,h);
        set(hText,tfont{:})
        axis square
        
        % Set thickness of axes, ticks and assign tick labels
        box on
        set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
        xlabel('Microtubule density', lfont{:});
        ylabel('Missed detections', lfont{:});
        set(gca,'LooseInset',get(gca,'TightInset'),'XScale','log')
        print('-dpng', '-r300', [analysisFolder filesep gapOutput.name{i} '.png']);
    end
end
