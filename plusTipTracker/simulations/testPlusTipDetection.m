% Script which simulates EB1 movies of varying shapes and signal-to-noise
% ratios, runs the chosen plusTip detection from the plusTipTracker,
% determines the number of false positives/negatives and create contour
% maps.
%
%
% Sebastien Besson, June 2011

clear
clc
close all

% Set up main analysis folder
dataFolder = '/home/sb286/Desktop/EB1Simulations/Detection';
detectionMethod =2;
generateMovie =0;

% Constant parameters
imSize=512;
roiYX = [1 1;1 imSize;imSize imSize; imSize 1];
pixelSize=.1;
mtDensity = .001;
ampAboveBG =5000;
samplingRate = 1;
totalTime = 40;
sigmay=1.5;
% MT dynamics parameters
mtSimParam = struct('growthSpeed',[20 5],'shrinkageSpeed',[30
5],'growthTime',[20 4],'shrinkageTime',[10 4]);

detectionFcn{1}=@(x)plusTipCometDetector(x,[1 totalTime],16,1);
detectionFcn{2}=@(x)plusTipGaussianCometDetector(x,sigmay,[1 totalTime],16,1,...
    'alpha',.01);

% Parameters sweep
n = 5; % Number of parameters
SNR_min=2;
SNR_max=10;
ecc_min = 1;
ecc_max = 5;
SNR = @(x) SNR_min+(x-1)*(SNR_max-SNR_min)/(n-1);
ecc = @(x) ecc_min+(x-1)*(ecc_max-ecc_min)/(n-1);
ebCometParam = @(x) [sigmay*ecc(x)  sigmay];

% Structure for saving image
saveDir =@(x,y) [dataFolder filesep 'SNR-' num2str(SNR(x)) '_ecc-'...
    num2str(ecc(y))];

%% Create project structure

projList(n)=struct('imDir','','anDir','');
% Initialize image and analysid directories
for i=1:n
    for j=1:n
        projList(i,j).imDir = [saveDir(i,j) filesep 'images'];
        projList(i,j).anDir = [saveDir(i,j) filesep 'roi'];
    end
end

% Create simulation anonuymous functino
saveStruct =@(i,j)struct('filenameBase','EB1_simulations',...
    'dir2save',projList(i,j).imDir);
paramSimulation = @(x,y) simEB1Images(imSize,pixelSize,mtDensity,SNR(x),...
    ampAboveBG,samplingRate,totalTime,mtSimParam,ebCometParam(y),...
    saveStruct(x,y));

%% Initialize simulations
mkClrDir(dataFolder);
for index=1:n^2
    [i,j]=ind2sub([n,n],index);
    % Clears the project directories if any
    mkClrDir(projList(i,j).imDir);
    mkClrDir(projList(i,j).anDir);
    % Save roi containing the whole image
    save([projList(i,j).anDir filesep 'roiYX.mat'],'roiYX')

end

% Generate the series of images
movieInfoGT = cell(n,n);
tracksGT = cell(n,n);
progressText(0,'Generating simulations');
for index=1:n^2
    [i,j]=ind2sub([n,n],index);
    [movieInfoGT{i,j},tracksGT{i,j}]=paramSimulation(i,j);
    progressText(index/n^2);
    
end
save([dataFolder filesep 'simData.mat'],'movieInfoGT','tracksGT');

%% Detect comets
progressText(0,'Detecting comets');
for index=1:n^2
    [i,j]=ind2sub([n,n],index);
    detectionFcn{detectionMethod}(projList(i,j));
    progressText(index/n^2);
end

%% Post-detection analysis
% Detect false positives
s=load([dataFolder filesep 'simData.mat'],'movieInfoGT','tracksGT');
movieInfoGT=s.movieInfoGT;

% Initialize analysis arrays
movieInfoPT = cell(n,n);
meanFN=zeros(n,n);
stdFN=zeros(n,n);
meanFP=zeros(n,n);
stdFP=zeros(n,n);

figure;
analysisFolder = [dataFolder filesep 'analysis'];
mkClrDir(analysisFolder);
% Set the threshold for particle matching
threshold=2;

progressText(0,'Calculating false positives/negatives');
for index=1:n^2
    [i,j]=ind2sub([n,n],index);
    % Load the detection movie Info
    s=load([projList(i,j).anDir filesep 'feat' filesep 'movieInfo.mat']);
    movieInfoPT{i,j} = s.movieInfo;
    
    % Initialize false negative/positive arrays
    fN=zeros(1,totalTime);
    fP=zeros(1,totalTime);
    % Run thorugh the movie and calculate the false positives/negatives
    for t=1:totalTime
        trueComets =[movieInfoGT{i,j}(t).xCoord(:,1) movieInfoGT{i,j}(t).yCoord(:,1)];
        detectedComets =[movieInfoPT{i,j}(t).xCoord(:,1) movieInfoPT{i,j}(t).yCoord(:,1)];
        D1=KDTreeBallQuery(detectedComets,trueComets,threshold*ones(size(trueComets,1),1));
        fN(t) = sum(cellfun(@isempty,D1))/size(trueComets,1);
        D2=KDTreeBallQuery(trueComets,detectedComets,threshold*ones(size(detectedComets,1),1));
        fP(t) = sum(cellfun(@isempty,D2))/size(detectedComets,1);
        
        % Show and save the first image with the false positives and false
        % negatives as overlays
        if t==1
            cla;
            imshow(imread([projList(i,j).imDir filesep 'EB1_simulations_001.tif']),[]);
            hold on
            plot(trueComets(:,1),trueComets(:,2),'ob','MarkerFaceColor','b','MarkerSize',4);
            plot(detectedComets(:,1),detectedComets(:,2),'or','MarkerFaceColor','r','MarkerSize',4);
            plot(trueComets(cellfun(@isempty,D1),1),trueComets(cellfun(@isempty,D1),2),'oc','MarkerSize',14);
            plot(detectedComets(cellfun(@isempty,D2),1),detectedComets(cellfun(@isempty,D2),2),'oy','MarkerSize',14);
            print('-dpng', '-r300', [analysisFolder filesep ...
                'SNR-' num2str(SNR(i)) '_ecc-' num2str(ecc(j)) 'png']);
        end
    end
    
    % Calculate the mean and standard deviation of the false
    % positives/negatives
    meanFN(i,j)=mean(fN);
    stdFN(i,j)=std(fN);
    meanFP(i,j)=mean(fP);
    stdFP(i,j)=std(fP);

    progressText(index/n^2);
end
save([analysisFolder filesep 'analysisData.mat'],'movieInfoPT','meanFN',...
    'meanFP','stdFN','stdFP');
%% Plot results

% define small and large fonts
tfont = {'FontName', 'Helvetica', 'FontSize', 14, 'FontAngle', 'italic'};
sfont = {'FontName', 'Helvetica', 'FontSize', 18};
lfont = {'FontName', 'Helvetica', 'FontSize', 22};

% plot

detectionOutput.var = {meanFN,meanFP};
detectionOutput.name = {'Detection false negatives','Detection false positives'};

for i=1:numel(detectionOutput.var);
    
    figure('PaperPositionMode', 'auto','Position',[50 50 500 500],...
        'Name',detectionOutput.name{i}); % enable resizing
    hold on;
    [C,h] = contour(SNR(1:n),ecc(1:n),detectionOutput.var{i}','LineWidth',2);    
    hText = clabel(C,h);
    set(hText,tfont{:})
    axis square
    
    % Set thickness of axes, ticks and assign tick labels
    box on
    set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
    xlabel('SNR', lfont{:});
    ylabel('Comet shape', lfont{:});
    set(gca,'LooseInset',get(gca,'TightInset'),'XScale','log')
    print('-dpng', '-r300', [analysisFolder filesep detectionOutput.name{i} '.png']);
end

%% Generate detection movie

if generateMovie, return; end
figure;
% Choose detection parameters
i = 5;
j = 4;
s=load([projList(i,j).anDir filesep 'feat' filesep 'movieInfo.mat']);
movieInfoPT{i,j} = s.movieInfo;
figure;
imageFileNames = imDir(projList(i,j).imDir);
MakeQTMovie('start',[dataFolder filesep 'SNR-' num2str(SNR(i)) ...
    '_ecc-' num2str(ecc(j)) '.mov']);
MakeQTMovie('quality',.9)
for t=1:totalTime
    trueComets =[movieInfoGT(t).xCoord(:,1) movieInfoGT(t).yCoord(:,1)];
    detectedComets =[movieInfoPT{i,j}(t).xCoord(:,1) movieInfoPT{i,j}(t).yCoord(:,1)];
    D1=KDTreeBallQuery(detectedComets,trueComets,threshold*ones(size(trueComets,1),1));
    D2=KDTreeBallQuery(trueComets,detectedComets,threshold*ones(size(detectedComets,1),1));
    
    cla;
    imshow(imread([projList(i,j).imDir filesep imageFileNames(t).name]),[]);
    hold on
    plot(trueComets(:,1),trueComets(:,2),'ob','MarkerFaceColor','b','MarkerSize',4);
    plot(detectedComets(:,1),detectedComets(:,2),'or','MarkerFaceColor','r','MarkerSize',4);
    plot(trueComets(cellfun(@isempty,D1),1),trueComets(cellfun(@isempty,D1),2),'oc','MarkerSize',16);
    plot(detectedComets(cellfun(@isempty,D2),1),detectedComets(cellfun(@isempty,D2),2),'oy','MarkerSize',16);
    MakeQTMovie('addfigure')
end
MakeQTMovie('finish')
