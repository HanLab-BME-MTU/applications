% simple script to pre-process the cells for Deep learning.

% Load image
% imNum = 2695;
% MD = load(cellDataSet{imNum}.cellMD);
% I = gpuArray(mat2gray(MD.MD.getChannel(1).loadImage(1)));
% imshow(imresize(I, [64 64]));


load('/work/bioinformatics/shared/dope/data/OMETIFF/Gen2n3_May15_ALL.mat', 'cellDataSet');

dataRootDir = '/work/bioinformatics/shared/dope/torch/test';
data_lowMet = fullfile(dataRootDir, 'lowMet');
data_highMet = fullfile(dataRootDir, 'highMet');
data_unMet = fullfile(dataRootDir, 'unMet');

newDim = [64 64];

parfor i = 1:length(cellDataSet)
    MD = load(cellDataSet{i}.cellMD,'MD');
%     I = gpuArray(mat2gray(MD.MD.getChannel(1).loadImage(1)));
    I = mat2gray(MD.MD.getChannel(1).loadImage(1));
    newI = imresize(I, newDim);
    MD = MD.MD;
    
    newFileOut = [MD.processes_{1}.funParams_.key '.png'];
    
    if MD.processes_{1}.funParams_.metEff == 0
        newFileOut = fullfile(data_lowMet, newFileOut);
    elseif MD.processes_{1}.funParams_.metEff == 1
        newFileOut = fullfile(data_highMet, newFileOut);
    else
        newFileOut = fullfile(data_unMet, newFileOut);
    end
    imwrite(newI, newFileOut);
end
   
imgs =imageDatastore(dataRootDir, 'IncludeSubfolders',true,'FileExtensions','.png','LabelSource','foldernames');
% get rid of unMet...
[imds1] = splitEachLabel(imgs,1500,'Exclude','unMet');
[trainData, valData] = splitEachLabel(imds1, 1000, 'randomize');


layers = [ ...
    imageInputLayer([64 64 1])
    convolution2dLayer(5,20)
    reluLayer
    maxPooling2dLayer(2,'Stride',2)
    fullyConnectedLayer(10)
    softmaxLayer
    classificationLayer]

options = trainingOptions('sgdm',...
    'MaxEpochs',3, ...
    'ValidationData',valData,...
    'ValidationFrequency',30,...
    'Verbose',false,...
    'Plots','training-progress');


% Steps
% Supress background ?

% Segment
% 

% load /work/bioinformatics/shared/dope/data/OMETIFF/Gen2n3_May15_ALL.mat
% load /work/bioinformatics/shared/dope/data/OMETIFF/Gen2n3_May16_metEffOnly.mat
% load '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/All/140604_MV3/140604_MV3_s01/140604_MV3_s01.mat'
% movieViewer(MD)
% ls
% I=MD.load();
% help MD.load
% help MD.getSeries
% I = MD.getSeries();
% MD.getChannel
% MD.getChannel()
% MD.getChannel(1)
% MD.getChannelNames
% MD.getChannel
% MD.getChannel)
% MD.getChannel(1)
% ch=MD.getChannel(1);
% ch.loadImage
% I=ch.loadImage(1);
% imshow(I):
% imshow(I);
% imshow(I,[]);
% I=mat2gray(ch.loadImage(1));
% imshow(I,[]);
% close all
% imshow(I,[]);
% imshowpair(stdfilt(I),I,';