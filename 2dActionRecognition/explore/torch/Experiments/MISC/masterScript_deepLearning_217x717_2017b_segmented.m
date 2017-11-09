
% simple script to pre-process the cells for Deep learning.
% Andrew R. Jamieson Oct 2017

% load data state 
clear;
clc;


%% define network
layers = [ ...
    imageInputLayer([217 217 1])
    convolution2dLayer(11, 10, 'Stride', 2,'Name','conv1')
    batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(3,'Stride',2,'Name','maxPool1')
    convolution2dLayer(5, 20,'Stride',2,'Name','conv2')
    batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(3,'Stride',2,'Name','maxPool2')    
    convolution2dLayer(3,50,'Name','conv3')
    reluLayer
    maxPooling2dLayer(3,'Stride',2,'Name','maxPool2')
    reluLayer
    dropoutLayer
%     convolution2dLayer(3,25,'Name','conv4')
%     reluLayer
%     convolution2dLayer(3,25,'Name','conv5','Padding','same')
%     reluLayer
%     maxPooling2dLayer(3,'Stride', 2,'Name','maxPool3','Padding','same')
    fullyConnectedLayer(384,'Name','fc1')
    reluLayer
    dropoutLayer
    fullyConnectedLayer(384,'Name','fc2')
    reluLayer
    dropoutLayer
    fullyConnectedLayer(3,'Name','fc3')
    softmaxLayer
    classificationLayer];


% 2017b
options = trainingOptions('sgdm',...
    'MaxEpochs',100000, ...
    'ValidationFrequency',30,...
    'MiniBatchSize',100,...
    'Verbose',false,...
    'Plots','training-progress',...
    'ValidationData',imdsValid,...
    'CheckpointPath',netTempDir);%fullfile(dataRootDir,'netTEMP'));

%% Training options
%2017a
% options = trainingOptions('sgdm','LearnRateSchedule','piecewise',...
%       'LearnRateDropFactor',0.2,'LearnRateDropPeriod',5,... 
%       'MaxEpochs',1000,'MiniBatchSize',75,...
%       'CheckpointPath',fullfile(dataRootDir,'TEMP_val'));    

trainedNet = trainNetwork(imdsTrain,layers,options);

% %% Resume training a network
% [filename, pathname] = uigetfile('Which file to resume from layers checkpoint?');
% load(fullfile(pathname,filename),'net');
% 
% [filename, pathname] = uigetfile('Which training set to use?');
% load(fullfile(pathname,filename));
% 
% options = trainingOptions('sgdm',...
%     'MaxEpochs',100000, ...
%     'ValidationFrequency',30,...
%     'MiniBatchSize',100,...
%     'Verbose',false,...
%     'Plots','training-progress',...
%     'ValidationData',imdsValid,...
%     'CheckpointPath',netTempDir);%fullfile(dataRootDir,'netTEMP'));
% 
% trainedNet = trainNetwork(imdsTrain,net.Layers,options);
