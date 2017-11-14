
% simple script to pre-process the cells for Deep learning.
% Andrew R. Jamieson Oct 2017


clear;
clc;


%% load data state 
[filename, pathname] = uigetfile('/work/bioinformatics/shared/dope/torch/test/217x217/segmented/','Which data to load?');
load(fullfile(pathname,filename));


%% Define Checkpoint
[checkPointDir] = uigetdir('/work/bioinformatics/shared/dope/torch/test/217x217/segmented/netTEMP',...
                            'Where to store net Checkpoints?');

if ~verLessThan('matlab', '9.3')

    %% define network
    disp('Note: using 2017b MATLAB layer definitions');
    layers = [ ...netTempDir
        imageInputLayer([217 217 1])
        convolution2dLayer(11,24, 'Stride', 4,'Name','conv1')
%         batchNormalizationLayer
        reluLayer
        maxPooling2dLayer(3,'Stride',2,'Name','maxPool1')
        convolution2dLayer(5, 48,'Stride',2,'Name','conv2')
%         batchNormalizationLayer
        reluLayer
        maxPooling2dLayer(3,'Stride',2,'Name','maxPool2')    
        convolution2dLayer(3,36,'Name','conv3')
        reluLayer
        maxPooling2dLayer(3,'Stride',2,'Name','maxPool3')
%         reluLayer
%         convolution2dLayer(3,48,'Name','conv4','Padding','same')
%         reluLayer
%         convolution2dLayer(3,48,'Name','conv5','Padding','same')
%         reluLayer
%         maxPooling2dLayer(3,'Stride', 2,'Name','maxPool4','Padding','same')
        fullyConnectedLayer(256,'Name','fc1')
        reluLayer
        dropoutLayer
        fullyConnectedLayer(256,'Name','fc2')
        reluLayer
        dropoutLayer
        fullyConnectedLayer(3,'Name','fc3')
        softmaxLayer
        classificationLayer]

else
    
    disp('Note: using 2017a MATLAB layer definitions');
    layers = [ ...
        imageInputLayer([217 217 1])
        convolution2dLayer(11, 10, 'Stride', 2,'Name','conv1')
        reluLayer
        maxPooling2dLayer(3,'Stride',2,'Name','maxPool1')
        convolution2dLayer(5, 20,'Stride',2,'Name','conv2')
        reluLayer
        maxPooling2dLayer(3,'Stride',2,'Name','maxPool2')    
        convolution2dLayer(3,50,'Name','conv3')
        reluLayer
        maxPooling2dLayer(3,'Stride',2,'Name','maxPool2')
        reluLayer
        dropoutLayer
        fullyConnectedLayer(384,'Name','fc1')
        reluLayer
        dropoutLayer
        fullyConnectedLayer(384,'Name','fc2')
        reluLayer
        dropoutLayer
        fullyConnectedLayer(3,'Name','fc3')
        softmaxLayer
        classificationLayer]
end

%% Training options
if ~verLessThan('matlab', '9.3')
                        % 2017b
    options = trainingOptions('sgdm',...
        'MaxEpochs',100000, ...
        'ValidationFrequency',100,...
        'MiniBatchSize',100,...
        'Verbose',false,...
        'Plots','training-progress',...
        'ValidationData',imdsValid,...
        'ValidationPatience', Inf,...
        'ExecutionEnvironment' , 'multi-gpu',...
        'CheckpointPath',checkPointDir);
else 
%     2017a
    options = trainingOptions('sgdm','LearnRateSchedule','piecewise',...
          'LearnRateDropFactor',0.2,'LearnRateDropPeriod',5,... 
          'MaxEpochs',1000,'MiniBatchSize',75,...
          'CheckpointPath',checkPointDir);    
end

%% Start Training
trainedNet = trainNetwork(imdsTrain,layers,options);
