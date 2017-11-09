
% simple script to pre-process the cells for Deep learning.
% Andrew R. Jamieson Oct 2017

% load data state 
clear;
clc;

[checkPointDir] = uigetdir('/work/bioinformatics/shared/dope/torch/test/217x217/segmented/netTEMP',...
                            'Where to store net Checkpoints?');

if ~verLessThan('matlab', '9.3')

    %% define network
    disp('Note: using 2017b MATLAB layer definitions');
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
        'ValidationFrequency',30,...
        'MiniBatchSize',100,...
        'Verbose',false,...
        'Plots','training-progress',...
        'ValidationData',imdsValid,...
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
