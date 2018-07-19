% Andrew R. Jamieson Nov. 2017

clear;
clc;

%% Define experiment dir
[exprDir] = uigetdir('/work/bioinformatics/shared/dope/torch/test/217x217/segmented/',...
                            'Where to store experimental output and confis?');


%% load data state 
% [filename, pathname] = uigetfile('/work/bioinformatics/shared/dope/torch/test/217x217/segmented/','Which data to load?');
% load(fullfile(pathname,filename));
load('/work/bioinformatics/shared/dope/torch/test/217x217/segmented/Experiments/High_Low_Melo/preTrainingDataStructure_217x217_segmented01-Nov-2017-0908.mat');

%  number of classes
numClass = 3;


%% Define Checkpoint
[checkPointDir] = uigetdir(exprDir,...
                            'Where to store net Checkpoints?');
                                                                      


if ~verLessThan('matlab', '9.3')

    %% define network
    disp('Note: using 2017b MATLAB layer definitions');
    layers = [ ...netTempDir
        imageInputLayer([217 217 1])
        convolution2dLayer(7, 24, 'Stride', 2,'Name','conv1')
        batchNormalizationLayer
        reluLayer
        maxPooling2dLayer(3,'Stride',2,'Name','maxPool1')
        convolution2dLayer(5, 48,'Stride',2,'Name','conv2')
        batchNormalizationLayer
        reluLayer
        maxPooling2dLayer(3,'Stride',2,'Name','maxPool11')
        convolution2dLayer(3,48,'Name','conv3')
        batchNormalizationLayer
        reluLayer
        maxPooling2dLayer(3,'Stride',2,'Name','maxPool3')
        convolution2dLayer(3,96,'Name','conv4','Padding','same')
        reluLayer
        convolution2dLayer(3,96,'Name','conv5','Padding','same')
        reluLayer
        maxPooling2dLayer(3,'Stride', 2,'Name','maxPool4')
        fullyConnectedLayer(1024,'Name','fc1')
        reluLayer
        dropoutLayer
        fullyConnectedLayer(1024,'Name','fc2')
        reluLayer
        dropoutLayer
        fullyConnectedLayer(numClass,'Name','fc3')
        softmaxLayer
        classificationLayer]

end

%% Training options
if ~verLessThan('matlab', '9.3')
    options = trainingOptions('sgdm',...
        'MaxEpochs',1000, ...
        'ValidationFrequency',5000,...
        'MiniBatchSize',50,...
        'Verbose',false,...
        'Plots','training-progress',...
        'ValidationData',imdsValid,...
        'ValidationPatience', Inf,...
        'ExecutionEnvironment' , 'multi-gpu',...
        'CheckpointPath', checkPointDir);
end

saveConfigFile = fullfile(exprDir,['preTrainConfig_' datestr(datetime,'dd-mmm-yyyy-hhMM') '.mat']);
disp(['Saving preRun configuration state....' saveConfigFile]);
save();

%% Start Training
trainedNet = trainNetwork(imdsTrain_Tumor,layers,options);
