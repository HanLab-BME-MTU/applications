
% simple script to pre-process the cells for Deep learning.
% Andrew R. Jamieson Oct 2017


clear;
clc;


%% load data state 
[filename, pathname] = uigetfile('/work/bioinformatics/shared/dope/torch/test/217x217/segmented/','Which data to load?');
load(fullfile(pathname,filename));

%  number of classes
numClass = 2;

imdsTrain_Tumor_v_Mel = imdsTrain;
imdsVal_Tumor_v_Mel = imdsValid;

tumorLabel = imdsTrain_Tumor_v_Mel.Labels == 'highMet' | ...
    imdsTrain_Tumor_v_Mel.Labels == 'lowMet';
imdsTrain_Tumor_v_Mel.Labels(tumorLabel) = 'tumor';

tumorLabel = imdsVal_Tumor_v_Mel.Labels == 'highMet' |...
    imdsVal_Tumor_v_Mel.Labels == 'lowMet';
imdsVal_Tumor_v_Mel.Labels(tumorLabel) = 'tumor';

% [imdsTrain_Tumor_v_Mel] = splitEachLabel(imdsTrain_Tumor_v_Mel,.9999999,...
%                             'Exclude',{'highMet','lowMet'});
% [imdsVal_Tumor_v_Mel] = splitEachLabel(imdsVal_Tumor_v_Mel,.9999999,...
%                             'Exclude',{'highMet','lowMet'});
uisave();

%% Define Checkpoint
[checkPointDir] = uigetdir('/work/bioinformatics/shared/dope/torch/test/217x217/segmented/netTEMP',...
                            'Where to store net Checkpoints?');
                                                                      
                        
if ~verLessThan('matlab', '9.3')
    % define network
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
        fullyConnectedLayer(2,'Name','fc3')
        softmaxLayer
        classificationLayer]
end

%% Training options
if ~verLessThan('matlab', '9.3')
                        % 2017b
    options = trainingOptions('sgdm',...
        'MaxEpochs',100000, ...
        'ValidationFrequency',2500,...
        'MiniBatchSize',50,...
        'Verbose',true,...
        'Plots','training-progress',...
        'ValidationData',imdsVal_Tumor_v_Mel,...
        'ValidationPatience', Inf,...
        'ExecutionEnvironment' , 'multi-gpu',...
        'CheckpointPath', checkPointDir);
else 
%     2017a
    options = trainingOptions('sgdm','LearnRateSchedule','piecewise',...
          'LearnRateDropFactor',0.2,'LearnRateDropPeriod',5,... 
          'MaxEpochs',1000,'MiniBatchSize',75,...
          'CheckpointPath',checkPointDir);    
end

%% Start Training
trainedNet = trainNetwork(imdsTrain_Tumor_v_Mel,layers,options);
