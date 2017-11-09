% simple script to pre-process the cells for Deep learning.

% Load image
% imNum = 2695;
% MD = load(cellDataSet{imNum}.cellMD);
% I = gpuArray(mat2gray(MD.MD.getChannel(1).loadImage(1)));
% imshow(imresize(I, [64 64]));

%% Gen Data
load('/work/bioinformatics/shared/dope/data/OMETIFF/Gen2n3_May15_ALL.mat', 'cellDataSet');

resizeOn = true;
newDim = [128 128];

% dataRootDir = '/work/bioinformatics/shared/dope/torch/test/217x217/single/';
dataRootDir = '/work/bioinformatics/shared/dope/torch/test/128x128/';
dataRootDirVal = '/work/bioinformatics/shared/dope/torch/test/128x128/validationSet';
dataRootDirTR = '/work/bioinformatics/shared/dope/torch/test/128x128/trainingSet';

randOrd = randperm(length(cellDataSet));
percentVal = .25;

parfor i = 1:length(cellDataSet)
    MD = load(cellDataSet{i}.cellMD,'MD');
    %     I = gpuArray(mat2gray(MD.MD.getChannel(1).loadImage(1)));
    MD = MD.MD;
    
    if rand > percentVal
        dataSetDir = 'trainingSet';
    else
        dataSetDir = 'validationSet';
    end

    data_lowMet = fullfile(dataRootDir,dataSetDir,'lowMet');
    data_highMet = fullfile(dataRootDir,dataSetDir,'highMet');
    data_unMet = fullfile(dataRootDir,dataSetDir,'unMet');        
    
    for fidx = 1:MD.nFrames_
        I = mat2gray(MD.getChannel(1).loadImage(fidx));
        if resizeOn
            I = imresize(I, newDim);
        end
        
        frameNum = num2str(fidx);
        newFileOut = [MD.processes_{1}.funParams_.key '_f' frameNum '.png'];

        if MD.processes_{1}.funParams_.metEff == 0
            newFileOut = fullfile(data_lowMet, newFileOut);
        elseif MD.processes_{1}.funParams_.metEff == 1
            newFileOut = fullfile(data_highMet, newFileOut);
        else
            newFileOut = fullfile(data_unMet, newFileOut);
        end
        imwrite(I, newFileOut);
    end
    
end


%% Gather image data
imdsTrain = imageDatastore(dataRootDirTR, 'IncludeSubfolders',true,'FileExtensions','.png','LabelSource','foldernames');
imdsValid = imageDatastore(dataRootDirVal, 'IncludeSubfolders',true,'FileExtensions','.png','LabelSource','foldernames');
% get rid of unMet...
% [imds1] = splitEachLabel(imgs,1500,'Exclude','unMet');
% [trainData, valData] = splitEachLabel(imgs, 1500, 'randomize');



%% define network
layers = [ ...
    imageInputLayer([128 128 1])
    convolution2dLayer(7, 6, 'Stride', 1, 'Name','conv1')
%     batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(3,'Stride',2,'Name','maxPool1')
    convolution2dLayer(5, 16,'Stride',2,'Name','conv2')
%     batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(3,'Stride',2,'Name','maxPool2')    
    convolution2dLayer(3,32,'Name','conv3')
    reluLayer
    maxPooling2dLayer(3,'Stride',2,'Name','maxPool2')
    reluLayer
    dropoutLayer
%     convolution2dLayer(3,25,'Name','conv4')
%     reluLayer
%     convolution2dLayer(3,25,'Name','conv5','Padding','same')
%     reluLayer
%     maxPooling2dLayer(3,'Stride', 2,'Name','maxPool3','Padding','same')
    fullyConnectedLayer(256,'Name','fc1')
    reluLayer
    dropoutLayer
    fullyConnectedLayer(256,'Name','fc2')
    reluLayer
    dropoutLayer
    fullyConnectedLayer(3,'Name','fc3')
    softmaxLayer
    classificationLayer];

% 2017b
% options = trainingOptions('sgdm',...
%     'MaxEpochs',1000, ...
%     'ValidationFrequency',30,...
%     'MiniBatchSize',500,...
%     'Verbose',false,...
%     'Plots','training-progress',...
%     'ValidationData',imdsValid,'CheckpointPath',fullfile(dataRootDir,'TEMP_val'));

%% Training options
%2017a
options = trainingOptions('sgdm','LearnRateSchedule','piecewise',...
      'LearnRateDropFactor',0.2,'LearnRateDropPeriod',5,... 
      'MaxEpochs',1000,'MiniBatchSize',100,'CheckpointPath',fullfile(dataRootDir,'TEMP_val'));    

trainedNet = trainNetwork(imdsTrain,layers,options);



%% Looks at results
  
dI = deepDreamImage(trainedNet,10,[1 2 3],'Verbose',false);  
        

imageIdx = 20030;
layerNum = 8;
numMontage = 8;
numCol = 20;
act1 = activations(trainedNet,imgs.readimage(imageIdx),layerNum,'OutputAs','channels');  
sz = size(act1);
act1 = reshape(act1,[sz(1) sz(2) 1 sz(3)]);
figure(5); montage(mat2gray(act1),'Size',[round(sz(3)/numCol)+1 numCol]);
figure(4);imshow(imgs.readimage(imageIdx),[]);
  

%% Evaluation & Visualization

% layer 2
w1 = net.Layers(2).Weights;
w1=mat2gray(w1);
rw1=imresize(w1,2);
figure; montage(rw1);


% layer 2
w1 = net.Layers(5).Weights;
w1=mat2gray(w1);
montage(imresize(w1(:,:,33,:), 5));

%% ROC Curve
plotroc((valData.Labels == 'lowMet')',score(:,2)');
[X,Y,T,AUC] = perfcurve((valData.Labels == 'lowMet'),score(:,2),1);

%% Training Autoencoders

I = trainData.readall;
auto1=trainAutoencoder(I, 200);
x_recon = predict(auto1,valData.readimage(1));
% Iv = valData.readall;
valData.ReadFcn = @customReadDouble
x_recon = predict(auto1,valData.readimage(10));
imshow(x_recon,[])

