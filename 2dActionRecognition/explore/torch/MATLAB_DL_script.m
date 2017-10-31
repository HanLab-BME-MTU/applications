% simple script to pre-process the cells for Deep learning.

% Load image
% imNum = 2695;
% MD = load(cellDataSet{imNum}.cellMD);
% I = gpuArray(mat2gray(MD.MD.getChannel(1).loadImage(1)));
% imshow(imresize(I, [64 64]));

%% Gen Data
load('/work/bioinformatics/shared/dope/data/OMETIFF/Gen2n3_May15_ALL.mat', 'cellDataSet');

resizeOn = false;
newDim = [128 128];

dataRootDir = '/work/bioinformatics/shared/dope/torch/test/217x217/single/';
% dataRootDir = '/work/bioinformatics/shared/dope/torch/test/128x128/';
data_lowMet = fullfile(dataRootDir, 'lowMet');
data_highMet = fullfile(dataRootDir, 'highMet');
data_unMet = fullfile(dataRootDir, 'unMet');

parfor i = 1:length(cellDataSet)
    MD = load(cellDataSet{i}.cellMD,'MD');
    %     I = gpuArray(mat2gray(MD.MD.getChannel(1).loadImage(1)));
    MD = MD.MD;
    
    for fidx = 1:1%MD.nFrames_
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

imgs =imageDatastore(dataRootDir, 'IncludeSubfolders',true,'FileExtensions','.png','LabelSource','foldernames');
% get rid of unMet...
% [imds1] = splitEachLabel(imgs,1500,'Exclude','unMet');
[trainData, valData] = splitEachLabel(imgs, 1500, 'randomize');



%% define network
layers = [ ...
    imageInputLayer([217 217 1])
    convolution2dLayer(7,96,'Name','conv1')
    reluLayer
    maxPooling2dLayer(3,'Stride',2,'Name','maxPool1')
    convolution2dLayer(5,128,'Name','conv2')
    reluLayer
    maxPooling2dLayer(3,'Stride',2,'Name','maxPool2')
    convolution2dLayer(3,256,'Name','conv3')
    reluLayer
    convolution2dLayer(3,256,'Name','conv4')
    reluLayer
    convolution2dLayer(3,128,'Name','conv5')
    reluLayer
    maxPooling2dLayer(3,'Stride',3,'Name','maxPool3')
    fullyConnectedLayer(512,'Name','fc1')
    reluLayer
    dropoutLayer
    fullyConnectedLayer(512,'Name','fc2')
    reluLayer
    dropoutLayer
    fullyConnectedLayer(3,'Name','fc3')
    softmaxLayer
    classificationLayer];

%   20x1 Layer array with layers:
% 
%      1   ''   Image Input             217x217x1 images with 'zerocenter' normalization
%      2   ''   Convolution             50 5x5 convolutions with stride [1  1] and padding [0  0]
%      3   ''   ReLU                    ReLU
%      4   ''   Max Pooling             3x3 max pooling with stride [2  2] and padding [0  0]
%      5   ''   Convolution             100 5x5 convolutions with stride [1  1] and padding [0  0]
%      6   ''   ReLU                    ReLU
%      7   ''   Max Pooling             3x3 max pooling with stride [2  2] and padding [0  0]
%      8   ''   Convolution             96 3x3 convolutions with stride [1  1] and padding [0  0]
%      9   ''   ReLU                    ReLU
%     10   ''   Convolution             96 3x3 convolutions with stride [1  1] and padding [0  0]
%     11   ''   ReLU                    ReLU
%     12   ''   Convolution             128 3x3 convolutions with stride [1  1] and padding [0  0]
%     13   ''   ReLU                    ReLU
%     14   ''   Max Pooling             3x3 max pooling with stride [3  3] and padding [0  0]
%     15   ''   Fully Connected         384 fully connected layer
%     16   ''   ReLU                    ReLU
%     17   ''   Dropout                 50% dropout
%     18   ''   Fully Connected         3 fully connected layer
%     19   ''   Softmax                 softmax
%     20   ''   Classification Output   crossentropyex

% 2017b
% options = trainingOptions('sgdm',...
%     'MaxEpochs',3, ...
%     'ValidationFrequency',30,...
%     'Verbose',false,...
%     'Plots','training-progress');
% %     'ValidationData',valData,...


%% Training options
%2017a
options = trainingOptions('sgdm',...
      'LearnRateSchedule','piecewise',...
      'LearnRateDropFactor',0.2,... 
      'LearnRateDropPeriod',5,... 
      'MaxEpochs',1000,... 
      'MiniBatchSize',100,...
      'CheckpointPath',fullfile(dataRootDir,'TEMP'));    

trainedNet = trainNetwork(trainData,layers,options);

s

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

