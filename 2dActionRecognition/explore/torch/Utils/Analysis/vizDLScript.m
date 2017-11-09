%% Looks at results

load('/work/bioinformatics/shared/dope/torch/test/217x217/TEMP/convnet_checkpoint__21600__2017_10_26__08_57_27.mat')
load('/work/bioinformatics/shared/dope/torch/test/217x217/preTrainingData_217x217_Oct25th2017.mat')



dI = deepDreamImage(trainedNet,10,[1 2 3],'Verbose',false);  
dI = deepDreamImage(trainedNet,22,[1:3],'Verbose',false);montage(dI);

imageIdx = 20030;
layerNum = 8;
numMontage = 8;
numCol = 20;
act1 = activations(trainedNet,imgs.readimage(imageIdx),layerNum,'OutputAs','channels');  
sz = size(act1);
act1 = reshape(act1,[sz(1) sz(2) 1 sz(3)]);
figure(5); montage(mat2gray(act1),'Size',[round(sz(3)/numCol)+1 numCol]);
figure(4);imshow(imgs.readimage(imageIdx),[]);