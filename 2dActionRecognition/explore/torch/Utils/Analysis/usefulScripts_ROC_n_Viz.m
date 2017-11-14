
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
w1 = net.Layers(6).Weights;
w1=mat2gray(w1);
rw1=imresize(w1,2);
figure; montage(rw1(:,:,2,:));
figure; montage(rw1(:,:,3,:));
figure; montage(rw1(:,:,4,:));


% layer 2
w1 = net.Layers(5).Weights;
w1=mat2gray(w1);
montage(imresize(w1(:,:,33,:), 5));

%% ROC Curve
plotroc((valData.Labels == 'lowMet')',score(:,2)');
[X,Y,T,AUC] = perfcurve((valData.Labels == 'lowMet'),score(:,2),1);