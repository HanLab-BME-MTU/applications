%% Looks at results

function [h_viz] = vizActivitationLater(net, layerNum, imds, imageIdx, varargin) 

numMontage = 8;
numCol = 20;

act1 = activations(net,imds.readimage(imageIdx),layerNum,'OutputAs','channels');  
sz = size(act1);
act1 = reshape(act1,[sz(1) sz(2) 1 sz(3)]);

h_viz = figure; montage(mat2gray(act1),'Size',[round(sz(3)/numCol)+1 numCol]);
figure; imshow(imds.readimage(imageIdx),[]);