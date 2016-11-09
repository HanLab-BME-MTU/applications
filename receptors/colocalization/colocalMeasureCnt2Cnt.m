function [r, randR] = colocalMeasureCnt2Cnt(imageCntA,imageCntB,maskingFile)
% COLOCALMEASURECNT2CNT measures colocalization for two channels where both are continuous
%

    % Read in continuum image
    ICntA = double(imageCntA); 
    ICntB = double(imageCntB); 
    
    %Correct images for non-uniform background
    compValue = mean(ICntA(maskingFile~=0));
    compMask = compValue*ones(size(ICntA,1),size(ICntA,2));
    nImage = filterGauss2D(ICntA,10);
    ICntA = ICntA-nImage;
    ICntA = ICntA+ compMask;
    
    compValue = mean(ICntB(maskingFile~=0));
    compMask = compValue*ones(size(ICntB,1),size(ICntB,2));
    nImage = filterGauss2D(ICntB,10);
    ICntB = ICntB-nImage;
    ICntB = ICntB+ compMask;
    
    
    r = corr2(ICntA,ICntB);
    [~,ImageRand] = randomizeImageBlocks(ICntA,maskingFile,[16,16]);
    randR = corr2(ImageRand,ICntB);
    
function [NewCellMask,imageRand]=randomizeImageBlocks(Image,CellMask,blockSize)
%% SYNOPSIS
%INPUT
%       Image; image to be randomized
%    CellMask; binary map of ROI 
%   blockSize; desired block size to be randomized, e.g if 4x4, input [4,4]
%   
%OUTPUT
%     NewCellMask; new cell mask binary fittable with desired blocks size within the input
%     CellMask
%     ImageRand;Randomized blocks of input Image within output NewCellMask
%Note: No matter how may times you loop this function, if blockSize is
%same, NewCellMask will always be the same for any given Image.
% John Maringa Githaka, January 2015;
%%
Image=double(Image);
CellMask=double(logical(CellMask)); %Convert to binary just incase input wasn't in binary.
[m,n]=size(Image);
%% divide CellMask into blocks
a=ones(1,m/blockSize(1,1))*blockSize(1,1);
b=ones(1,n/blockSize(1,2))*blockSize(1,2);

cellMaskBlocks = mat2cell(CellMask,a,b);

A=cellMaskBlocks(:);
BlockArea=(blockSize(1,1)*blockSize(1,2));
Rand=1:length(A);
%Find index of cell 'blocks' within CellMask and get the NewCellMask 
for i=1:length(A);
    if sum(A{i,1}(:))<BlockArea; 
        A{i,1}(:)=0; 
        Rand(Rand==i)=[];
    end
end
NewCellMask=reshape(A,size(cellMaskBlocks));
NewCellMask=cell2mat(NewCellMask);
% Randomize Image Blocks within NewCellMask within NewCellMask
imageRand={};
idxRand=Rand(randperm(length(Rand))); %randomize the indices
imageBlocks = mat2cell(Image,a,b);
b=imageBlocks(:);
reshuffle=1;
for i=1:length(A);
    if sum(A{i,1}(:))==BlockArea; 
    imageRand{i,1}=b{idxRand(reshuffle),1};
    reshuffle=reshuffle+1;
    else
    imageRand{i,1}=zeros(blockSize(1,1),blockSize(1,2));
    end
end
imageRand=reshape(imageRand,size(cellMaskBlocks));
imageRand=cell2mat(imageRand);
    

