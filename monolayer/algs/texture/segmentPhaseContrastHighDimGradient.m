%% Input:
% I - input image
% patchSize - size of patches
% 
% Output: gradient image based on high dimensional gray-level histograms

function [GRAD] = segmentPhaseContrastHighDimGradient(I,patchSize)

if size(I,3) == 3
    I = rgb2gray(I);
end

p1 = prctile(I(:),1);
p99 = prctile(I(:),99);
I(I < p1) = p1;
I(I > p99) = p99;
centers = p1 : (p99-p1)/18 : p99;

% get histogram for each patch in the image
ncenters = length(centers);
for l = 1 : ncenters
    if l == 1
        minVal = -inf;
        maxVal = centers(l) + (centers(l+1) - centers(l))/2;
    else
        if l == ncenters
            minVal = centers(l-1) + (centers(l) - centers(l-1))/2;
            maxVal = inf;            
        else
            minVal = centers(l-1) + (centers(l) - centers(l-1))/2;
            maxVal = centers(l) + (centers(l+1) - centers(l))/2;            
        end
    end    
    fun = @(block_struct) ...
        sum(sum(block_struct.data >= minVal & block_struct.data < maxVal))/(patchSize^2);
    tmp = blockproc(I,[patchSize patchSize],fun);
    if ~exist('Ihist','var')
        Ihist = nan(size(tmp,1),size(tmp,2),ncenters);
    end
    Ihist(:,:,l) = tmp;
    %     figure; title(sprintf('bin #%d',l)); imagesc(tmp);
    %     eval(sprintf('print -dbmp16m %s', ['~/tmp/segmentation/bin' int2str(l) '.bmp']));
end

mask = [[1 2 1];[0 0 0];[-1 -2 -1]];
Deriv = highDimensionDerivative(Ihist,mask);

GRAD = imresize(Deriv,size(I));
end

function [DERIV] = highDimensionDerivative(I,mask)
[sizeY,sizeX,nDim] = size(I);
allDeriv  = nan(sizeY,sizeX,nDim);
for d = 1 : nDim
    Id = I(:,:,d);
    energy1 = abs(imfilter(Id,mask,'same'));
    energy2 = abs(imfilter(Id,mask','same'));
    allDeriv(:,:,d) = sqrt(energy1.^2+energy2.^2);
end
DERIV = sum(allDeriv,3);
end