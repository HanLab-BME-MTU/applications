function [matPosMIP,matNegMIP] = getMIP(matImStack,MOTION_THRESHOLD)
%takes an image stack (binary conversion of uint16) and calculates 16 MIP
%channels

%initialize the MIP matricies (dim: M x N x 8 directions x 8 alphas)
matPosMIP = zeros(length(matImStack(:,1,1)),length(matImStack(1,:,1)),8,8);
matNegMIP = zeros(length(matImStack(:,1,1)),length(matImStack(1,:,1)),8,8);

%loops through the entire video skipping 1st and last frames
for idFrame = 2:length(matImStack(1,1,1:end-1))
    
    %find the current, previous, and next frames
    matPrevF    = matImStack(:,:,idFrame - 1);
    matCurrentF = matImStack(:,:,idFrame);
    matNextF    = matImStack(:,:,idFrame + 1);
    
    %loop through 8 alphas
    %0:7 because next function uses circrotate and we want no rotation for
    %the first channel
    for alpha = 0:7
        temp_matMIP = calculateMIP(matPrevF, matCurrentF, matNextF, alpha, MOTION_THRESHOLD);
        [temp_matPos,temp_matNeg] = mapValMIP(temp_matMIP);
        
        %sum across all time points in given stack
        matPosMIP(:,:,:,alpha+1) = matPosMIP(:,:,:,alpha+1) + temp_matPos;
        matNegMIP(:,:,:,alpha+1) = matNegMIP(:,:,:,alpha+1) + temp_matNeg;
    end
end

%sum across the directions to form 8 rotation invariant channels
matPosMIP = squeeze(sum(matAllMIP,3));
matNegMIP = squeeze(sum(matNegMIP,3));
end


function [matPos,matNeg] = mapValMIP(matMIP)
%splits the aggregate MIP channel into positive and negative MIP channels

matPos = (matMIP > 0);
matNeg = (matMIP < 0);
end


function matMIP = calculateMIP(matPrevF, matCurrentF, matNextF, ALPHA ,MOTION_THRESHOLD)
%Computes a MIP image given frames to compare, a threshold, and a
%comparison angle ALPHA

%Given MOTION_THRESHOLD is intensity between 2 pixels
%Square and multiply by 9 since we compute SSD between 3x3 patches
MOTION_THRESHOLD = MOTION_THRESHOLD*MOTION_THRESHOLD*9;

%Define the ring of 6 positions to compare for the previous and next frames
vecPatchPositions = [0 4; 3 3; 4 0; 3 -3; 0 -4; -3 -3; -4 0; -3 3];
vecPrevFPositions = vecPatchPositions;
vecNextFPositions = circshift(vecPatchPositions,ALPHA);

%Initialize the value storage matricies (dim = M x N x 8 directions)
matMIP = zeros(size(matPrevF,1),size(matPrevF,2),8);
matReverseMIP = zeros(size(matMIP));

for patch = 1:8
    matMIP(:,:,patch) = calculateOnePatch(matPrevF, matCurrentF, matNextF, vecPrevFPositions(patch,:), vecNextFPositions(patch,:), MOTION_THRESHOLD);
end

% Switch patches position for static edge suppression
vecPrevFPositions = vecNextFPositions;
vecNextFPositions = vecPatchPositions;
for patch = 1:8
    matReverseMIP(:,:,patch) = calculateOnePatch(matPrevF, matCurrentF, matNextF, vecPrevFPositions(patch,:), vecNextFPositions(patch,:), MOTION_THRESHOLD);
end

% Suppress votes where different signs and both not 0
matMIP(matMIP ~= 0 & matReverseMIP ~= 0  & matMIP ~=  matReverseMIP ) = 0;
end


function matOnePatch = calculateOnePatch(matPrevF, matCurrentF, matNextF, pointPrevFPosition, pointNextFPosition, MOTION_THRESHOLD)
%computes a MIP image for 1 direction (8 directions in 1 alpha channel)

%Calculates the sum of squared difference (SSD) between the previous
%frame & current frame and the next frame & current frame

%Rotates the previous and next frames by their rotation to get patch
%comparisons
SSD1 = ( circshift(matPrevF,pointPrevFPosition) - matCurrentF).^2;
SSD2 = ( circshift(matNextF,pointNextFPosition) - matCurrentF).^2;

%Sums the 3x3 grid around each pixel
%Edges wrap around
PA_SSD1 = 0;
PA_SSD2 = 0;
for i = -1:1
    for j = -1:1
        PA_SSD1 = PA_SSD1 + circshift(SSD1, [i j]);
        PA_SSD2 = PA_SSD2 + circshift(SSD2, [i j]);
    end
end

%Finds the SSD difference
PA_SSD = PA_SSD1 - PA_SSD2;

%Calculate where ||SSD1 - SSD2|| > treshold
matPosMot = double(PA_SSD > MOTION_THRESHOLD);
matNegMot = double(PA_SSD < -MOTION_THRESHOLD);

%Groups positive and negative channels
matOnePatch = matPosMot - matNegMot;
end
