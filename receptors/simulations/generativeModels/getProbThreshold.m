function [thresholds, sensitivity, specificity, probMerged] = getProbThreshold(numFrames)

%% Parameters
locError = 0.02; %measurement error
diffConst = .1; %diffusion constant
dT = 0.01; %time between frames
maxLagFrames = 3;
maxWinSize = 20;
minWinSize = 14;
boxSize = 10;

%% Simulation
%get red and green trajectories and add error; G1 is bound to R, G2 not
[paths,~] = brownianMotionBox(2,diffConst,round(numFrames*dT),dT,boxSize,2,[.5 .5]*boxSize);
pathR = squeeze(paths(:,:,1)) + randn(size(squeeze(paths(:,:,1))))*locError;
pathG1 = squeeze(paths(:,:,1)) + randn(size(pathR))*locError;
pathG2 = squeeze(paths(:,:,2)) + randn(size(pathR))*locError;

separations1 = sqrt(sum((pathR-pathG1).^2,2));
separations2 = sqrt(sum((pathR-pathG2).^2,2));

numFrames = length(pathR);

%% Get Best Threshold
%initialize threshold and statistic matrices
thresholds = zeros(maxLagFrames,maxWinSize-minWinSize+1);
sensitivity = zeros(size(thresholds));
specificity = zeros(size(thresholds));

%loop through available window sizes and lag times for vectors
for winSize = minWinSize:maxWinSize
    winIndx = winSize-minWinSize+1; % done so that a minimum window size can be set
    for lagFrames = 1:maxLagFrames
        fprintf('%g : %g\n',winSize,lagFrames); %to see progress
        
        %calculates the correlation statistic across the trajectories; zero
        %if outside the range (N is much bigger than this, so only neglible
        %effect
        corr1 = zeros(numFrames,1);
        corr2 = zeros(numFrames,1);
        
        %corr1(winSize:(numFrames-winSize-lagFrames)) = trajCorrelation(pathR,pathG1,winSize:(numFrames-winSize-lagFrames),winSize,lagFrames);
        %corr2(winSize:(numFrames-winSize-lagFrames)) = trajCorrelation(pathR,pathG2,winSize:(numFrames-winSize-lagFrames),winSize,lagFrames);
        
        corr1(1:(numFrames-lagFrames)) = probMergedByDist(boxSize,diffConst,dT,separations1((lagFrames+1):numFrames),separations1(1:(numFrames-lagFrames)));
        corr2(1:(numFrames-lagFrames)) = probMergedByDist(boxSize,diffConst,dT,separations2((lagFrames+1):numFrames),separations2(1:(numFrames-lagFrames)));
        
        %determines which threshold value maximizes the sum of sensitivity
        %and specificity
        minError = 5;
        for thresh = max(corr2):-0.001:min(corr1)
            error = mean(corr1<thresh) + mean(corr2>thresh);
            if minError > error
                thresholds(lagFrames,winIndx) = thresh;
                minError = error;
                sensitivity(lagFrames,winIndx) = mean(corr1>thresh);
                specificity(lagFrames,winIndx) = mean(corr2<thresh);
            end
        end
        
    end
end

%pick the best set of window size and frame lag to maximize sensitivity and
%specificity
sumMatrix = sensitivity+specificity;
[~,indx] = max(sumMatrix(:));
[bestLagFrames, bestWinSizeIndx] = ind2sub(size(sensitivity),indx);
bestWinSize = bestWinSizeIndx + minWinSize-1;

%since the correlation values aren't stored, this recalculates it
corr1 = zeros(numFrames,1);
corr2 = zeros(numFrames,1);
corr1(1:(numFrames-bestLagFrames)) = probMergedByDist(boxSize,diffConst,dT,separations1((bestLagFrames+1):numFrames),separations1(1:(numFrames-bestLagFrames)));
corr2(1:(numFrames-bestLagFrames)) = probMergedByDist(boxSize,diffConst,dT,separations2((bestLagFrames+1):numFrames),separations2(1:(numFrames-bestLagFrames)));

%get the sensitivity, specificity, and threshold from the best set
bestSensitivity = sensitivity(bestLagFrames,bestWinSizeIndx);
bestSpecificity = specificity(bestLagFrames,bestWinSizeIndx);
bestThreshold = thresholds(bestLagFrames,bestWinSizeIndx);

%Create vector of probabilities that a given correlation coefficient is
%merged vs. unmerged. Resolution of 0.005
probMerged = zeros(200,2);
probMerged(:,1) = (1:1:200)/200;

%do the first case because 0 isn't included in probMerged
numMerged = sum(sum(and(corr1>0,corr1<probMerged(1,1))));
numIndep = sum(sum(and(corr2>0,corr2<probMerged(1,1))));
if numIndep+numMerged == 0
    probMerged(1,2) = 0;
else
    probMerged(1,2) = numMerged/(numIndep+numMerged);
end

for i = 2:200
    numMerged = sum(sum(and(corr1>probMerged(i-1,1),corr1<probMerged(i,1))));
    numIndep = sum(sum(and(corr2>probMerged(i-1,1),corr2<probMerged(i,1))));
    if numIndep+numMerged == 0
        probMerged(i,2) = 0;
    else
        probMerged(i,2) = numMerged/(numIndep+numMerged);
    end
end

%plot the correlation coefficient over time for the best set of parameters
%and place in the title the important statistics; it also draws the
%threshold determined;
figure
plot(1:numFrames,corr1','.',1:numFrames,corr2','.',[1 numFrames],[bestThreshold bestThreshold],'-')
title(sprintf('LagFrames = %g, WinSize = %g, Sensitivity = %g, Specificity = %g',bestLagFrames,bestWinSize,bestSensitivity,bestSpecificity));

figure
plot(1:numFrames,separations1,1:numFrames,separations2);
xlim([0 numFrames]);
ylim([0 boxSize*sqrt(2)]);

%%%%%% the end %%%%%%