function [thresholds, sensitivity, specificity] = getCorrCoeffThreshold()

%% Parameters
locError = 0.02; %measurement error
diffConst = .1; %diffusion constant
dT = 0.01; %time between frames
maxLagFrames = 3;
maxWinSize = 20;
minWinSize = 14;
boxSize = 10;
numFrames = 500000; %time steps

%% Simulation
%get red and green trajectories and add error; G1 is bound to R, G2 not
[paths,~] = brownianMotionBox(2,diffConst,round(numFrames*dT),dT,boxSize,2,[.5 .5]*boxSize);
pathR = squeeze(paths(:,:,1)) + randn(size(squeeze(paths(:,:,1))))*locError;
pathG1 = squeeze(paths(:,:,1)) + randn(size(pathR))*locError;
pathG2 = squeeze(paths(:,:,2)) + randn(size(pathR))*locError;

[numFrames, ~] = size(pathR);

separations1 = sqrt(sum((pathR-pathG1).^2,2));
separations2 = sqrt(sum((pathR-pathG2).^2,2));

%% Get Best Threshold
%initialize threshold and statistic matrices
thresholds = zeros(maxLagFrames,maxWinSize-minWinSize+1,3);
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
        corrDist1 = zeros(numFrames,1);
        corrDist2 = zeros(numFrames,1);
                
        corr1(winSize:(numFrames-winSize-lagFrames)) = trajCorrelation(pathR,pathG1,winSize:(numFrames-winSize-lagFrames),winSize,lagFrames);
        corr2(winSize:(numFrames-winSize-lagFrames)) = trajCorrelation(pathR,pathG2,winSize:(numFrames-winSize-lagFrames),winSize,lagFrames);
       
        corr1CutOff = corr1.*(separations1<6*locError);
        corr2CutOff = corr2.*(separations2<6*locError);
        
        corrDist1(1:(numFrames-lagFrames)) = probMergedByDist(boxSize,diffConst,dT,separations1((lagFrames+1):numFrames),separations1(1:(numFrames-lagFrames)));
        corrDist2(1:(numFrames-lagFrames)) = probMergedByDist(boxSize,diffConst,dT,separations2((lagFrames+1):numFrames),separations2(1:(numFrames-lagFrames)));
                    
        %determines which threshold value maximizes the sum of sensitivity
        %and specificity
        minError = realmax('double');
        for thresh = max(corr2):-0.001:min(corr1)
            error = mean(corr1<thresh)*sum(corr2~=0) + mean(corr2>thresh)*sum(corr1~=0);
            if minError > error
                thresholds(lagFrames,winIndx,1) = thresh;
                minError = error;
                sensitivity(lagFrames,winIndx,1) = mean(corr1>thresh);
                specificity(lagFrames,winIndx,1) = mean(corr2<thresh);
            end
        end
        
        minError = realmax('double');
        for thresh = max(corr2CutOff):-0.001:min(corr1CutOff)
            error = mean(corr1CutOff<thresh)*sum(corr2CutOff~=0) + mean(corr2CutOff>thresh)*sum(corr1CutOff~=0);
            if minError > error
                thresholds(lagFrames,winIndx,2) = thresh;
                minError = error;
                sensitivity(lagFrames,winIndx,2) = mean(corr1CutOff>thresh);
                specificity(lagFrames,winIndx,2) = mean(corr2CutOff<thresh);
            end
        end
        
        minError = realmax('double');
        for thresh = max(corrDist2):-0.001:min(corrDist1)
            error = mean(corrDist1<thresh) + mean(corrDist2>thresh);
            if minError > error
                thresholds(lagFrames,winIndx,3) = thresh;
                minError = error;
                sensitivity(lagFrames,winIndx,3) = mean(corrDist1>thresh);
                specificity(lagFrames,winIndx,3) = mean(corrDist2<thresh);
            end
        end
        
    end
end

%pick the best set of window size and frame lag to maximize sensitivity and
%specificity
sumMatrix = squeeze(sensitivity(:,:,1)+specificity(:,:,1));
[~,indx] = max(sumMatrix(:));
[bestLagFrames1, bestWinSizeIndx1, ~] = ind2sub(size(sensitivity(:,:,1)),indx);
bestWinSize1 = bestWinSizeIndx1 + minWinSize-1;

sumMatrix = squeeze(sensitivity(:,:,2)+specificity(:,:,2));
[~,indx] = max(sumMatrix(:));
[bestLagFrames2, bestWinSizeIndx2, ~] = ind2sub(size(sensitivity(:,:,2)),indx);
bestWinSize2 = bestWinSizeIndx2 + minWinSize-1;

sumMatrix = squeeze(sensitivity(:,:,3)+specificity(:,:,3));
[~,indx] = max(sumMatrix(:));
[bestLagFrames3, bestWinSizeIndx3, ~] = ind2sub(size(sensitivity(:,:,3)),indx);
bestWinSize3 = bestWinSizeIndx3 + minWinSize-1;

%since the correlation values aren't stored, this recalculates it
corr1 = zeros(numFrames,1);
corr2 = zeros(numFrames,1);
corrDist1 = zeros(numFrames,1);
corrDist2 = zeros(numFrames,1);

corr1(bestWinSize1:(numFrames-bestWinSize1-bestLagFrames1)) = trajCorrelation(pathR,pathG1,bestWinSize1:(numFrames-bestWinSize1-bestLagFrames1),bestWinSize1,bestLagFrames1);
corr2(bestWinSize1:(numFrames-bestWinSize1-bestLagFrames1)) = trajCorrelation(pathR,pathG2,bestWinSize1:(numFrames-bestWinSize1-bestLagFrames1),bestWinSize1,bestLagFrames1);

corr1CutOff = corr1.*(separations1<6*locError);
corr2CutOff = corr2.*(separations2<6*locError);

corrDist1(1:(numFrames-bestLagFrames3)) = probMergedByDist(boxSize,diffConst,dT,separations1((bestLagFrames3+1):numFrames),separations1(1:(numFrames-bestLagFrames3)));
corrDist2(1:(numFrames-bestLagFrames3)) = probMergedByDist(boxSize,diffConst,dT,separations2((bestLagFrames3+1):numFrames),separations2(1:(numFrames-bestLagFrames3)));

%get the sensitivity, specificity, and threshold from the best set
bestSensitivity1 = sensitivity(bestLagFrames1,bestWinSizeIndx1,1);
bestSpecificity1 = specificity(bestLagFrames1,bestWinSizeIndx1,1);
bestThreshold1 = thresholds(bestLagFrames1,bestWinSizeIndx1,1);
bestSensitivity2 = sensitivity(bestLagFrames2,bestWinSizeIndx2,2);
bestSpecificity2 = specificity(bestLagFrames2,bestWinSizeIndx2,2);
bestThreshold2 = thresholds(bestLagFrames2,bestWinSizeIndx2,2);
bestSensitivity3 = sensitivity(bestLagFrames3,bestWinSizeIndx3,3);
bestSpecificity3 = specificity(bestLagFrames3,bestWinSizeIndx3,3);
bestThreshold3 = thresholds(bestLagFrames3,bestWinSizeIndx3,3);


%plot the correlation coefficient over time for the best set of parameters
%and place in the title the important statistics; it also draws the
%threshold determined;
figure
plot(1:numFrames,corr1','.',1:numFrames,corr2','.',[1 numFrames],[bestThreshold1 bestThreshold1],'-')
title(sprintf('Corr Coeff, LagFrames = %g, WinSize = %g, Sensitivity = %g, Specificity = %g',bestLagFrames1,bestWinSize1,bestSensitivity1,bestSpecificity1));

figure
plot(1:numFrames,corr1CutOff','.',1:numFrames,corr2CutOff','.',[1 numFrames],[bestThreshold2 bestThreshold2],'-')
title(sprintf('Corr Coeff with Separation Cutoff, LagFrames = %g, WinSize = %g, Sensitivity = %g, Specificity = %g',bestLagFrames2,bestWinSize2,bestSensitivity2,bestSpecificity2));

figure
plot(1:numFrames,corrDist1','.',1:numFrames,corrDist2','.',[1 numFrames],[bestThreshold3 bestThreshold3],'-')
title(sprintf('Dist Prob, LagFrames = %g, WinSize = %g, Sensitivity = %g, Specificity = %g',bestLagFrames3,bestWinSize3,bestSensitivity3,bestSpecificity3));

%%%%%% the end %%%%%%