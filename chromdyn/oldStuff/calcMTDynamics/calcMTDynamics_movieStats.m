function [statisticsStruct,statisticsCell] = calcMTDynamics_movieStats(distanceData,stateList,verbose,fileName)
%function to calculate the parameters describing the trajectories in calcMTDynamics


%we take the inverses of arrays that can be zero: turn off warnings
warning off MATLAB:divideByZero;

%mean distance. Std is calculated for mean and for the sample (the latter is independent of n!)
[distanceMean, distanceMeanStd, distanceStd] = weightedStats(distanceData(:,1),...
    distanceData(:,2),'s');

%find max and min distances
%we will take the max and the min of the distance, and the mean of the 5 largest
%and 5 smallest distances

sortedDistance = sortrows(distanceData(:,1:2), 1);
nMean = min(5,size(sortedDistance,1));

%max and min
maxDistance = sortedDistance(end,1);
maxDistanceStd = sortedDistance(end,2);
minDistance = sortedDistance(1,1);
minDistanceStd = sortedDistance(1,2);

%mean max and min
[maxDistanceM5,dummy,maxDistanceM5Std] = weightedStats(sortedDistance(end-nMean+1:end,1), sortedDistance(end-nMean+1:end,2));
[minDistanceM5,dummy,minDistanceM5Std] = weightedStats(sortedDistance(1:nMean,1), sortedDistance(1:nMean,2));


%growth statistics

%init vars in case we do not find any growth
growthSpeedMean = [];
growthSpeedStd = [];
invGrowthTimeMean = [];
invGrowthTimeStd = [];
growthTimeTotal = 0;
growthDistanceMean = [];
growthDistanceStd = [];    
growth2ShrinkageNum = [];
invGrowth2ShrinkageTimeMean = [];
invGrowth2ShrinkageTimeStd = [];
growth2UndeterminedIdx = [];
growth2UndeterminedNum = [];
invGrowth2UndeterminedTimeMean = [];
invGrowth2UndeterminedTimeStd = [];
growth2ShrinkageRatio = [];
growth2UndeterminedRatio = [];
growthDistanceTotal = [];

growthIdx = find(stateList(:,2)>0);

if ~isempty(growthIdx)
    %speed
    [growthSpeedMean, growthSpeedStd] = weightedStats(stateList(growthIdx,7),...
        stateList(growthIdx,8),'s');
    %time (inv needed for frequencies, total for %undetermined)
    %rescue and catastrophe are considered poisson processes - hence the
    %catastrophe rate (="half life" of a growing MT) is the inverse mean growth time
    invGrowthTimeMean = 1/mean(stateList(growthIdx,5));
    %std from gauss (sigma = sqrt((df/dm*sigmaM)^2))
    invGrowthTimeStd = std(stateList(growthIdx,5))*invGrowthTimeMean^2;
    growthTimeTotal = sum(stateList(growthIdx,5));
    %distance
    growthDistanceMean = mean(stateList(growthIdx,6));
    growthDistanceStd = std(stateList(growthIdx,6));
    growthDistanceTotal = sum(stateList(growthIdx,6));
    
    %calculate detailed transition frequencies
    
    %if the last block is growth: classify as g2u
    if growthIdx(end) == size(stateList,1)
        growth2UndeterminedIdx = growthIdx(end);
        growthIdx(end) = [];
    end
    
    %calculate g2s: look for <0 blocks after >0 blocks, then calculate
    %stats as above
    growth2ShrinkageIdx = growthIdx(find( stateList(growthIdx+1,2)<0 ));
    growth2ShrinkageNum = length(growth2ShrinkageIdx);
    invGrowth2ShrinkageTimeMean = 1/mean(stateList(growth2ShrinkageIdx,5));
    invGrowth2ShrinkageTimeStd = std(stateList(growth2ShrinkageIdx,5))*invGrowth2ShrinkageTimeMean^2;
    
    %calculate g2u: look for ==0 blocks after >0 blocks, then calculate
    %stats as above. Don't forget the g2uIdx assigned above!
    growth2UndeterminedIdx = [growthIdx(find( stateList(growthIdx+1,2)==0 ));growth2UndeterminedIdx];
    growth2UndeterminedNum = length(growth2UndeterminedIdx);
    invGrowth2UndeterminedTimeMean = 1/mean(stateList(growth2UndeterminedIdx,5));
    invGrowth2UndeterminedTimeStd = std(stateList(growth2UndeterminedIdx,5))*invGrowth2UndeterminedTimeMean^2;
    
    %calculate ratios in %
    growth2ShrinkageRatio = 100*growth2ShrinkageNum/(growth2ShrinkageNum+growth2UndeterminedNum);
    growth2UndeterminedRatio = 100*growth2UndeterminedNum/(growth2ShrinkageNum+growth2UndeterminedNum);        
end




%shrinkage statistics

%init vars
shrinkageSpeedMean = [];
shrinkageSpeedStd = [];
invShrinkageTimeMean = [];
invShrinkageTimeStd = [];
shrinkageTimeTotal = 0;
shrinkageDistanceMean = [];
shrinkageDistanceStd = [];
shrinkage2GrowthNum = [];
invShrinkage2GrowthTimeMean = [];
invShrinkage2GrowthTimeStd = [];
shrinkage2UndeterminedIdx = [];
shrinkage2UndeterminedNum = [];
invShrinkage2UndeterminedTimeMean = [];
invShrinkage2UndeterminedTimeStd = [];
shrinkage2GrowthRatio = [];
shrinkage2UndeterminedRatio = [];
shrinkageDistanceTotal = [];

shrinkageIdx = find(stateList(:,2)<0);

if ~isempty(shrinkageIdx)
    %speed
    [shrinkageSpeedMean, shrinkageSpeedStd] = weightedStats(stateList...
        (shrinkageIdx,7),stateList(shrinkageIdx,8),'s');
    %time
    invShrinkageTimeMean = 1/mean(stateList(shrinkageIdx,5));
    invShrinkageTimeStd = std(stateList(shrinkageIdx,5))*invShrinkageTimeMean^2;
    shrinkageTimeTotal = sum(stateList(shrinkageIdx,5));
    %distance
    shrinkageDistanceMean = mean(stateList(shrinkageIdx,6));
    shrinkageDistanceStd = std(stateList(shrinkageIdx,6));
    shrinkageDistanceTotal = sum(stateList(shrinkageIdx,6));
    
    %calculate detailed transition frequencies
    
    %if the last block is shrinkage: classify as s2u
    if shrinkageIdx(end) == size(stateList,1)
        shrinkage2UndeterminedIdx = shrinkageIdx(end);
        shrinkageIdx(end) = [];
    end
    
    %calculate s2g: look for >0 blocks after <0 blocks, then calculate
    %stats as above
    shrinkage2GrowthIdx = shrinkageIdx(find( stateList(shrinkageIdx+1,2)>0 ));
    shrinkage2GrowthNum = length(shrinkage2GrowthIdx);
    invShrinkage2GrowthTimeMean = 1/mean(stateList(shrinkage2GrowthIdx,5));
    invShrinkage2GrowthTimeStd = std(stateList(shrinkage2GrowthIdx,5))*invShrinkage2GrowthTimeMean^2;
    
    %calculate g2u: look for ==0 blocks after >0 blocks, then calculate
    %stats as above. Don't forget the g2uIdx assigned above!
    shrinkage2UndeterminedIdx = [shrinkageIdx(find( stateList(shrinkageIdx+1,2)==0 ));shrinkage2UndeterminedIdx];
    shrinkage2UndeterminedNum = length(shrinkage2UndeterminedIdx);
    invShrinkage2UndeterminedTimeMean = 1/mean(stateList(shrinkage2UndeterminedIdx,5));
    invShrinkage2UndeterminedTimeStd = std(stateList(shrinkage2UndeterminedIdx,5))*invShrinkage2UndeterminedTimeMean^2;
    
    %calculate ratios in %
    shrinkage2GrowthRatio = 100*shrinkage2GrowthNum/(shrinkage2GrowthNum+shrinkage2UndeterminedNum);
    shrinkage2UndeterminedRatio = 100*shrinkage2UndeterminedNum/(shrinkage2GrowthNum+shrinkage2UndeterminedNum);        
end




%undetermined statistics 

%init vars
undeterminedTimeTotal = 0;
undeterminedDistanceMean = [];
undeterminedDistanceStd = [];
undetermined2GrowthNum = [];
invUndetermined2GrowthTimeMean = [];
invUndetermined2GrowthTimeStd = [];
undetermined2ShrinkageNum = [];
invUndetermined2ShrinkageTimeMean = [];
invUndetermined2ShrinkageTimeStd = [];
undetermined2ShrinkageRatio = [];
undetermined2GrowthRatio = [];

undeterminedIdx = find(stateList(:,2)==0);
if ~isempty(undeterminedIdx)
    %time
    undeterminedTimeTotal = sum(stateList(undeterminedIdx,5));
    %distance
    undeterminedDistanceMean = mean(stateList(undeterminedIdx,6));
    undeterminedDistanceStd = std(stateList(undeterminedIdx,6));
    
    
    %calculate detailed transition frequencies
    
    %if the last block is undetermined: don't consider
    if undeterminedIdx(end) == size(stateList,1)
        undeterminedIdx(end) = [];
    end
    
    %calculate u2g: look for >0 blocks after ==0 blocks, then calculate
    %stats as above
    undetermined2GrowthIdx = undeterminedIdx(find( stateList(undeterminedIdx+1,2)>0 ));
    undetermined2GrowthNum = length(undetermined2GrowthIdx);
    invUndetermined2GrowthTimeMean = 1/mean(stateList(undetermined2GrowthIdx,5));
    invUndetermined2GrowthTimeStd = std(stateList(undetermined2GrowthIdx,5))*invUndetermined2GrowthTimeMean^2;
    
    %calculate u2s: look for <0 blocks after ==0 blocks, then calculate
    %stats as above. 
    undetermined2ShrinkageIdx = undeterminedIdx(find( stateList(undeterminedIdx+1,2)<0 ));
    undetermined2ShrinkageNum = length(undetermined2ShrinkageIdx);
    invUndetermined2ShrinkageTimeMean = 1/mean(stateList(undetermined2ShrinkageIdx,5));
    invUndetermined2ShrinkageTimeStd = std(stateList(undetermined2ShrinkageIdx,5))*invUndetermined2ShrinkageTimeMean^2;
    
    %calculate ratios in %
    undetermined2ShrinkageRatio = 100*undetermined2ShrinkageNum/(undetermined2ShrinkageNum+undetermined2GrowthNum);
    undetermined2GrowthRatio = 100*undetermined2GrowthNum/(undetermined2ShrinkageNum+undetermined2GrowthNum);        
end


%calculate additional stats
undeterminedTimeRatio = 100*undeterminedTimeTotal/(undeterminedTimeTotal+...
    growthTimeTotal+shrinkageTimeTotal);

%write statistics cell. Cat/Resc: as it is considered to be a poisson
%process, the frequencies of catastrophe and rescue are equal to the
%inverse of the time spent in growth/shrinkage
statisticsCell = [{'distanceMean'},{distanceMean},{distanceMeanStd};...
        {'distanceStd'},{distanceStd},{distanceStd/distanceMean};...
        {'minDistance'},{minDistance},{minDistanceStd};...
        {'minDistanceM5'},{minDistanceM5},{minDistanceM5Std};...
        {'maxDistance'},{maxDistance},{maxDistanceStd};...
        {'maxDistanceM5'},{maxDistanceM5},{maxDistanceM5Std};...
        
    {'catFreq'},{invGrowthTimeMean},{invGrowthTimeStd};...
        {'g2sFreq'},{invGrowth2ShrinkageTimeMean},{invGrowth2ShrinkageTimeStd};...
        {'g2sRatio%'},{growth2ShrinkageRatio},{growth2ShrinkageNum};...
        {'g2uFreq'},{invGrowth2UndeterminedTimeMean},{invGrowth2UndeterminedTimeStd};...
        {'g2uRatio%'},{growth2UndeterminedRatio},{growth2UndeterminedNum};...
        
    {'resFreq'},{invShrinkageTimeMean},{invShrinkageTimeStd};...
        {'s2gFreq'},{invShrinkage2GrowthTimeMean},{invShrinkage2GrowthTimeStd};...
        {'s2gRatio%'},{shrinkage2GrowthRatio},{shrinkage2GrowthNum};...
        {'s2uFreq'},{invShrinkage2UndeterminedTimeMean},{invShrinkage2UndeterminedTimeStd};...
        {'s2uRatio%'},{shrinkage2UndeterminedRatio},{shrinkage2UndeterminedNum};...
        
    {'u2gFreq'},{invUndetermined2GrowthTimeMean},{invUndetermined2GrowthTimeStd};...
        {'u2gRatio%'},{undetermined2GrowthRatio},{undetermined2GrowthNum};...
        {'u2sFreq'},{invUndetermined2ShrinkageTimeMean},{invUndetermined2ShrinkageTimeStd};...
        {'u2sRatio%'},{undetermined2ShrinkageRatio},{undetermined2ShrinkageNum};...
        
    {'growthSpeed'},{growthSpeedMean},{growthSpeedStd};...
        {'shrinkageSpeed'},{shrinkageSpeedMean},{shrinkageSpeedStd};...
        
    {'growthDistance'},{growthDistanceMean},{growthDistanceStd};...
        {'growthDistanceTotal'},{growthDistanceTotal},{[]};...
        {'shrinkageDistance'},{shrinkageDistanceMean},{shrinkageDistanceStd};...
        {'shrinkageDistanceTotal'},{shrinkageDistanceTotal},{[]};...
        {'undetDistance'},{undeterminedDistanceMean},{undeterminedDistanceStd};...
        
    {'undeterminedTime%'},{undeterminedTimeRatio},{[]}];

statisticsStruct = struct('distanceMean',[distanceMean,distanceMeanStd],...
    'distanceStd',[distanceStd, distanceStd/distanceMean],...
        'minDistance', [minDistance , minDistanceStd],...
         'minDistanceM5' , [minDistanceM5 , minDistanceM5Std],...
         'maxDistance' , [maxDistance , maxDistanceStd],...
         'maxDistanceM5' , [maxDistanceM5 , maxDistanceM5Std],...
     'catFreq' , [invGrowthTimeMean , invGrowthTimeStd],...
         'g2sFreq' , [invGrowth2ShrinkageTimeMean , invGrowth2ShrinkageTimeStd],...
         'g2sRatioPrc' , [growth2ShrinkageRatio , growth2ShrinkageNum],...
         'g2uFreq' , [invGrowth2UndeterminedTimeMean , invGrowth2UndeterminedTimeStd],...
         'g2uRatioPrc' , [growth2UndeterminedRatio , growth2UndeterminedNum],...
     'resFreq' , [invShrinkageTimeMean , invShrinkageTimeStd],...
         's2gFreq' , [invShrinkage2GrowthTimeMean , invShrinkage2GrowthTimeStd],...
         's2gRatioPrc' , [shrinkage2GrowthRatio , shrinkage2GrowthNum],...
         's2uFreq' , [invShrinkage2UndeterminedTimeMean , invShrinkage2UndeterminedTimeStd],...
         's2uRatioPrc' , [shrinkage2UndeterminedRatio , shrinkage2UndeterminedNum],...
     'u2gFreq' , [invUndetermined2GrowthTimeMean , invUndetermined2GrowthTimeStd],...
         'u2gRatioPrc' , [undetermined2GrowthRatio , undetermined2GrowthNum],...
         'u2sFreq' , [invUndetermined2ShrinkageTimeMean , invUndetermined2ShrinkageTimeStd],...
         'u2sRatioPrc' , [undetermined2ShrinkageRatio , undetermined2ShrinkageNum],...
     'growthSpeed' , [growthSpeedMean , growthSpeedStd],...
         'shrinkageSpeed' , [shrinkageSpeedMean , shrinkageSpeedStd],...
     'growthDistance' , [growthDistanceMean , growthDistanceStd],...
         'growthDistanceTotal' , [growthDistanceTotal],...
         'shrinkageDistance' , [shrinkageDistanceMean , shrinkageDistanceStd],...
         'shrinkageDistanceTotal' , [shrinkageDistanceTotal],...
         'undetDistance' , [undeterminedDistanceMean , undeterminedDistanceStd],...
     'undeterminedTimePrc' , [undeterminedTimeRatio]);

if verbose
    disp(' ');
    disp(['Statistics for: ',fileName,' :']);
    disp(statisticsCell);
end

%don't forget to turn the warnings back on
warning on MATLAB:divideByZero