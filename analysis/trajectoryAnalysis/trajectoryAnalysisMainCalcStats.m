function [statisticsStruct, distributionStruct, clusterStruct] = trajectoryAnalysisMainCalcStats(dataListG,distance,time,verbose,fileName,constants,individualMeanDistances)
%TRAJECTORYANALYSISMAINCALCSTATS calculates statistics for trajectoryAnalysis
%
% OUTPUT statistics structure with fields
%     'ap2tpFreq' ,             catastrophe frequency
%     'tp2apFreq' ,             rescue frequency
%     'antipolewardSpeed' ,         growth speed
%     'polewardSpeed' ,        shrinkage speed
%           for individual stats
%     'distanceMean',             std is std of sample!
%           for overall stats
%     'distanceMean' mean +/- std
%     'distanceStd' mean of sample stds +/- std
%     'minDistance',              
%     'minDistanceM5' ,           mean of 5 smallest distances
%     'maxDistance' ,             
%     'maxDistanceM5' ,           mean of 5 largest distances
%     'antipolewardDistance' ,      mean distance per growth event
%     'antipolewardDistanceTotal' , 
%     'polewardDistance' ,     mean distance per shrinkage event
%     'polewardDistanceTotal', 
%     'undeterminedDistance' ,    
%     'antipolewardTime' ,          
%     'polewardTime' ,         
%     'pauseTime' ,               
%     'undeterminedTime' ,        
%     'deletedTime' ,             
%     'nTimepoints', 
%
%      and distribution structure with distributions from contHisto for
%      both speeds and the distance
%
%     distributionStruct.antipolewardSpeedDistribution = [growthSpeedDistX,growthSpeedDistY];
%     distributionStruct.polewardSpeedDistribution = [shrinkageSpeedDistX,shrinkageSpeedDistY];
%     distributionStruct.distanceDistribution = [distanceDistX,distanceDistY];
%
%      and clusterStruct:
%       clusterStruct.numberOfDistributions = # of clusters
%       clusterStruct.centersOfDistributions = #oD means
%       clusterStruct.?
%
% c: 1/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distributionStruct = [];
clusterStruct      = [];
statisticsStruct   = [];

%we take the inverses of arrays that can be zero: turn off warnings
warning off MATLAB:divideByZero;

%------- INIT VARIABLES

%GROWTH
growthSpeedMean = [];
growthSpeedStd = [];
growthSpeedMeanNW = [];
growthSpeedStdNW = [];
growthSpeedMeanPW = [];
growthSpeedStdPW = [];
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
growthTimeRatio = [];

%SHRINKAGE
shrinkageSpeedMean = [];
shrinkageSpeedStd = [];
shrinkageSpeedMeanNW = [];
shrinkageSpeedStdNW = [];
shrinkageSpeedMeanPW = [];
shrinkageSpeedStdPW = [];
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
shrinkageTimeRatio = [];

%PAUSE
pauseTimeTotal = 0;
pauseTimeRatio = [];

%UNDETERMINED
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

%DELETED
deletedTimeTotal = 0;



%----END INIT VARIABLES

%--- OVERALL STATISTICS
if nargin < 7 | isempty(individualMeanDistances)
    % we're calculating individual data. Therefore: Mean+/-stdOfMean, std
    % of sample
    
    % mean distance. Std is calculated for the sample (it is independent of n!)
    [distanceMean, distanceMeanStd, distanceStd] = weightedStats(distance(:,1),...
        distance(:,2),'s');
    distanceStdStd = NaN;
    
    % assign average timepoints
    avgTimePoints = NaN;
else
    % we're calculating overall stats - take mean of distances and stds
    [distanceMean, distanceMeanStd] = weightedStats(distance(:,1),...
        distance(:,2),'s');
    % mean of std is mean of individual std, weighted with number of timepoints
    [distanceStd, distanceStdStd] = weightedStats(individualMeanDistances(:,2),individualMeanDistances(:,3),'w');
    
    % assign average timepoints
    avgTimePoints = mean(individualMeanDistances(:,3));
end

%find max and min distances
%we will take the max and the min of the distance, and the mean of the 5 largest
%and 5 smallest distances

sortedDistance = sortrows(distance, [1,2]);
nMean = min(5,size(sortedDistance,1));

%max and min
maxDistance = sortedDistance(end,1);
maxDistanceStd = sortedDistance(end,2);
minDistance = sortedDistance(1,1);
minDistanceStd = sortedDistance(1,2);

%mean max and min
[maxDistanceM5,dummy,maxDistanceM5Std] = weightedStats(sortedDistance(end-nMean+1:end,1), sortedDistance(end-nMean+1:end,2));
[minDistanceM5,dummy,minDistanceM5Std] = weightedStats(sortedDistance(1:nMean,1), sortedDistance(1:nMean,2));

%number of valid time points in this movie
nTimePoints = size(distance,1);

%--- END OVERALL DISTANCE


%------- STATISTICS

dlgSize = size(dataListG);


%CALC INDICES
%we have to find two kinds of indices for g&s: one that points to the
%single events, a second that indicates groups of events (i.e. ggpg as one g)

%---growth---
growthIdx = find(dataListG(:,3)==1);

%growth belonging to the same phase are either direcly one after the other
%or separated by a pause -> make list of ones and zeros and find the groups


%try to find groups
gGroups = findGroups(dataListG(:,3),1);

%search where the groups are separated by just one other unit
%then check whether these are pauses
%then look through the list and fuse entries

gGroupsDiff  = gGroups(2:end,2)-gGroups(1:end-1,1);
gGroupsDiff2 = find(gGroupsDiff==2);

%lookfor pauses (gpgStates points to gGroupsDiff2, which points to gGroups, which points to dataListG)
gpgStates = find(dataListG(gGroups(gGroupsDiff2,2)+1,3)==3);

%take again list for growthIdx
stateListGrowth = dataListG(:,3);
%change pauses to growth
stateListGrowth(gGroups(gGroupsDiff2(gpgStates),2)+1) = 1;

%find groups again
growthGroups = findGroups(stateListGrowth,1);


%---shrinkage---
shrinkageIdx = find(dataListG(:,3)==2);

%shrinkage belonging to the same phase are either direcly one after the other
%or separated by a pause -> make list of ones and zeros and find the groups

%try to find groups
sGroups = findGroups(dataListG(:,3),2);

%search where the groups are separated by just one other unit
%then check whether these are pauses
%then look through the list and fuse entries

sGroupsDiff  = sGroups(2:end,2)-sGroups(1:end-1,1);
sGroupsDiff2 = find(sGroupsDiff==2);

%lookfor pauses (spsStates points to sGroupsDiff2, which points to sGroups, which points to dataListG)
spsStates = find(dataListG(sGroups(sGroupsDiff2,2)+1,3)==3);

%take again list for shrinkageIdx
stateListShrinkage = dataListG(:,3);
%change pauses to shrinkage
stateListShrinkage(sGroups(sGroupsDiff2(spsStates),2)+1) = 2;

%find groups again
shrinkageGroups = findGroups(stateListShrinkage,2);

%pause & undetermined & deleted
pauseIdx = find(dataListG(:,3)==3);
undeterminedIdx = find(dataListG(:,3)<=0 & dataListG(:,3) >= -constants.MAXDELETED);

% try to find undetermined groups
undeterminedGroups = findGroups(dataListG(:,3),[0:-1:constants.MAXDELETED]);

% deleted index. don't count -99
deletedIdx = find(dataListG(:,3) < -constants.MAXDELETED & dataListG(:,3)>-99);

sigmaMax  = max(dataListG(:,5:6),[],2);
sigmaMax2 = sigmaMax.^2; 

% 

%growth statistics

if ~isempty(growthIdx)
    %speed -- weight sigma and time
    
    % old sigma - if used, change distributions below!!
    %indGrowthSpeedSigma = sigmaMax(growthIdx)./sqrt(dataListG(growthIdx,7));
    % new sigma: less speed dep., better time dependence. (use relative sigma, weight with time)
    % We could multiply by growthSpeedMean (we have to to plot distribution of individual
    % speeds!), but it does not change anything here: the Std is calculated
    % for the difference between values, not for the individual stds!
    
    % 2/23: Do not use sigma-weighting any more: too unstable!
    % indGrowthSpeedSigma = sigmaMax(growthIdx)./(abs(dataListG(growthIdx,4)).*sqrt(dataListG(growthIdx,7)));
    %     [growthSpeedMean, growthSpeedStd] = weightedStats(dataListG(growthIdx,4),...
    %         indGrowthSpeedSigma,'s');
    %     [growthSpeedMeanPW, growthSpeedStdPW] = weightedStats(dataListG(growthIdx,4),...
    %         sigmaMax(growthIdx),'s');
    
    % use normal mean, weigh only with time
    indGrowthSpeedSigma = 1./dataListG(growthIdx,7);
    [growthSpeedMeanNW,growthSpeedStdNW] = weightedStats(dataListG(growthIdx,4)*60,dataListG(growthIdx,7),'w');
    
    %transform to um/min
    %     growthSpeedMean = growthSpeedMean*60;
    %     growthSpeedStd  = growthSpeedStd *60;
    %     growthSpeedMeanPW = growthSpeedMeanPW*60;
    %     growthSpeedStdPW  = growthSpeedStdPW *60;
    
    %time (inv needed for frequencies, total for %undetermined)
    %rescue and catastrophe are considered poisson processes - hence the
    %catastrophe rate (="half life" of a growing MT) is the inverse mean growth time
    %time is calculated with groups
    deltaT = time(dataListG(growthGroups(:,2),2),1) - time(dataListG(growthGroups(:,1),1),1);
    
    invGrowthTimeMean = 1/mean(deltaT);
    %std from gauss (sigma = sqrt((df/dm*sigmaM)^2))
    invGrowthTimeStd = (std(deltaT)/sqrt(size(growthGroups,1)))*invGrowthTimeMean^2;
    growthTimeTotal = sum(deltaT);
    
    %distance - again, use groups. unfortunately, we need the loop
    groupDist = zeros(size(growthGroups));
    for i = 1:size(growthGroups,1)
        groupDist(i,1) = sum(dataListG(growthGroups(i,1):growthGroups(i,2),9));
        groupDist(i,2) = sqrt(sum(sigmaMax2(growthGroups(i,1):growthGroups(i,2))));
    end
    [growthDistanceMean, growthDistanceStd] = weightedStats(groupDist(:,1),groupDist(:,2),'s');
    growthDistanceTotal = sum(groupDist(:,1));
    
    %     %calculate detailed transition frequencies
    %     
    %     %if the last block is growth: classify as g2u
    %     if growthIdx(end) == size(stateList,1)
    %         growth2UndeterminedIdx = growthIdx(end);
    %         growthIdx(end) = [];
    %     end
    %     
    %     %calculate g2s: look for <0 blocks after >0 blocks, then calculate
    %     %stats as above
    %     growth2ShrinkageIdx = growthIdx(find( stateList(growthIdx+1,2)<0 ));
    %     growth2ShrinkageNum = length(growth2ShrinkageIdx);
    %     invGrowth2ShrinkageTimeMean = 1/mean(stateList(growth2ShrinkageIdx,5));
    %     invGrowth2ShrinkageTimeStd = std(stateList(growth2ShrinkageIdx,5))*invGrowth2ShrinkageTimeMean^2;
    %     
    %     %calculate g2u: look for ==0 blocks after >0 blocks, then calculate
    %     %stats as above. Don't forget the g2uIdx assigned above!
    %     growth2UndeterminedIdx = [growthIdx(find( stateList(growthIdx+1,2)==0 ));growth2UndeterminedIdx];
    %     growth2UndeterminedNum = length(growth2UndeterminedIdx);
    %     invGrowth2UndeterminedTimeMean = 1/mean(stateList(growth2UndeterminedIdx,5));
    %     invGrowth2UndeterminedTimeStd = std(stateList(growth2UndeterminedIdx,5))*invGrowth2UndeterminedTimeMean^2;
    %     
    %     %calculate ratios in %
    %     growth2ShrinkageRatio = 100*growth2ShrinkageNum/(growth2ShrinkageNum+growth2UndeterminedNum);
    %     growth2UndeterminedRatio = 100*growth2UndeterminedNum/(growth2ShrinkageNum+growth2UndeterminedNum);        
end
% 
% 
% 
% 
% %shrinkage statistics
if ~isempty(shrinkageIdx)
    
    % old sigma - if used, change distributions below!!
    %indShrinkageSpeedSigma = sigmaMax(shrinkageIdx)./sqrt(dataListG(shrinkageIdx,7));
    % new sigma: less speed dep., better time dependence. (use relative sigma, weight with time)
    % We could multiply by shrinkageSpeedMean (we have to to plot distribution of individual
    % speeds!), but it does not change anything here: the Std is calculated
    % for the difference between values, not for the individual stds!
    
    % 2/23: Do not use sigma-weighting any more: too unstable!
    % indShrinkageSpeedSigma = sigmaMax(shrinkageIdx)./(abs(dataListG(shrinkageIdx,4)).*sqrt(dataListG(shrinkageIdx,7)));
    % 
    %     [shrinkageSpeedMean, shrinkageSpeedStd] = weightedStats(dataListG(shrinkageIdx,4),...
    %         indShrinkageSpeedSigma,'s');
    %     [shrinkageSpeedMeanPW, shrinkageSpeedStdPW] = weightedStats(dataListG(shrinkageIdx,4),...
    %         sigmaMax(shrinkageIdx),'s');
    
    %non-weighted stats
    indShrinkageSpeedSigma = 1./dataListG(shrinkageIdx,7);
    [shrinkageSpeedMeanNW,shrinkageSpeedStdNW] = weightedStats(dataListG(shrinkageIdx,4)*60,dataListG(shrinkageIdx,7),'w');
    
    %transform to um/min
    %     shrinkageSpeedMean = shrinkageSpeedMean*60;
    %     shrinkageSpeedStd  = shrinkageSpeedStd *60;
    %     shrinkageSpeedMeanPW = shrinkageSpeedMeanPW*60;
    %     shrinkageSpeedStdPW  = shrinkageSpeedStdPW *60;
    
    %time (inv needed for frequencies, total for %undetermined)
    %rescue and catastrophe are considered poisson processes - hence the
    %catastrophe rate (="half life" of a growing MT) is the inverse mean shrinkage time
    %time is calculated with groups
    deltaT = time(dataListG(shrinkageGroups(:,2),2),1) - time(dataListG(shrinkageGroups(:,1),1),1);
    
    invShrinkageTimeMean = 1/mean(deltaT);
    %std from gauss (sigma = sqrt((df/dm*sigmaM)^2))
    invShrinkageTimeStd = (std(deltaT)/sqrt(size(shrinkageGroups,1)))*invShrinkageTimeMean^2;
    shrinkageTimeTotal = sum(deltaT);
    
    %distance - again, use groups. unfortunately, we need the loop
    groupDist = zeros(size(shrinkageGroups));
    for i = 1:size(shrinkageGroups,1)
        groupDist(i,1) = sum(dataListG(shrinkageGroups(i,1):shrinkageGroups(i,2),9));
        groupDist(i,2) = sqrt(sum(sigmaMax2(shrinkageGroups(i,1):shrinkageGroups(i,2))));
    end
    [shrinkageDistanceMean, shrinkageDistanceStd] = weightedStats(groupDist(:,1),groupDist(:,2),'s');
    shrinkageDistanceTotal = sum(groupDist(:,1));
end


%UNDETERMINED
if ~isempty(undeterminedIdx)
    %time
    undeterminedTimeTotal = sum(dataListG(undeterminedIdx,7));
    %distance - use groups
    for i = 1:size(undeterminedGroups,1)
        groupDist(i,1) = sum(dataListG(undeterminedGroups(i,1):undeterminedGroups(i,2),9));
        groupDist(i,2) = sum(dataListG([undeterminedGroups(i,1),undeterminedGroups(i,2)],10));
    end
    %use absolute values for undetermined distance!
    [undeterminedDistanceMean,undeterminedDistanceStd] = weightedStats(abs(groupDist(:,1)),groupDist(:,2),'s');
end

%PAUSE
if ~isempty(pauseIdx)
    %time
    pauseTimeTotal = sum(dataListG(pauseIdx,7));
    
end

%DELETED
if ~isempty(deletedIdx)
    deletedTimeTotal = sum(dataListG(deletedIdx,7));
end



%calculate additional stats
divisor = (undeterminedTimeTotal+...
    growthTimeTotal+shrinkageTimeTotal+pauseTimeTotal);
undeterminedTimeRatio = 100*undeterminedTimeTotal/divisor;
growthTimeRatio = 100*growthTimeTotal/divisor;
shrinkageTimeRatio = 100*shrinkageTimeTotal/divisor;
pauseTimeRatio = 100*pauseTimeTotal/divisor;

%write statistics structure
% 
statisticsStruct = struct(...
    'ap2tpFreq__cat' ,            [invGrowthTimeMean , invGrowthTimeStd],...
    'tp2apFreq__res' ,            [invShrinkageTimeMean , invShrinkageTimeStd],...
    'antipolewardSpeed' ,         [growthSpeedMeanNW , growthSpeedStdNW],...
    'polewardSpeed' ,             [shrinkageSpeedMeanNW , shrinkageSpeedStdNW],...
    'distanceMean',               [distanceMean,distanceMeanStd],...
    'distanceStd',                [distanceStd,distanceStdStd],...
    'minDistance',                [minDistance , minDistanceStd],...
    'minDistanceM5' ,             [minDistanceM5 , minDistanceM5Std],...
    'maxDistance' ,               [maxDistance , maxDistanceStd],...
    'maxDistanceM5' ,             [maxDistanceM5 , maxDistanceM5Std],...
    'antipolewardDistance' ,      [growthDistanceMean , growthDistanceStd],...
    'polewardDistance' ,          [shrinkageDistanceMean , shrinkageDistanceStd],...
    'undeterminedDistance' ,      [undeterminedDistanceMean , undeterminedDistanceStd],...
    'antipolewardTime' ,          [growthTimeTotal,growthTimeRatio],...
    'polewardTime' ,              [shrinkageTimeTotal,shrinkageTimeRatio],...
    'pauseTime' ,                 [pauseTimeTotal,pauseTimeRatio],...
    'undeterminedTime' ,          [undeterminedTimeTotal,undeterminedTimeRatio],...
    'deletedTime' ,               [deletedTimeTotal],...
    'nTimepoints',                [nTimePoints,avgTimePoints]);  

%    'distanceStd',[distanceStd, distanceStd/distanceMean],...
%          'g2sFreq' , [invGrowth2ShrinkageTimeMean , invGrowth2ShrinkageTimeStd],...
%          'g2sRatioPrc' , [growth2ShrinkageRatio , growth2ShrinkageNum],...
%          'g2uFreq' , [invGrowth2UndeterminedTimeMean , invGrowth2UndeterminedTimeStd],...
%          'g2uRatioPrc' , [growth2UndeterminedRatio , growth2UndeterminedNum],...
%      
%          's2gFreq' , [invShrinkage2GrowthTimeMean , invShrinkage2GrowthTimeStd],...
%          's2gRatioPrc' , [shrinkage2GrowthRatio , shrinkage2GrowthNum],...
%          's2uFreq' , [invShrinkage2UndeterminedTimeMean , invShrinkage2UndeterminedTimeStd],...
%          's2uRatioPrc' , [shrinkage2UndeterminedRatio , shrinkage2UndeterminedNum],...
%      'u2gFreq' , [invUndetermined2GrowthTimeMean , invUndetermined2GrowthTimeStd],...
%          'u2gRatioPrc' , [undetermined2GrowthRatio , undetermined2GrowthNum],...
%          'u2sFreq' , [invUndetermined2ShrinkageTimeMean , invUndetermined2ShrinkageTimeStd],...
%          'u2sRatioPrc' , [undetermined2ShrinkageRatio , undetermined2ShrinkageNum],...
%     'antipolewardDistanceTotal' , [growthDistanceTotal],...
%    'polewardDistanceTotal',      [shrinkageDistanceTotal],...
% 'growthSpeed' , [growthSpeedMean , growthSpeedStd],...
%     'shrinkageSpeed' ,
if any(verbose == 1)
    titles = fieldnames(statisticsStruct);
    values = struct2cell(statisticsStruct);
    disp(' ');
    disp(['Statistics for  ',fileName,' :']);
    statstr = '';
    for nTit = 1:length(titles)
        statstr = [statstr,sprintf(['\n',titles{nTit}])];
        currVal = values{nTit};
        for nVal = 1:length(currVal)
            statstr = [statstr,sprintf('\t%f',currVal(nVal))];
        end
    end
    disp(statstr)
end


%--------------CALCULATE DISTRIBUTIONS
if nargout > 1
    % initialize to avoid problems with empty indexLists
    distributionStruct.antipolewardSpeedDistribution = [];
    distributionStruct.polewardSpeedDistribution = [];
    distributionStruct.distanceDistribution = [];
    
    
    
    
    if ~isempty(growthIdx)
        
        % calculate distributions
        [growthSpeedDistY,growthSpeedDistX] = contHisto([60*dataListG(growthIdx,4),...
                indGrowthSpeedSigma],'norm',1,0); %indGrowthSpeedSigma*growthSpeedMeanNW
        
        
        
        
        
    else
        growthSpeedDistY = [];
        growthSpeedDistX = [];
    end
    
    if ~isempty(shrinkageIdx)
        
        % calculate distributions
        [shrinkageSpeedDistY,shrinkageSpeedDistX] = contHisto([60*dataListG(shrinkageIdx,4),...
                indShrinkageSpeedSigma],'norm',1,0); % indShrinkageSpeedSigma*abs(shrinkageSpeedMeanNW)
        
        
    else
        shrinkageSpeedDistY = [];
        shrinkageSpeedDistX = [];
    end
    
    % there will always be distance...
    [distanceDistY,distanceDistX] = contHisto(distance,'norm',1,0);
    
    % plot without detail. Would cost too much memory
    switch any(verbose == 3)*(~isempty(growthIdx)+2*~isempty(shrinkageIdx))
        case 1 % only growth
            figure('Name','speed distribution'),plot(growthSpeedDistX,growthSpeedDistY/max(growthSpeedDistY),'-g'...
                );
        case 2 % only shrinkage
            figure('Name','speed distribution'),plot(...
                -shrinkageSpeedDistX,shrinkageSpeedDistY/max(shrinkageSpeedDistY),'-r');
        case 3 % both growth and shrinkage
            figure('Name','speed distribution'),plot(growthSpeedDistX,growthSpeedDistY/max(growthSpeedDistY),'-g',...
                -shrinkageSpeedDistX,shrinkageSpeedDistY/max(shrinkageSpeedDistY),'-r');
        otherwise % don't plot
    end
    if any(verbose == 3) % don't forget distance!!
        figure('Name','distance distribution'),area(distanceDistX,distanceDistY);
    end
    
    
    distributionStruct.antipolewardSpeedDistribution = [growthSpeedDistX,growthSpeedDistY];
    distributionStruct.polewardSpeedDistribution = [shrinkageSpeedDistX,shrinkageSpeedDistY];
    distributionStruct.distanceDistribution = [distanceDistX,distanceDistY];
end
%-------------end calc dist
% if requested and possible, cluster speeds
if nargout > 2 
    % init cluster, make sure we get the right size (i.e. do we need 1 or 5 fields?)
    initRange = [1:(constants.CLUSTERIND-1)*constants.CLUSTERMAX+1];
    clusterStruct(initRange) = struct(...
        'antipolewardSpeed',[],'polewardSpeed',[],'apIndividual',[],...
        'tpIndividual',[],'apNum',[],'tpNum',[]);
    
    if length(growthIdx) > constants.CLUSTERMAX
        
        % ======== Cluster ============            
        [clusterStruct(1).antipolewardSpeed,clusterStruct(1).apNum,clusterStruct(initRange).apIndividual] =...
            trajectoryAnalysisMainCalcStatsCluster(dataListG,growthIdx,constants);
        
    end
    
    % if requested and possible, cluster speeds
    if length(shrinkageIdx) > constants.CLUSTERMAX
        
        % ======== Cluster ============
        [clusterStruct(1).polewardSpeed,clusterStruct(1).tpNum,clusterStruct(initRange).tpIndividual] =...
            trajectoryAnalysisMainCalcStatsCluster(dataListG,shrinkageIdx,constants);
    end
    
end

%don't forget to turn the warnings back on
warning on MATLAB:divideByZero