function [statisticsStruct, addedStats, distributionStruct, clusterStruct] = trajectoryAnalysisMainCalcStats(dataListG,distance,time,verbose,fileName,constants,individualMeanDistances)
%TRAJECTORYANALYSISMAINCALCSTATS calculates statistics for trajectoryAnalysis
%
% OUTPUT statistics structure with fields
%        fields of statistics-structure 
%           SEM: standard error of the mean
%           STD: standard deviation of the sample (SEM*sqrt(n))
%
%     'ap2tpFreq__cat' ,        catastrophe frequency (m,sem,n; [s^-1])
%     'tp2apFreq__res' ,        rescue frequency      (m,sem,n; [s^-1])
%     'antipolewardSpeed' ,     growth speed          (m,sem,n; [um/min])
%     'polewardSpeed' ,         shrinkage speed       (m,sem,n; [um/min])
%     'distanceMean',           mean spb-cen distance (m,sem; [um])
%     'distanceStd'             std of distance       (m,sem; [um])
%           for individual statistics, there can be no sem
%     'minDistance',            global minimum distance (m,std; [um])
%     'minDistanceM5' ,         mean of 5 smallest distances (m,sem; [um])
%     'maxDistance' ,           global maximum distance (m,std; [um])  
%     'maxDistanceM5' ,         mean of 5 largest distances (m,sem; [um])
%     'pauseNumber',            number of pause events
%     'avgApDistance' ,         mean distance per growth event (m,sem; [um])
%     'avgTpDistance' ,         mean distance per shrinkage event (m,sem; [um])
%     'avgUndetDistance' ,      avg of absolute distance in undet. intervals (m,sem; [um])
%     'antipolewardTime' ,      total AP time [s]; % of total traj. time     
%     'polewardTime' ,          total TP time [s]; % of total traj. time
%     'pauseTime' ,             total pause time [s]; % of total traj. time 
%     'undeterminedTime' ,      total undet. time [s]; % of total traj. time 
%     'deletedTime' ,           total not analyzed time [s] - not counting
%                               deletion at the end of the trajectory
%     'nTimepoints',            number of total timepoints; avg per trajectory
%        
%
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
statisticsStruct   = [];

%we take the inverses of arrays that can be zero: turn off warnings
warningState = warning;
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
growthNumber = 0;
growthGroupNumber = 0;

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
shrinkageNumber = 0;
shrinkageGroupNumber = 0;

%PAUSE
pauseTimeTotal = 0;
pauseTimeRatio = [];
pauseNumber = 0;

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
[growthIdx,growthGroups] = trajectoryAnalysisMainCalcStatsFindGroups(dataListG,1);

%---shrinkage---
[shrinkageIdx,shrinkageGroups] = trajectoryAnalysisMainCalcStatsFindGroups(dataListG,2);

%------------

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
%====================
% growth statistics
%====================

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
    growthGroupsDeltaT = time(dataListG(growthGroups(:,2),2),1) - time(dataListG(growthGroups(:,1),1),1);
    
    invGrowthTimeMean = 1/mean(growthGroupsDeltaT);
    %std from gauss (sigma = sqrt((df/dm*sigmaM)^2))
    invGrowthTimeStd = (std(growthGroupsDeltaT)/sqrt(size(growthGroups,1)))*invGrowthTimeMean^2;
    growthTimeTotal = sum(growthGroupsDeltaT);
    
    %distance - again, use groups. unfortunately, we need the loop
    groupDist = zeros(size(growthGroups));
    for i = 1:size(growthGroups,1)
        groupDist(i,1) = sum(dataListG(growthGroups(i,1):growthGroups(i,2),9));
        groupDist(i,2) = sqrt(sum(dataListG(growthGroups(i,1):growthGroups(i,2),10)));
    end
    [growthDistanceMean, growthDistanceStd] = weightedStats(groupDist(:,1),groupDist(:,2),'s');
    growthDistanceTotal = sum(groupDist(:,1));
    
    % calculate numbers
    growthNumber = length(growthIdx);
    growthGroupNumber = size(growthGroups,1);
end
% 
% 
% 
% 
%====================
% %shrinkage statistics
%====================

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
    shrinkageGroupsDeltaT = time(dataListG(shrinkageGroups(:,2),2),1) - time(dataListG(shrinkageGroups(:,1),1),1);
    
    invShrinkageTimeMean = 1/mean(shrinkageGroupsDeltaT);
    %std from gauss (sigma = sqrt((df/dm*sigmaM)^2))
    invShrinkageTimeStd = (std(shrinkageGroupsDeltaT)/sqrt(size(shrinkageGroups,1)))*invShrinkageTimeMean^2;
    shrinkageTimeTotal = sum(shrinkageGroupsDeltaT);
    
    %distance - again, use groups. unfortunately, we need the loop
    groupDist = zeros(size(shrinkageGroups));
    for i = 1:size(shrinkageGroups,1)
        groupDist(i,1) = sum(dataListG(shrinkageGroups(i,1):shrinkageGroups(i,2),9));
        groupDist(i,2) = sqrt(sum(dataListG(shrinkageGroups(i,1):shrinkageGroups(i,2),10)));
    end
    [shrinkageDistanceMean, shrinkageDistanceStd] = weightedStats(groupDist(:,1),groupDist(:,2),'s');
    shrinkageDistanceTotal = sum(groupDist(:,1));
    
    % calculate numbers
    shrinkageNumber = length(shrinkageIdx);
    shrinkageGroupNumber = size(shrinkageGroups,1);
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
    % number
    pauseNumber = length(pauseIdx);
end

%DELETED
if ~isempty(deletedIdx)
    deletedTimeTotal = sum(dataListG(deletedIdx,7));
end



%====================
% calculate additional stats
%====================
divisor = (undeterminedTimeTotal+...
    growthTimeTotal+shrinkageTimeTotal+pauseTimeTotal);
undeterminedTimeRatio = 100*undeterminedTimeTotal/divisor;
growthTimeRatio = 100*growthTimeTotal/divisor;
shrinkageTimeRatio = 100*shrinkageTimeTotal/divisor;
pauseTimeRatio = 100*pauseTimeTotal/divisor;

% J and L (see Verde et al, 1992)
if ~isempty(growthSpeedMeanNW) & ~isempty(shrinkageSpeedMeanNW) & ~isempty(invGrowthTimeMean) & ~isempty(invShrinkageTimeMean)
   
    fres = invShrinkageTimeMean*60; % transform into min^-1
    fcat = invGrowthTimeMean*60;    % same
    
    % J = [vg*fres-vs*fcat]/[fcat+fres]
    jVerdeEtAl = (growthSpeedMeanNW * fres - abs(shrinkageSpeedMeanNW) * fcat)/(fcat + fres);
    % L = [vs*vg]/[vs*fcat-vg*fres]
    lVerdeEtAl = (growthSpeedMeanNW * abs(shrinkageSpeedMeanNW)) / (abs(shrinkageSpeedMeanNW) * fcat - growthSpeedMeanNW * fres);
    
else
    jVerdeEtAl = NaN;
    lVerdeEtAl = NaN;
end

% FRAP
% - here is where it would be -



%write statistics structure
% 
statisticsStruct = struct(...
    'ap2tpFreq__cat' ,            [invGrowthTimeMean , invGrowthTimeStd, growthGroupNumber],...
    'tp2apFreq__res' ,            [invShrinkageTimeMean , invShrinkageTimeStd, shrinkageGroupNumber],...
    'antipolewardSpeed' ,         [growthSpeedMeanNW , growthSpeedStdNW, growthNumber],...
    'polewardSpeed' ,             [shrinkageSpeedMeanNW , shrinkageSpeedStdNW, shrinkageNumber],...
    'distanceMean',               [distanceMean,distanceMeanStd],...
    'distanceStd',                [distanceStd,distanceStdStd],...
    'minDistance',                [minDistance , minDistanceStd],...
    'minDistanceM5' ,             [minDistanceM5 , minDistanceM5Std],...
    'maxDistance' ,               [maxDistance , maxDistanceStd],...
    'maxDistanceM5' ,             [maxDistanceM5 , maxDistanceM5Std],...
    'pauseNumber',                [pauseNumber],...
    'avgApDistance' ,             [growthDistanceMean , growthDistanceStd],...
    'avgTpDistance' ,             [shrinkageDistanceMean , shrinkageDistanceStd],...
    'avgUndetDistance' ,          [undeterminedDistanceMean , undeterminedDistanceStd],...
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
if nargout > 3 
    % init cluster, make sure we get the right size (i.e. do we need 1 or 5 fields?)
    initRange = [1:(constants.CLUSTERIND)*constants.CLUSTERMAX+1];
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

%======================
% additional statistics

addedStats.distance = distance;

%don't forget to turn the warnings back on
warning(warningState);