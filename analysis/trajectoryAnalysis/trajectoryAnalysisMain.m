function trajectoryDescription = trajectoryAnalysisMain(data,constants,showDetails,doConvergence,verbose,fileNameList);
%TRAJECTORYANALYSISMAIN calls the individual programs to analyze the trajectories
%
% The trajectory is first scanned for significant inidividual changes of
% the distance between the two tags. This creates dataListS. Then, the
% software looks for pauses and eventually for longer-than-1-interval
% growth and shrinkabe/separation and congression events, by trying to get
% the longest possible fits. This yields dataListG. Finally, the statistics
% are calculated.
%
% for the output structure, type help trajectoryAnalysisMainCalcStats.
%
% dataList has the following cols
%     1:startIdx, 2:endIdx, 3:state, 4:slope, 5:slopeSigma, 6:slopeSigmaAPR,
%     7:deltaT, 8:(deltaTSigma), 9:deltaD, 10:deltaDSigma, 11:startDistance
%
% c: 1/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%init dataListAll
tpAll = cat(1,data.timePoints);
lengthData = length(data);
dataAll = zeros(length(tpAll)+lengthData,11);
dataAllFirstIdx = 1;
distanceAll = cat(1,data.distance);
timeAll = cat(1,data.time);
tpInd = 0;
dlaLength = 0;
individualMeanDistances = [];

%for plotting: calc extreme coords
yAxisLimits = [0,max(sum(distanceAll(:,[1,2,2]),2))];
xAxisLimits = [0,max(sum(timeAll,2))];

% init trajectoryDescription
trajectoryDescription = struct('overallStatistics',[],'individualStatistics',[],'convergenceStatistics',[],...
    'convergenceClusters',[],'overallClusters',[]);

%--------LOOP THROUGH EVERYTHING
for numData = 1:length(data)
    
    %get data
    distance   = data(numData).distance;
    time       = data(numData).time;
    timePoints = data(numData).timePoints;
    
    %number of timepoints
    nTimePoints = length(timePoints);    
    
    %initialize storage matrix: 1:startIdx, 2:endIdx, 3:state, 4:slope, 5:slopeSigma, 6:slopeSigmaAPR,
    %7:deltaT, 8:(deltaTSigma), 9:deltaD, 10:deltaDSigma, 11:startDistance
    dataList = zeros(nTimePoints-1,11);
    
    %-------- FILL IN DATALIST--------
    
    %fill in dataIdx, time, distances
    dataList(:,[1,2])    = [[1:nTimePoints-1]',[2:nTimePoints]'];
    dataList(:,[7])      = [diff(time(:,1),1)];
    dataList(:,[9,11])   = [diff(distance(:,1),1),distance(1:end-1,1)];
    
    %calculate deltaDistanceSigma
    dataList(:,10) = sqrt(distance(1:end-1,2).^2+distance(2:end,2).^2);
    
    %-----END FILL IN DATALIST--------
    
    
    
    %-------- SEED STATES
    %this module does an initial fit of the smallest allowed unit and
    %classifies them in the following way
    %
    %growth = 1
    %shrinkage = 2
    %pause = 3
    %undetermined = 0
    %frame deleted = -1, -2 etc, depending on how many frames have been
    %deleted in between two timepoints
    %
    %currently, only a initial classification of doublets is implemented;
    %for triplet-fitting, we would proceed as below with fitting all
    %allowed and then resolving conflicts
    
    dataListS = trajectoryAnalysisMainSeedStates(dataList,distance,time,timePoints,constants);
    
    %-----END SEED STATES
    
    
    %---------PLOT
    if any(verbose == 2)
        plotFigH = figure('Name',fileNameList{numData});
        axes('NextPlot','add','YLim',yAxisLimits,'XLim',xAxisLimits);
        plot(time(:,1),distance(:,1),'-dk');
        myErrorbar(time(:,1),distance(:,1),distance(:,2));
        
        
        %---plot deleted intervals with '--w' here
   
        %---growth
        %read indices
        growthIdx = find(dataListS(:,3)==1);
        growthStartDataIdx = (dataListS(growthIdx,1));
        growthEndDataIdx = (dataListS(growthIdx,2));
        
        %build NaN-separated lists for time and distance
        %time
        growthTime = repmat(NaN,[3,length(growthIdx)]);
        growthTime(1,:) = time(growthStartDataIdx,1)';
        growthTime(2,:) = time(growthEndDataIdx,1)';
        growthTime = growthTime(:);
        %distance
        growthDistance = repmat(NaN,[3,length(growthIdx)]);
        growthDistance(1,:) = distance(growthStartDataIdx,1)';
        growthDistance(2,:) = distance(growthEndDataIdx,1)';
        growthDistance = growthDistance(:);
        
        %plot data
        plot(growthTime,growthDistance,'--g');
        
        %----shrinkage
        %read indices
        shrinkageIdx = find(dataListS(:,3)==2);
        shrinkageStartDataIdx = (dataListS(shrinkageIdx,1));
        shrinkageEndDataIdx = (dataListS(shrinkageIdx,2));
        
        %build NaN-separated lists for time and distance
        %time
        shrinkageTime = repmat(NaN,[3,length(shrinkageIdx)]);
        shrinkageTime(1,:) = time(shrinkageStartDataIdx,1)';
        shrinkageTime(2,:) = time(shrinkageEndDataIdx,1)';
        shrinkageTime = shrinkageTime(:);
        %distance
        shrinkageDistance = repmat(NaN,[3,length(shrinkageIdx)]);
        shrinkageDistance(1,:) = distance(shrinkageStartDataIdx,1)';
        shrinkageDistance(2,:) = distance(shrinkageEndDataIdx,1)';
        shrinkageDistance = shrinkageDistance(:);
        
        %plot data
        plot(shrinkageTime,shrinkageDistance,'--r');
        
        %---deleted
        %read indices
        deletedIdx = find(dataListS(:,3)<0);
        deletedStartDataIdx = (dataListS(deletedIdx,1));
        deletedEndDataIdx = (dataListS(deletedIdx,2));
        
        %build NaN-separated lists for time and distance
        %time
        deletedTime = repmat(NaN,[3,length(deletedIdx)]);
        deletedTime(1,:) = time(deletedStartDataIdx,1)';
        deletedTime(2,:) = time(deletedEndDataIdx,1)';
        deletedTime = deletedTime(:);
        %distance
        deletedDistance = repmat(NaN,[3,length(deletedIdx)]);
        deletedDistance(1,:) = distance(deletedStartDataIdx,1)';
        deletedDistance(2,:) = distance(deletedEndDataIdx,1)';
        deletedDistance = deletedDistance(:);
        
        %plot data
        plot(deletedTime,deletedDistance,'-b');
    end
    %-----END PLOT
    
    
    %--------FIND PAUSES
    %now that we have found the initial growth or shrinkage (or pauses), we
    %go and look for pauses specifically. We could do the fitting jointly,
    %but then we might include pauses in larger growth phases if we use
    %the topDown-strategy
    
    compatible = [0,3,-constants.MAXDELETED]; %define compatible states
    dataListG = trajectoryAnalysisMainGroupUnits(dataListS,distance,time,constants,compatible);
    
    %-----END FIND PAUSES
    
    %---------PLOT
    if any(verbose == 2)
        figure(plotFigH)
        %----pause
        %read indices
        pauseIdx = find(dataListG(:,3)==3);
        pauseStartDataIdx = (dataListG(pauseIdx,1));
        pauseEndDataIdx = (dataListG(pauseIdx,2));
        
        %build NaN-separated lists for time and distance
        %time
        pauseTime = repmat(NaN,[3,length(pauseIdx)]);
        pauseTime(1,:) = time(pauseStartDataIdx,1)';
        pauseTime(2,:) = time(pauseEndDataIdx,1)';
        pauseTime = pauseTime(:);
        
        %distance
        pauseDistance = repmat(NaN,[3,length(pauseIdx)]);
        pauseDistance(1,:) = dataListG(pauseIdx,11)';
        pauseDistance(2,:) = dataListG(pauseIdx,11)';
        pauseDistance = pauseDistance(:);
        
        
        %plot data
        plot(pauseTime,pauseDistance,'-c');
    end
    %-----END PLOT
    
    
    %-------- GROW STATES
    %proceed to increase the length of the individual units until we have
    %completely described the trajectory. Thanks to the pauses, this
    %process is somewhat sped up, because we do not have as many conflicts
    
    compatible = [0,1,-constants.MAXDELETED;0,2,-constants.MAXDELETED]; %define compatible states
    dataListG = trajectoryAnalysisMainGroupUnits(dataListG,distance,time,constants,compatible);
    
    %-----END GROW STATES
    
    %---------PLOT
    if any(verbose == 2)
        %---growth
        %read indices
        growthIdx = find(dataListG(:,3)==1);
        growthStartDataIdx = (dataListG(growthIdx,1));
        growthEndDataIdx = (dataListG(growthIdx,2));
        
        %build NaN-separated lists for time and distance
        %time
        growthTime = repmat(NaN,[3,length(growthIdx)]);
        growthTime(1,:) = time(growthStartDataIdx,1)';
        growthTime(2,:) = time(growthEndDataIdx,1)';
        growthTime = growthTime(:);
        %distance
        
        %use start/end distance of fit
        growthDistance = repmat(NaN,[3,length(growthIdx)]);
        growthDistance(1,:) = dataListG(growthIdx,11)'; %startDistance
        growthDistance(2,:) = (dataListG(growthIdx,11)+dataListG(growthIdx,9))'; %startDist.+deltaD
        growthDistance = growthDistance(:);
        
        %plot data
        plot(growthTime,growthDistance,'-g');
        
        %----shrinkage
        %read indices
        shrinkageIdx = find(dataListG(:,3)==2);
        shrinkageStartDataIdx = (dataListG(shrinkageIdx,1));
        shrinkageEndDataIdx = (dataListG(shrinkageIdx,2));
        
        %build NaN-separated lists for time and distance
        %time
        shrinkageTime = repmat(NaN,[3,length(shrinkageIdx)]);
        shrinkageTime(1,:) = time(shrinkageStartDataIdx,1)';
        shrinkageTime(2,:) = time(shrinkageEndDataIdx,1)';
        shrinkageTime = shrinkageTime(:);
        %distance
        
        %use start/end distance of fit
        shrinkageDistance = repmat(NaN,[3,length(shrinkageIdx)]);
        shrinkageDistance(1,:) = dataListG(shrinkageIdx,11)'; %startDistance
        shrinkageDistance(2,:) = (dataListG(shrinkageIdx,11)+dataListG(shrinkageIdx,9))'; %startDist.+deltaD
        shrinkageDistance = shrinkageDistance(:);
        
        %plot data
        plot(shrinkageTime,shrinkageDistance,'-r');
    end
    %-----END PLOT
    
    
    %---------CALCULATE STATISTICS
    if constants.DOCLUSTER == 3 % check if we have to do individual clustering
        [trajectoryDescription.individualStatistics(numData).summary, ...
                distributions,...
                trajectoryDescription.individualStatistics(numData).speedClustering] =...
            trajectoryAnalysisMainCalcStats(dataListG,distance,time,verbose,fileNameList{numData},constants);
    else
        [trajectoryDescription.individualStatistics(numData).summary, distributions] =...
            trajectoryAnalysisMainCalcStats(dataListG,distance,time,verbose,fileNameList{numData},constants);
    end
    %store fileName
    trajectoryDescription.individualStatistics(numData).name = fileNameList{numData};
    
    if showDetails
        trajectoryDescription.individualStatistics(numData).dataListSeed = dataListS;
        trajectoryDescription.individualStatistics(numData).dataListGroup = dataListG;
        trajectoryDescription.individualStatistics(numData).distributions = distributions;
    end
    
    %store data for overall trajectory
    dlgLength = size(dataListG,1);
    dataListG(:,1:2) = dataListG(:,1:2) + tpInd; %update indices, so that they point to distanceAll, timeAll
    dataListAll(dlaLength+1:dlaLength+dlgLength,:) = dataListG; %add dataListG
    dataListAll(dlaLength + dlgLength +1, 3) = -99; %add separator between lists
    dlaLength = dlaLength + dlgLength + 1; %update length
    tpInd = tpInd + nTimePoints; %update individual length
    
    % store data for convergence/overall statistics
    individualMeanDistances = [individualMeanDistances;...
            trajectoryDescription.individualStatistics(numData).summary.distanceMean(1),...
            trajectoryDescription.individualStatistics(numData).summary.distanceStd(1),...
            trajectoryDescription.individualStatistics(numData).summary.nTimepoints(1)];
    
    
    if doConvergence
        if constants.DOCLUSTER > 1 %cluster for convergence only, and if individuals are selected
            [trajectoryDescription.convergenceStatistics(numData),dummy,...
                    trajectoryDescription.convergenceClusters(numData)] = ...
                trajectoryAnalysisMainCalcStats(dataListAll(1:dlaLength,:),distanceAll(1:tpInd,:),timeAll(1:tpInd,:),verbose,...
                ['convergence statistics for ',num2str(numData),' trajectories'],constants,individualMeanDistances);
        else
            trajectoryDescription.convergenceStatistics(numData) = ...
                trajectoryAnalysisMainCalcStats(dataListAll(1:dlaLength,:),distanceAll(1:tpInd,:),timeAll(1:tpInd,:),verbose,...
                ['convergence statistics for ',num2str(numData),' trajectories'],constants,individualMeanDistances);
        end
    end
    %-----END CALCULATE STATISTICS
    
end %for numData = 1:length(data)

%---------CALCULATE OVERALL STATISTICS
%if we did convergence, this will just be the nth entry, else we have a
%field overallStatistics
if constants.DOCLUSTER
[overallStatistics,overallDistribution,overallCluster] = ...
        trajectoryAnalysisMainCalcStats(dataListAll(1:dlaLength,:),distanceAll(1:tpInd,:),timeAll(1:tpInd,:),verbose,...
        ['overall statistics for ',num2str(numData),' trajectories'],constants,individualMeanDistances);
else
    [overallStatistics,overallDistribution] = ...
        trajectoryAnalysisMainCalcStats(dataListAll(1:dlaLength,:),distanceAll(1:tpInd,:),timeAll(1:tpInd,:),verbose,...
        ['overall statistics for ',num2str(numData),' trajectories'],constants,individualMeanDistances);
    overallCluster = [];
end

if doConvergence
    %create a field but leave it empty
    trajectoryDescription.overallStatistics = [];
else
    %calculate overallStatistics
    trajectoryDescription.overallStatistics = overallStatistics;
    %assign empty for convergence
    trajectoryDescription.convergenceStatistics = [];
end

trajectoryDescription.overallDistribution = overallDistribution;
trajectoryDescription.overallClusters     = overallCluster;