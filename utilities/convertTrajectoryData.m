function convertedData = convertTrajectoryData(data,join)
%CONVERTTRAJECTORYDATA converts trajectory data into Khulouds Nan-separated lists and vice versa
%
% SYNOPSIS convertedData = convertTrajectoryData(data)
%
% INPUT    data: either trajectoryAnalysis data structure with fields distance (ntpx2), time(ntpx2),
%                timepoints (ntpx1) or simulation data structure with fields
%                distance (timeIntx2), timeInterval (1x1)
%          join: (opt) whether to join the individual elements of the struct
%                array data. Default is 0 for sim2TA and 1 for TA2sim
%          
%
% OUTPUT   convertedData : the other data type
%
%c: jonas, 04/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==============
% TEST INPUT
%==============

if isempty(data) | ~isstruct(data)
    error('input for CONVERTTRAJECTORYDATA has to be a structure')
end

% test for taData
if isfield(data,'distance') & isfield(data,'time') & isfield(data,'timePoints')
    dataType = 1;
elseif isfield(data,'distance') & isfield(data,'timeInterval')
    dataType = 2;
else
    error('input structure for CONVERTTRAJECTORYDATA does not contain the necessary fields. Please seek help.')
end

% test join. Default = 0 for dataType==1 and 1 for dataType==2
if nargin < 2 | isempty(join)
    join = dataType-1;
end

%----------------


%================
% CONVERT
%================

switch dataType*2 + join
    %--------TrajectoryAnalysisData to Nan-separated list---------
    case 1 % TA2sim, no join
        
        % get the length of data
        numData = length(data);
        
        % find the number of timepoints and prepare output data
        nanCell = cell(numData,1);
        for iData = 1:numData
            % prepare Nan-arrays for distance
            nanCell{iData} = repmat(NaN,[data(iData).timePoints(end),2]);
        end
        
        % init array. Thanks to the nanCell, it already has the right size
        convertedData = struct('distance',nanCell,'timeInterval',0);
        
        % now loop through data and fill in distance and calculate average
        % time interval
        for iData = 1:numData
            
            % read distance into output
            convertedData(iData).distance(data(iData).timePoints,:) = data(iData).distance;
            
            % find mean time by finding where we had just one timepoint
            % difference
            timeAndTp = [data(iData).time,data(iData).timePoints];
            deltaTaT = diff(timeAndTp,1,1);
            
            convertedData(iData).timeInterval = mean(deltaTaT(find(deltaTaT(:,3)==1),1));
            
        end %for iData = 1:numData
        
        %END CASE 1
        
    case 2 %TA2sim, join
        
        % get the length of data
        numData = length(data);
        
        
        % since we will not need to loop to collect the distances, cat
        % already now
        allDist = cat(1,data.distance);
        
        %--- find the number of timepoints and prepare output data
        nanNum = -1;
        tpList = zeros(size(allDist,1),1);
        tpListIdx = 1;
        
        for iData = 1:numData
            % the individual trajectories have to be nan-separated
            nanNum = nanNum + 1;
            
            % fill in timepointList
            tpListIdxNew = tpListIdx + length(data(iData).timePoints)-1;
            tpList(tpListIdx:tpListIdxNew) = data(iData).timePoints + nanNum;
            tpListIdx = tpListIdxNew + 1;
            
            % count timepoints to build nan-array
            nanNum = nanNum + data(iData).timePoints(end);
            
        end
        
        % init array. Thanks to the nanCell, it already has the right size
        convertedData = struct('distance',repmat(NaN,[nanNum,2]),'timeInterval',0);
        
        
        % read distance into output
        convertedData.distance(tpList,:) = allDist;
        
        % find mean time by finding where we had just one timepoint
        % difference
        timeAndTp = [cat(1,data.time),cat(1,data.timePoints)];
        deltaTaT = diff(timeAndTp,1,1);
        
        convertedData.timeInterval = mean(deltaTaT(find(deltaTaT(:,3)==1),1));
        
        
        %END CASE 2
        
        %--------------- sim2TA -------------------    
    case 4 % sim2TA, no join
        
        % get the length of data
        numData = length(data);
        
        % init output data
        convertedData(1:numData) = struct('distance',[],'time',[],'timePoints',[]);
        
        % loop through data and assign output
        for iData = 1:numData
            
            % timepoints
            convertedData(iData).timePoints = find(isfinite(data(iData).distance(:,1)));
            
            % distance
            convertedData(iData).distance = data(iData).distance(convertedData(iData).timePoints,:);
            
            % time. Starts at min 1 time interval
            convertedData(iData).time = zeros(size(convertedData(iData).distance));
            convertedData(iData).time(:,1) = convertedData(iData).timePoints * data(iData).timeInterval;
            
        end
        
        %END CASE 3
        
    case 5 % sim2TA, join
        
        % get the length of data
        numData = length(data);
        
        % test sampling. if not all equal, exit. We do not allow uneven sampling in one and the same trajectory
        allIntervals = cat(1,data.timeInterval);
        if any(allIntervals-allIntervals(1))
            error('you are not allowed to combine trajectories with different sampling intervals')
        end
        
        % init output data
        convertedData(1:numData) = struct('distance',[],'time',[],'timePoints',[]);
        
        allDist = cat(1,data.distance);
        
        % timepoints
        convertedData.timePoints = find(isfinite(allDist(:,1)));
        
        % distance
        convertedData.distance = allDist(convertedData.timePoints,:);
        
        % time. Starts at min 1 time interval
        convertedData.time = zeros(size(convertedData.distance));
        convertedData.time(:,1) = convertedData.timePoints * data(1).timeInterval;
        
        %END CASE 4
        
    otherwise
        error('unhandled case')
        
end % switch