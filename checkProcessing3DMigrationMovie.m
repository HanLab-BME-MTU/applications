function [procStatus,statInfo] = checkProcessing3DMigrationMovie(movieList)


%TEMP - still using my stupid hard-coded process list and index
processNames = {'SegmentationProcess3D',...             % 1
                'MaskGeometry3DProcess',...             % 2
                'SkeletonizationProcess',...            % 3
                'SkeletonPruningProcess',...            % 4
                'MaskObjectTrackingProcess'};           % 5            
procInd = 1;            

nProc = numel(processNames);            
            
nMov = numel(movieList.movies_);

procStatus = ones(nMov,nProc)* -1;%Negative one indicates no process present
statInfo(1:nMov,1:nProc) = struct('CompleteTime',NaN,'CompleteTimeString','','CompleteTimeNum',0);

for j = 1:nMov
    
    for k = 1:nProc
        
        currProcInd = movieList.movies_{j}.getProcessIndex(processNames{k},1,0);
        
        if ~isempty(currProcInd)
            
            procStatus(j,k) = movieList.movies_{j}.processes_{currProcInd}.checkChannelOutput(procInd);
            if procStatus(j,k)
                statInfo(j,k).CompleteTime = movieList.movies_{j}.processes_{currProcInd}.finishTime_;
                statInfo(j,k).CompleteTimeString = datestr(movieList.movies_{j}.processes_{currProcInd}.finishTime_);
                statInfo(j,k).CompleteTimeNum = datenum(movieList.movies_{j}.processes_{currProcInd}.finishTime_);
            end
        end
    end
    
end