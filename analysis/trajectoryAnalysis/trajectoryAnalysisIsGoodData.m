function [data,rmIdx,isRun] = trajectoryAnalysisIsGoodData(data,constants,standardTag1,standardTag2)
%TRAJECTORYANALYSISISGOODDATA tests if the data structure is built correctly
%
% The function will return an error if there is a problem. It will also add
% the field .info.tags if necessary
% 
%
% SYNOPSIS [data,rmIdx,isRun] = trajectoryAnalysisIsGoodData(data,constants,standardTag1,standardTag2)
% 
% INPUT    data structure "inputData" (see help trajectoryAnalysis for more info)
%              or "run" if there are several data sets
%          constants    : needs field constants.MINTRAJECTORYLENGTH
%          standardTagX : (opt) strings with standard tag names
%                               standards: standardTag1, standardTag2
%
% OUTPUT   data structure "data" (see help trajectoryAnalysis for more
%              info) or "run" 
%          rmIdx  : index of trajectories that have been removed because of
%                   too few data points (only interesting for non-run)
%          isRun  : is the input given as run?
%
% pulled out of trajectoryAnalysis 3/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=================
% Test input
%=================

if nargin == 0 | isempty(data)
    error('empty or missing data structure ''data''');
end

if nargin < 2 | isempty(constants)
    error('please specify constants.MINTRAJECTORYLENGTH');
end

if nargin < 3 | isempty(standardTag1) | ~isstr(standardTag1)
    % set default standard tag 1
    standardTag1 = 'spb1';
end
if nargin < 3 | isempty(standardTag2) | ~isstr(standardTag2)
    % set default standard tag 2
    standardTag2 = 'cen1';
end

%================================




%====================
% Test data (main)
%====================

% find out whether we work with run or data here
if isfield(data,'data')
    % we have a run
    isRun = 1;
    run = data;
else
    % it's data. turn it into a run for now, so that we have to write less
    % code below
    isRun = 0;
    run.data = data;
end

% check for all runs
for iRun = 1:length(run)
    
    % find the missing mandatory fieldnames
    dataFieldNames = fieldnames(run(iRun).data);
    
    
    requiredDataFieldNames = {'distance';...
            'time';...
            'timePoints'};
    % setdiff(a,b) returns what is in a but not in b
    missingDataNames = setdiff(requiredDataFieldNames, dataFieldNames);
    
    if ~isempty(missingDataNames)
        errorMsg = sprintf('data structure is missing field(s) %s %s %s',missingDataNames{1:end});
        error(errorMsg);
    end
    
    % loop through data. add the info.tags field if it's missing. Complain if
    % time goes in the wrong way. Complain if inf, nan, or imaginary numbers
    
    if ~isfield(run(iRun).data,'info') 
        run(iRun).data(1).info = [];
    end
    
    rmIdx = [];
    nData = length(run(iRun).data);
    for iData = 1:nData
        
        
        distSize = size(run(iRun).data(iData).distance);
        timeSize = size(run(iRun).data(iData).time);
        tpSize   = size(run(iRun).data(iData).timePoints);
        
        if ~isequal(distSize(1),timeSize(1),tpSize(1))
            errorMsg = sprintf('in data(%i): distance and time have to be tx2, and timePoints tx1 arrays, respectively!',iData);
            error(errorMsg)
        end
        
        % check real and nan and inf distances
        if any(~isfinite(run(iRun).data(iData).distance)) | ~isreal(run(iRun).data(iData).distance)
            errorMsg = sprintf('run(%i).data(%i).distance is either nan or inf or imaginary!',iRun,iData);
            error(errorMsg)
        end
        
        if any(diff(run(iRun).data(iData).time,1,1)<=0) | any(diff(run(iRun).data(iData).timePoints,1,1)<=0)
            errorMsg = sprintf('run(%i).data(%i).time or .timePoints is not strictly increasing!',iRun,iData);
            error(errorMsg)
        end
        
        % check field info.tags
        if ~isfield(run(iRun).data(iData).info,'tags') | isempty(run(iRun).data(iData).info.tags)
            run(iRun).data(iData).info.tags = {standardTag1,standardTag2};
        end
        
        if distSize(1) < constants.MINTRAJECTORYLENGTH
            warningMsg = sprintf('run(%i).data(%i) contains less than %i data points. removed from list',iRun,iData,constants.MINTRAJECTORYLENGTH);
            warning('TRAJECTORYANALYSIS:checkInput',warningMsg)
            rmIdx = [rmIdx;iData];
        end
    end
    
    % remove data that contains not enough points
    run(iRun).data(rmIdx) = [];
    
    if isempty(run(iRun).data)
        errorMsg = sprintf('no valid data points for run %i after removal',iRun);
        error(errorMsg)
    end
    
    if isRun
        % test fileNameList
        if ~isfield(run(iRun), 'fileNameList') | ~iscell(run(iRun).fileNameList) | length(run(iRun).fileNameList)~=nData-length(rmIdx)
            errorMsg = sprintf('run(%i) must have a field ''fileNameList'' of the same length as data!',iRun);
            error(errorMsg)
        else
            % remove the entries that have not enough fileNames
            run(iRun).fileNameList(rmIdx) = [];
        end
    end
    
end % for iRun = 1:length(iRun)