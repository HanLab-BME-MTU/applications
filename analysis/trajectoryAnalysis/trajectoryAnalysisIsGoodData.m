function [data,rmIdx] = trajectoryAnalysisIsGoodData(data,constants,standardTag1,standardTag2)
%TRAJECTORYANALYSISISGOODDATA tests if the data structure is built correctly
%
% The function will return an error if there is a problem. It will also add
% the field .info.tags if necessary
% 
%
% SYNOPSIS data = trajectoryAnalysisIsGoodData(data)
% 
% INPUT    data structure "data" (see help trajectoryAnalysis for more info)
%          constants    : needs field constants.MINTRAJECTORYLENGTH
%          standardTagX : (opt) strings with standard tag names
%                               standards: standardTag1, standardTag2
%
% OUTPUT   data structure "data" (see help trajectoryAnalysis for more info)
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

% find the missing mandatory fieldnames
dataFieldNames = fieldnames(data);

if any(strcmp(dataFieldNames,'data'))
    error('run structure not supported as input to trajectoryAnalysis yet')
end

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

if ~isfield(data,'info') 
    data(1).info = [];
end

rmIdx = [];
for iData = 1:length(data)
    
    
    distSize = size(data(iData).distance);
    timeSize = size(data(iData).time);
    tpSize   = size(data(iData).timePoints);
    
    if ~isequal(distSize(1),timeSize(1),tpSize(1))
        errorMsg = sprintf('in data(%i): distance and time have to be tx2, and timePoints tx1 arrays, respectively!',iData);
        error(errorMsg)
    end
    
    % check real and nan and inf distances
    if any(~isfinite(data(iData).distance)) | ~isreal(data(iData).distance)
        errorMsg = sprintf('data(%i).distance is either nan or inf or imaginary!',iData);
        error(errorMsg)
    end
    
    if any(diff(data(iData).time)<=0) | any(diff(data(iData).timePoints)<=0)
        errorMsg = sprintf('data(%i).time or .timePoints is not strictly increasing!',iData);
        error(errorMsg)
    end
    
    % check field info.tags
    if ~isfield(data(iData).info,'tags') | isempty(data(iData).info.tags)
        data(iData).info.tags = {standardTag1,standardTag2};
    end
    
    if distSize(1) < constants.MINTRAJECTORYLENGTH
        warningMsg = sprintf('data(%i) contains less than %i data points. removed from list',iData,constants.MINTRAJECTORYLENGTH);
        warning('TRAJECTORYANALYSIS:checkInput',warningMsg)
        rmIdx = [rmIdx;iData];
    end
end

% remove data that contains not enough points
data(rmIdx) = [];

if isempty(data)
    error('no valid data points after removal')
end