% runTrackAnalysis
%
% INPUT     data        : array of experiment structures
%           {overwrite} : '1' to overwrite previous results

% Francois Aguet, November 2010

function runTrackAnalysis(data, overwrite, buffer)

if nargin<2 || isempty(overwrite)
    overwrite = 0;
end
if nargin < 3
    buffer = [];
end


nExp = length(data);
for i = 1:nExp
    
    if exist([data(i).source filesep 'Tracking' filesep 'trackAnalysis.mat'],'file') ~= 2 || overwrite
        trackAnalysis(data(i), buffer, []);
    else 
        display(['skipping movie ' num2str(i)])
    end
end