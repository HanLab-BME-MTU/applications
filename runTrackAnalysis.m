% runTrackAnalysis
%
% INPUT     data        : array of experiment structures
%           {overwrite} : '1' to overwrite previous results

% Francois Aguet, November 2010

function runTrackAnalysis(data, overwrite)

if nargin<2
    overwrite = 0;
end

nExp = length(data);
parfor i = 1:nExp   
    trackAnalysis(data(i), overwrite);
end