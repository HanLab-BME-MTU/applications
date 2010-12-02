% runTrackAnalysis
%
% INPUT     data        : array of experiment structures
%           {buffer}    : length of buffer before and after tracks
%           {overwrite} : 1 to overwrite previous results
%
% Usage example: runTrackAnalysis(data, 'buffer', 10);

% Francois Aguet, November 2010

function runTrackAnalysis(data, varargin)

% defaults:
buffer = [];
overwrite = 0;

nv = length(varargin);
if mod(nv,2)~=0
    error('Incorrect format for optional inputs.');
end
for k = 1:nv/2
    switch lower(varargin{2*k-1})
        case 'buffer'
            buffer = varargin{2*k};
        case 'overwrite'
            overwrite = 1;
        otherwise
            error('Unrecognized option');
    end
end

nExp = length(data);
parfor i = 1:nExp
    if exist([data(i).source filesep 'Tracking' filesep 'trackAnalysis.mat'],'file') ~= 2 || overwrite
        trackAnalysis(data(i), buffer, []);
    else
        fprintf('TrackAnalysis: movie %d has already been analyzed.\n', i);
    end
end