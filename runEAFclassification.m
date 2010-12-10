% runEAFclassification(data, varargin)
%
% INPUT     data        : array of experiment structures
%           {overwrite} : 1 to overwrite previous results
%
% Usage example: runEAFclassification(data, 'overwrite', 1);

% Francois Aguet, December 2010

function data = runEAFclassification(data, lifetimeCutoff, varargin)

% defaults:
if nargin<2
    lifetimeCutoff = 5;
end
overwrite = false;


nv = length(varargin);
if mod(nv,2)~=0
    error('Incorrect format for optional inputs.');
end
for k = 1:nv/2
    switch lower(varargin{2*k-1})
        case 'overwrite'
            overwrite = logical(varargin{2*k});
        otherwise
            error('Unrecognized option');
    end
end


nExp = length(data);
for i = 1:nExp
    data(i).tracks = [];
    data(i).EAFposPct = [];
end

parfor i = 1:nExp
    if exist([data(i).source filesep 'Tracking' filesep 'trackAnalysis.mat'],'file') == 2
        tdata = load([data(i).source 'Tracking' filesep 'trackAnalysis.mat']);
        tracks = tdata.tracks;
        % select only valid tracks
        tracks = tracks([tracks.lifetime]>=lifetimeCutoff & [tracks.valid]==1);
        % run classification
        tracks = classifySlaveChannels(data(i), tracks);
        data(i).tracks = tracks;

        nPos = sum(arrayfun(@(x) x.cStatus(2), tracks));
        data(i).EAFposPct = 100*nPos/length(tracks);
    else
        fprintf('runEAFclassification: no ''trackAnalysis.mat'' found for movie %d.\n', i);
    end
end