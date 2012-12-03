%assignSlaveTracksToMaster(dataSlave, dataMaster, varargin) assigns a set of 'slave' tracks to a set of 'master' tracks
%
% Inputs:
%     dataSlave : slave data sets, or slave tracks
%    dataMaster : master data sets, or master tracks
%
% Options:
%          InputMode : Specifies whether the inputs are data ('data') or track ('tracks', default) sets.
%        MaxDistance : Maximum distance between slave and master track points, in pixels (default: 3).
%   MaxSlaveOverhang : Maximum length of the slave track before or after the master track, in frames (default: 5).
%         MinOverlap : Minimum admissible length of the overlap between the slave and master tracks, in frames (default: 1).

% Francois Aguet (last modified 09/2012)

function [map, masterIdx, slaveIdx, unassignedMasterIdx, unassignedSlaveIdx] = assignSlaveTracksToMaster(dataSlave, dataMaster, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('FileName', 'ProcessedTracks.mat', @ischar);
ip.addParamValue('MaxDistance', 3, @isscalar);
ip.addParamValue('MinOverlap', 1, @isscalar);
ip.addParamValue('MaxSlaveOverhang', 5, @isscalar);
ip.addParamValue('InputMode', 'tracks', @(x) any(strcmpi(x, {'tracks', 'data'})));
ip.parse(varargin{:});
minOverlap = ip.Results.MinOverlap;
mode = ip.Results.InputMode;

if strcmpi(mode, 'data')
    nd = numel(dataMaster);
else
    nd = 1;
end
map = cell(1,nd);
masterIdx = cell(1,nd);
slaveIdx = cell(1,nd);
unassignedMasterIdx = cell(1,nd);
unassignedSlaveIdx = cell(1,nd);

fprintf('Slave->master assignment: processing track set    ');
for i = 1:nd
    fprintf('\b\b\b%d/%d', i, nd);
    % load tracks
    if strcmpi(mode, 'data')
        % category Ia only, to match the index in lftData loaded by runLifetimeAnalysis()
        masterTracks = loadTracks(dataMaster(i), 'Mask', true, 'Category', 'Ia', 'Cutoff_f', 0, 'FileName', ip.Results.FileName);
        slaveTracks = loadTracks(dataSlave(i), 'Mask', true, 'Category', 'Ia', 'Cutoff_f', 0, 'FileName', ip.Results.FileName);
    else
        masterTracks = dataMaster;
        slaveTracks = dataSlave;
    end
    % mean positions
    mMeans = arrayfun(@(t) [nanmean(t.x(1,:)) nanmean(t.y(1,:))], masterTracks, 'UniformOutput', false);
    mMeans = vertcat(mMeans{:});
    sMeans = arrayfun(@(t) [nanmean(t.x(1,:)) nanmean(t.y(1,:))], slaveTracks, 'UniformOutput', false);
    sMeans = vertcat(sMeans{:});
    
    % slave tracks in spatial vicinity to master tracks (ignoring time overlap for now)
    % look for slave tracks close to master
    candIdx = KDTreeBallQuery(sMeans, mMeans, ip.Results.MaxDistance);
    
    %masterLengths = arrayfun(@(i) numel(i.t), masterTracks);
    slaveLengths = arrayfun(@(i) numel(i.t), slaveTracks);
    
    nM = numel(masterTracks);
    for k = 1:nM
        if ~isempty(candIdx{k})
            % overlaps
            overlapStarts = max([slaveTracks(candIdx{k}).start], masterTracks(k).start);
            overlapEnds = min([slaveTracks(candIdx{k}).end], masterTracks(k).end);
            overlaps = overlapEnds - overlapStarts + 1;
            % retain slave tracks with a minimal overlap value and maximal
            sel = find(overlaps>=minOverlap & slaveLengths(candIdx{k})-overlaps <= ip.Results.MaxSlaveOverhang);
            candIdx{k} = candIdx{k}(sel);
            if ~isempty(candIdx{k})
                dist = zeros(1,numel(sel));
                for s = 1:numel(sel)
                    mFrames = unique(masterTracks(k).f);
                    mFrames(isnan(mFrames)) = [];
                    % binary mask of master frames overlapping with current slave segment
                    imaster = mFrames>=overlapStarts(sel(s)) & mFrames<=overlapEnds(sel(s));
                    
                    sFrames = unique(slaveTracks(candIdx{k}(s)).f);
                    sFrames(isnan(sFrames)) = [];
                    % binary mask of current slave segment overlapping with master
                    islave = sFrames>=overlapStarts(sel(s)) & sFrames<=overlapEnds(sel(s));
                    
                    % max distance between the two tracks
                    dist(s) = max(sqrt( (masterTracks(k).x(1,imaster)-slaveTracks(candIdx{k}(s)).x(1,islave)).^2 +...
                        (masterTracks(k).y(1,imaster)-slaveTracks(candIdx{k}(s)).y(1,islave)).^2));
                end
                % retain slave segments that are within 'MaxDistance' of master
                candIdx{k} = candIdx{k}(dist<=ip.Results.MaxDistance);
            end
            %candIdx{k} = candIdx{k}((masterLengths(k) >= slaveLengths(candIdx{k})));
        end
    end
    map{i} = candIdx;
    masterIdx{i} = find(cellfun(@(x) ~isempty(x), candIdx));
    slaveIdx{i} = unique(vertcat(candIdx{:}));
    unassignedMasterIdx{i} = setdiff(1:numel(masterTracks), masterIdx{i});
    unassignedSlaveIdx{i} = setdiff(1:numel(slaveTracks), slaveIdx{i});
end
fprintf('\n');
