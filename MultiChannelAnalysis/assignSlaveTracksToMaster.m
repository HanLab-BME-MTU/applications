function [map masterIdx slaveIdx unassignedMasterIdx unassignedSlaveIdx] = assignSlaveTracksToMaster(dataSlave, dataMaster, varargin)

ip = inputParser;
ip.CaseSensitive = false;
% ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('FileName', 'ProcessedTracks.mat', @ischar);
ip.addParamValue('MaxDistance', 5, @isscalar);
ip.addParamValue('MinOverlap', 1, @isscalar);
ip.addParamValue('MaxSlaveDelta_f', 10, @isscalar);
ip.parse(varargin{:});
minOverlap = ip.Results.MinOverlap;

nd = numel(dataMaster);
map = cell(1,nd);
masterIdx = cell(1,nd);
slaveIdx = cell(1,nd);
unassignedMasterIdx = cell(1,nd);
unassignedSlaveIdx = cell(1,nd);

fprintf('Slave/master assignment: processing track set ');
for i = 1:nd
    fprintf('%d/%d', i, nd);
    % load tracks
    %(category Ia only, to match the index in lftData loaded by runLifetimeAnalysis())
    masterTracks = loadTracks(dataMaster(i), 'Mask', true, 'Category', 'Ia', 'Cutoff_f', 0, 'FileName', ip.Results.FileName);
    slaveTracks = loadTracks(dataSlave(i), 'Mask', true, 'Category', 'Ia', 'Cutoff_f', 0, 'FileName', ip.Results.FileName);
    
    % mean positions
    mMeans = arrayfun(@(t) [nanmean(t.x(1,:)) nanmean(t.y(1,:))], masterTracks, 'UniformOutput', false);
    mMeans = vertcat(mMeans{:});
    sMeans = arrayfun(@(t) [nanmean(t.x(1,:)) nanmean(t.y(1,:))], slaveTracks, 'UniformOutput', false);
    sMeans = vertcat(sMeans{:});
    
    % slave tracks in spatial vicinity to master tracks (ignoring time overlap for now)
    % look for slave tracks close to master
    candIdx = KDTreeBallQuery(sMeans, mMeans, ip.Results.MaxDistance);
    
    % masterLengths = arrayfun(@(i) numel(i.t), masterTracks);
    slaveLengths = arrayfun(@(i) numel(i.t), slaveTracks);
    
    nM = numel(masterTracks);
    for k = 1:nM
        if ~isempty(candIdx{k})
            % overlaps
              overlapStarts = max([slaveTracks(candIdx{k}).start], masterTracks(k).start);
              overlapEnds = min([slaveTracks(candIdx{k}).end], masterTracks(k).end);
              overlaps = overlapEnds - overlapStarts + 1;
              % retain slave tracks with a minimal overlap value and maximal 
              sel = find(overlaps>=minOverlap & slaveLengths(candIdx{k})-overlaps <= ip.Results.MaxSlaveDelta_f);
              candIdx{k} = candIdx{k}(sel);
              if ~isempty(candIdx{k})
                  dist = zeros(1,numel(sel));
                  for s = 1:numel(sel)
                      mFrames = unique(masterTracks(k).f);
                      mFrames(isnan(mFrames)) = [];
                      imaster = mFrames>=overlapStarts(sel(s)) & mFrames<=overlapEnds(sel(s));
                      
                      sFrames = unique(slaveTracks(candIdx{k}(s)).f);
                      sFrames(isnan(sFrames)) = [];
                      islave = sFrames>=overlapStarts(sel(s)) & sFrames<=overlapEnds(sel(s));
                      
                       dist(s) = min(sqrt( (masterTracks(k).x(1,imaster)-slaveTracks(candIdx{k}(s)).x(1,islave)).^2 +...
                           (masterTracks(k).y(1,imaster)-slaveTracks(candIdx{k}(s)).y(1,islave)).^2));
                  end
                  candIdx{k} = candIdx{k}(dist<=5);
              end
              %candIdx{k} = candIdx{k}((masterLengths(k) >= slaveLengths(candIdx{k})));
        end
    end
    map{i} = candIdx;
    masterIdx{i} = find(cellfun(@(x) ~isempty(x), candIdx));
    slaveIdx{i} = unique(vertcat(candIdx{:})); % temporary fix, query above needs to be the other way around
    unassignedMasterIdx{i} = setdiff(1:numel(masterTracks), masterIdx{i});
    unassignedSlaveIdx{i} = setdiff(1:numel(slaveTracks), slaveIdx{i});
    
    fprintf('\b\b\b');
end
fprintf('\n');
