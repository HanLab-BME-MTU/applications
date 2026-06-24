function adhesionTracks = tracksFinal2adhesionTracks(tracksFinal, pixelSize, timeInterval)
%TRACKSFINAL2ADHESIONTRACKS  Convert u-track output (tracksFinal compound
%track struct from trackCloseGapsKalmanSparse) into the FilopodiaForcePackage
%adhesionTracks struct used by Processes 4-6.
%
%   adhesionTracks = tracksFinal2adhesionTracks(tracksFinal, pixelSize, timeInterval)
%
% Each output element has:
%   trackId, frames (1xN), pos (Nx2 [x y]), amp (1xN),
%   lifetime (#frames), lifetime_s (s)
% Note: 'dist'/'dist_nm' (signed distance to the cell-body edge) are NOT set
% here; Process 4 recomputes those against the body mask. Gap frames (NaN
% coordinates after gap closing) are dropped so pos contains only real
% detections, matching the old nearest-neighbour tracker's output.
% Sangyoon J. Han / 2026

if nargin < 2 || isempty(pixelSize),    pixelSize = 1;    end
if nargin < 3 || isempty(timeInterval), timeInterval = 1; end

adhesionTracks = struct('trackId',{},'frames',{},'featIdx',{},'pos',{}, ...
    'amp',{},'dist',{},'dist_nm',{},'lifetime',{},'lifetime_s',{});

if isempty(tracksFinal), return; end

% Flatten compound tracks into a matrix (ignore merge/split branches) plus a
% per-row start frame, using u-track's own utility when available.
useConv = exist('convStruct2MatIgnoreMS','file')==2;

kept = 0;
for i = 1:numel(tracksFinal)
    tr = tracksFinal(i);
    cg = tr.tracksCoordAmpCG;            % [nSeg x 8*nFrameSpan]
    if isempty(cg), continue; end

    % start frame of this compound track
    soe = tr.seqOfEvents;
    startFrame = min(soe(:,1));

    for s = 1:size(cg,1)                 % each segment -> one adhesionTrack
        row = cg(s,:);
        nSpan = numel(row)/8;
        xs = row(1:8:end);
        ys = row(2:8:end);
        as = row(4:8:end);
        fr = startFrame + (0:nSpan-1);

        valid = ~isnan(xs) & ~isnan(ys); % drop gap-closed (interpolated) frames
        if nnz(valid) < 1, continue; end
        xs = xs(valid); ys = ys(valid); as = as(valid); fr = fr(valid);

        kept = kept + 1;
        adhesionTracks(kept).trackId    = kept;
        adhesionTracks(kept).frames     = fr(:)';
        adhesionTracks(kept).featIdx    = nan(1,numel(fr));   % not tracked back to detection index
        adhesionTracks(kept).pos        = [xs(:), ys(:)];
        adhesionTracks(kept).amp        = as(:)';
        adhesionTracks(kept).dist       = nan(1,numel(fr));   % Process 4 fills these
        adhesionTracks(kept).dist_nm    = nan(1,numel(fr));
        adhesionTracks(kept).lifetime   = numel(fr);
        adhesionTracks(kept).lifetime_s = (fr(end)-fr(1)) * timeInterval;
    end
end
end
