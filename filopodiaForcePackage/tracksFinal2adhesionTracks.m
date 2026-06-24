function adhesionTracks = tracksFinal2adhesionTracks(tracksFinal, pixelSize, timeInterval, adhesionInfo)
%TRACKSFINAL2ADHESIONTRACKS  Convert u-track output (tracksFinal compound
%track struct from trackCloseGapsKalmanSparse) into the FilopodiaForcePackage
%adhesionTracks struct used by Processes 4-6.
%
%   adhesionTracks = tracksFinal2adhesionTracks(tracksFinal, pixelSize, ...
%                                               timeInterval, adhesionInfo)
%
% Each output element has:
%   trackId, frames (1xN), featIdx (1xN), pos (Nx2 [x y]), amp (1xN),
%   dist (1xN), dist_nm (1xN), lifetime, lifetime_s
%
% IMPORTANT: when adhesionInfo (the per-frame detection cell array saved by
% Process 2) is supplied, each tracked position is matched back to the
% nearest detection in that frame so that featIdx (index into
% adhesionInfo{t}) and dist (signed distance to the body edge) are recovered.
% Process 4 relies on featIdx to map adhesions to tracks and on dist to test
% tip reach; without them P4 cannot use the P2/P3 detections and falls back
% to re-estimating shafts near the body. So always pass adhesionInfo.
% Sangyoon J. Han / 2026

if nargin < 2 || isempty(pixelSize),    pixelSize = 1;    end
if nargin < 3 || isempty(timeInterval), timeInterval = 1; end
if nargin < 4, adhesionInfo = {}; end
haveAdh = ~isempty(adhesionInfo);

adhesionTracks = struct('trackId',{},'frames',{},'featIdx',{},'pos',{}, ...
    'amp',{},'dist',{},'dist_nm',{},'isInterp',{},'lifetime',{},'lifetime_s',{});
if isempty(tracksFinal), return; end

matchTol = 2;   % px; tracked pos vs detection must be within this to match

kept = 0;
for i = 1:numel(tracksFinal)
    tr = tracksFinal(i);
    cg = tr.tracksCoordAmpCG;
    if isempty(cg), continue; end
    startFrame = min(tr.seqOfEvents(:,1));

    for s = 1:size(cg,1)
        row = cg(s,:);
        nSpan = numel(row)/8;
        xs = row(1:8:end); ys = row(2:8:end); as = row(4:8:end);
        fr = startFrame + (0:nSpan-1);

        real = ~isnan(xs) & ~isnan(ys);   % frames with a real detection
        if nnz(real) < 1, continue; end

        % Keep the FULL span from first to last real detection and fill the
        % gap-closed frames in between by linear interpolation, so the track
        % stays continuous (this is the whole point of u-track gap closing).
        % Dropping NaN frames here is what made filopodia blink in/out.
        fi = find(real); a = fi(1); b = fi(end);
        keep = a:b;
        fr = fr(keep); xs = xs(keep); ys = ys(keep); as = as(keep);
        isInterp = isnan(xs) | isnan(ys);
        gi = find(~isInterp);
        if numel(gi) >= 2
            xs = interp1(fr(gi), xs(gi), fr, 'linear');
            ys = interp1(fr(gi), ys(gi), fr, 'linear');
            as = interp1(fr(gi), as(gi), fr, 'linear');
        end

        nF = numel(fr);
        featIdx = nan(1,nF); dist = nan(1,nF);
        if haveAdh
            for j = 1:nF
                if isInterp(j), continue; end       % no detection on gap frames
                t = fr(j);
                if t<1 || t>numel(adhesionInfo) || isempty(adhesionInfo{t}), continue; end
                A = adhesionInfo{t};
                ap = reshape([A.pos],2,[])';        % nDet x 2
                d2 = (ap(:,1)-xs(j)).^2 + (ap(:,2)-ys(j)).^2;
                [dm, kmin] = min(d2);
                if sqrt(dm) <= matchTol
                    featIdx(j) = kmin;
                    if isfield(A,'dist'), dist(j) = A(kmin).dist; end
                end
            end
        end
        % fill dist on interpolated frames so tip-reach tests stay smooth
        gd = find(~isnan(dist));
        if numel(gd) >= 2
            dist = interp1(fr(gd), dist(gd), fr, 'linear', 'extrap');
        end

        kept = kept + 1;
        adhesionTracks(kept).trackId    = kept;
        adhesionTracks(kept).frames     = fr(:)';
        adhesionTracks(kept).featIdx    = featIdx;
        adhesionTracks(kept).pos        = [xs(:), ys(:)];
        adhesionTracks(kept).amp        = as(:)';
        adhesionTracks(kept).dist       = dist;
        adhesionTracks(kept).dist_nm    = dist * pixelSize;
        adhesionTracks(kept).isInterp   = isInterp;
        adhesionTracks(kept).lifetime   = nF;
        adhesionTracks(kept).lifetime_s = (fr(end)-fr(1)) * timeInterval;
    end
end
end