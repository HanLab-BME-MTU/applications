

function tracks = cutTrack(track, cutIdx)

np = numel(track.t);
ns = numel(cutIdx)+1;

cutGaps = cellfun(@(i) any(ismember(i, cutIdx)), track.gapIdx);

ub = [cellfun(@(i) i(1)-1, track.gapIdx(cutGaps)) np];
lb = [1 cellfun(@(i) i(end)+1, track.gapIdx(cutGaps))];


fnames = {'t', 'f', 'x', 'y', 'A', 'c', 'x_pstd', 'y_pstd', 'A_pstd', 'c_pstd',...
    'sigma_r', 'SE_sigma_r', 'pval_Ar', 'isPSF', 'tracksFeatIndxCG', 'gapVect',...
    'maskA', 'maskN', 'mask_Ar', 'hval_Ar', 'hval_AD'};

bnames = {'x', 'y', 'A', 'c', 'A_pstd', 'sigma_r', 'SE_sigma_r', 'pval_Ar'};

framerate = track.t(2)-track.t(1);

tracks(1:ns) = struct(track); 
for s = 1:ns
    range = lb(s):ub(s);
    for f = 1:numel(fnames)
        tracks(s).(fnames{f}) = track.(fnames{f})(:,range);
    end
    tracks(s).start = tracks(s).f(1);
    tracks(s).end = tracks(s).f(end);
    tracks(s).lifetime_s = (tracks(s).end-tracks(s).start+1)*framerate;
    tracks(s).seqOfEvents = [tracks(s).start 1 1 NaN; tracks(s).end 2 1 NaN];

    
    if s==1
        tracks(s).startBuffer = track.startBuffer;
    else
        nb = numel(track.startBuffer.t);
        tracks(s).startBuffer.t = tracks(s).t(1) + (-nb:-1)*framerate;
        for f = 1:numel(bnames)
            % within track
            nbIn = min(cutIdx(s-1),nb);
            nbOut = nb-nbIn;
            tracks(s).startBuffer.(bnames{f})(nbOut+1:nb) = track.(bnames{f})(:,cutIdx(s-1)-(nbIn-1:-1:0));
            % within parent start buffer
            if nbOut>0
                tracks(s).startBuffer.(bnames{f})(1:nbOut) = track.startBuffer.(bnames{f})(:,end-nbOut+1:end);
            end
        end
    end
    
    if s==ns
        tracks(s).endBuffer = track.endBuffer;
    else
        nb = numel(track.endBuffer.t);
        tracks(s).endBuffer.t = tracks(s).t(end) + (1:nb)*framerate;
        for f = 1:numel(bnames)
            tracks(s).endBuffer.(bnames{f}) = track.(bnames{f})(:,cutIdx(s)+(0:nb-1));
        end
    end
    % gaps in current segment
    idx = cellfun(@(i) all(lb(s)<=i & i<=ub(s)), track.gapIdx);
    tracks(s).gapStatus = track.gapStatus(idx);
    tracks(s).gapIdx = track.gapIdx(idx);
end