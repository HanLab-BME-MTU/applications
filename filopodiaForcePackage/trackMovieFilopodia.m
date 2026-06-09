function trackMovieFilopodia(movieData)
%TRACKMOVIEFILOPODIA  Process 3 wrapper (skeleton).
% Loads filoInfo (P2), links tip(+base) across frames by LAP, gap-closes,
% then derives L(t), signed dL/dt, state and fluctuation frequency.
% Sangyoon J. Han / 2026

iProc = movieData.getProcessIndex('FilopodiaTrackingProcess', 1, 0);
proc  = movieData.processes_{iProc};
p     = parseProcessParams(proc);
iChan = p.ChannelIndex;

iDet = p.DetProcessIndex;
if isempty(iDet), iDet = movieData.getProcessIndex('FilopodiaDetectionProcess', 1, 0); end
detProc = movieData.processes_{iDet};

inFilePaths = cell(1, numel(movieData.channels_));
inFilePaths{1, iChan} = detProc.outFilePaths_{1, iChan};
proc.setInFilePaths(inFilePaths);

outDir = p.OutputDirectory; mkClrDir(outDir);
outFile = fullfile(outDir, 'filoTracks.mat');
outFilePaths = cell(1, numel(movieData.channels_));
outFilePaths{1, iChan} = outFile;
proc.setOutFilePaths(outFilePaths);

filoInfo = detProc.loadChannelOutput(iChan, 'output', 'filoInfo'); %#ok<NASGU>

% TODO ------------------------------------------------------------------
%  Frame-to-frame: build tip cost matrix (gated by MaxLinkDist; optionally
%  add base-displacement term when LinkUseBase), solve with lap.m, then gap
%  close up to MaxGapFrames. For each surviving track:
%    L(t)   = arc(end) of the matched filoInfo centerline (geodesic base->tip)
%    vel    = smooth(diff(L), VelSmoothWin)         % + protrusion, - retraction
%    state  = sign(vel) thresholded by PauseThreshVel
%    fluctFreq = dominant frequency of detrended L(t) per FreqMethod
%  Assemble filoTracks struct array (see data-structures doc).
filoTracks = struct([]);   % placeholder
% ----------------------------------------------------------------------

save(outFile, 'filoTracks', '-v7.3');
disp('Filopodia tracking done (skeleton).');
end
