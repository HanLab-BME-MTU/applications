function windowMovieFilopodia(movieData)
%WINDOWMOVIEFILOPODIA  Process 4 wrapper (skeleton).
% Lays arclength windows on each tracked centerline; tip window re-anchored to
% the current distal endpoint per frame; stores narrow (intensity) and wide
% (force) footprints plus the local tangent.
% Sangyoon J. Han / 2026

iProc = movieData.getProcessIndex('FilopodiaWindowingProcess', 1, 0);
proc  = movieData.processes_{iProc};
p     = parseProcessParams(proc);
iChan = p.ChannelIndex;

iTrk = p.TrackProcessIndex;
if isempty(iTrk), iTrk = movieData.getProcessIndex('FilopodiaTrackingProcess', 1, 0); end
trkProc = movieData.processes_{iTrk};

inFilePaths = cell(1, numel(movieData.channels_));
inFilePaths{1, iChan} = trkProc.outFilePaths_{1, iChan};
proc.setInFilePaths(inFilePaths);

outDir = p.OutputDirectory; mkClrDir(outDir);
outFile = fullfile(outDir, 'filoWindows.mat');
outFilePaths = cell(1, numel(movieData.channels_));
outFilePaths{1, iChan} = outFile;
proc.setOutFilePaths(outFilePaths);

filoTracks = trkProc.loadChannelOutput(iChan, 'output', 'filoTracks'); %#ok<NASGU>

% TODO ------------------------------------------------------------------
%  For each track/frame: resample centerline by arclength (NumWindows in
%  'normalized' mode, or WindowLength px in 'fixed' mode). For each window
%  store winIndex, region (BaseFraction/TipFraction), sCenter, unit tangent,
%  narrow footprint (LateralWidth) and wide footprint (ForceLateralWidth or
%  auto from the TFM regularization scale). Build filoWindows.
filoWindows = struct([]);   % placeholder
% ----------------------------------------------------------------------

save(outFile, 'filoWindows', '-v7.3');
disp('Filopodia windowing done (skeleton).');
end
