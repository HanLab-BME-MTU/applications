function computeMovieFilopodiaStats(movieData)
%COMPUTEMOVIEFILOPODIASTATS  Process 6 wrapper (skeleton).
% Aggregates tracks/windows/samples into per-movie statistics and correlations.
% Sangyoon J. Han / 2026

iProc = movieData.getProcessIndex('FilopodiaStatisticsProcess', 1, 0);
proc  = movieData.processes_{iProc};
p     = parseProcessParams(proc);
iChan = p.ChannelIndex;

iTrk = movieData.getProcessIndex('FilopodiaTrackingProcess', 1, 0);
iWin = movieData.getProcessIndex('FilopodiaWindowingProcess', 1, 0);
iSmp = movieData.getProcessIndex('FilopodiaSamplingProcess', 1, 0);
trkProc = movieData.processes_{iTrk};
winProc = movieData.processes_{iWin};
smpProc = movieData.processes_{iSmp};

inFilePaths = cell(1, numel(movieData.channels_));
inFilePaths{1, iChan} = trkProc.outFilePaths_{1, iChan};
proc.setInFilePaths(inFilePaths);

outDir = p.OutputDirectory; mkClrDir(outDir);
outFile = fullfile(outDir, 'filoStats.mat');
outFilePaths = cell(1, numel(movieData.channels_));
outFilePaths{1, iChan} = outFile;
proc.setOutFilePaths(outFilePaths);

filoTracks  = trkProc.loadChannelOutput(iChan, 'output', 'filoTracks');  %#ok<NASGU>
filoWindows = winProc.loadChannelOutput(iChan, 'output', 'filoWindows'); %#ok<NASGU>
filoSamples = smpProc.loadChannelOutput(iChan, 'output', 'filoSamples'); %#ok<NASGU>

% TODO ------------------------------------------------------------------
%  count per frame & mean; length / velocity / fluctuation-frequency dists;
%  lifetime; region-resolved force & talin dists (tip/shaft/base);
%  correlations in p.Correlations; density vs body perimeter/area.
%  Optionally export CSV (p.ExportCSV) and figures (p.MakeFigures).
filoStats = struct();   % placeholder
% ----------------------------------------------------------------------

save(outFile, 'filoStats', '-v7.3');
disp('Filopodia statistics done (skeleton).');
end
