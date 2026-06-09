function sampleMovieFilopodia(movieData)
%SAMPLEMOVIEFILOPODIA  Process 5 wrapper (skeleton).
% Samples per-window traction (magnitude + axial/lateral) and talin intensity.
% Traction is read from the TFMPackage's ForceFieldCalculationProcess.
% Sangyoon J. Han / 2026

iProc = movieData.getProcessIndex('FilopodiaSamplingProcess', 1, 0);
proc  = movieData.processes_{iProc};
p     = parseProcessParams(proc);
iChan = p.ChannelIndex;

iWin = p.WindowProcessIndex;
if isempty(iWin), iWin = movieData.getProcessIndex('FilopodiaWindowingProcess', 1, 0); end
winProc = movieData.processes_{iWin};

% --- cross-package: locate the traction field ---
iPkg = movieData.getPackageIndex(p.ForcePackageName);
assert(~isempty(iPkg), 'TFMPackage not found on this MovieData.');
tfmPack = movieData.getPackage(iPkg);
iForce = tfmPack.getProcessIndex(p.ForceProcessName);
assert(~isempty(iForce), 'ForceFieldCalculationProcess not found / not run.');
forceProc = tfmPack.getProcess(iForce); %#ok<NASGU>
% Per frame, fetch the traction magnitude map, e.g.:
%   tMap = forceProc.loadChannelOutput('output', p.ForceOutput, 'iFrame', t);

inFilePaths = cell(1, numel(movieData.channels_));
inFilePaths{1, iChan} = winProc.outFilePaths_{1, iChan};
proc.setInFilePaths(inFilePaths);

outDir = p.OutputDirectory; mkClrDir(outDir);
outFile = fullfile(outDir, 'filoSamples.mat');
outFilePaths = cell(1, numel(movieData.channels_));
outFilePaths{1, iChan} = outFile;
proc.setOutFilePaths(outFilePaths);

filoWindows = winProc.loadChannelOutput(iChan, 'output', 'filoWindows'); %#ok<NASGU>

% TODO ------------------------------------------------------------------
%  Apply ChannelShift (green<->red registration). For each window/frame:
%    forceMag     = SampleStat over pixIdxForce on |traction|
%    forceAxial   = mean(traction . tangent)
%    forceLateral = mean(traction . normal)
%    talinInt     = SampleStat over pixIdxInt on the talin channel
%  Build filoSamples indexed [trackId, frame, winIndex].
filoSamples = struct([]);   % placeholder
% ----------------------------------------------------------------------

save(outFile, 'filoSamples', '-v7.3');
disp('Filopodia sampling done (skeleton).');
end
