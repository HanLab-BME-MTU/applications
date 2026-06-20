function runFilopodiaForceOnMD(MD, iChanTal, varargin)
%RUNFILOPODIAFORCEONMD  Run the full FilopodiaForcePackage (P1-P6) on one
%MovieData, using the default settings from runFilopodiaForce_debug (no QC,
%no figures except P6 normalized-profile plots).
%
%   runFilopodiaForceOnMD(MD)                 % iChanTal defaults to 1
%   runFilopodiaForceOnMD(MD, iChanTal)
%   runFilopodiaForceOnMD(MD, iChanTal, 'Overwrite', true)
%
% Requires the TFMPackage ForceFieldCalculationProcess to have already been
% run on MD (P5 reads forceField.mat). If it is missing, P5 traction is NaN.
% Sangyoon J. Han / 2026

if nargin < 2 || isempty(iChanTal), iChanTal = 1; end
ip = inputParser;
ip.addParameter('Overwrite', false, @islogical);
ip.parse(varargin{:});
overwrite = ip.Results.Overwrite;

tFrame = [];   % [] = all frames

MD.sanityCheck;
fprintf('  Movie: %d ch, %d frames, pix=%.4g nm, dt=%.4g s\n', ...
    numel(MD.channels_), MD.nFrames_, MD.pixelSize_, MD.timeInterval_);

% register package (idempotent)
setupFilopodiaForcePackage(MD);

%% P1: SEGMENTATION
segProc = getOrAddProc(MD,'FilopodiaSegmentationProcess');
if needRun(segProc, overwrite)
    pp = segProc.funParams_;
    pp.ChannelIndex      = iChanTal;
    pp.SteerableOrder    = 4;
    pp.SigmaArray        = [1 2];
    pp.BodyThreshold     = 'otsu';
    pp.GaussianBlurSigma = 2;
    pp.BodyOpenRadius    = 8;
    pp.BodyClosingRadius = 8;
    pp.ProcessFrames     = tFrame;
    segProc.setPara(pp);
    tic; segProc.run(); fprintf('  P1 segmentation %.0fs\n',toc);
else
    fprintf('  P1 already done, skipping.\n');
end

%% P2: DETECTION
detProc = getOrAddProc(MD,'FilopodiaDetectionProcess');
if needRun(detProc, overwrite)
    defP2 = FilopodiaDetectionProcess.getDefaultParams(MD);
    pp2 = detProc.funParams_;
    fn = fieldnames(defP2);
    for ii = 1:numel(fn)
        if ~isfield(pp2, fn{ii}), pp2.(fn{ii}) = defP2.(fn{ii}); end
    end
    pp2.ChannelIndex       = iChanTal;
    pp2.PSFsigma           = 2.1;
    pp2.Alpha              = 0.05;
    pp2.TipMaxDistFromBody = 70;
    pp2.MaxTipBaseDist     = 90;
    pp2.BaseInsideBand     = 4;
    pp2.DetectMode         = 'auto';
    pp2.ProcessFrames      = tFrame;
    detProc.setPara(pp2);
    tic; detProc.run(); fprintf('  P2 detection %.0fs\n',toc);
else
    fprintf('  P2 already done, skipping.\n');
end

% confirm detectMode is 'all' (required for P3-P6)
detFile = fullfile(detProc.funParams_.OutputDirectory,'filoDetection.mat');
S = load(detFile,'detectMode');
if ~strcmp(S.detectMode,'all')
    warning('detectMode=%s (not ''all''); P3-P6 require ''all''. Skipping rest of this movie.', S.detectMode);
    return;
end

%% P3: TRACKING
trkProc = getOrAddProc(MD,'FilopodiaTrackingProcess');
if needRun(trkProc, overwrite)
    pp3 = trkProc.funParams_;
    pp3.MaxLinkDist     = 8;
    pp3.MaxGapFrames    = 3;
    pp3.MinTrackLength  = 3;
    pp3.OutputDirectory = fullfile(MD.outputDirectory_,'FilopodiaForcePackage','FilopodiaTracking');
    trkProc.setPara(pp3);
    tic; trkProc.run(); fprintf('  P3 tracking %.0fs\n',toc);
else
    fprintf('  P3 already done, skipping.\n');
end

%% P4: CLASSIFICATION
clsProc = getOrAddProc(MD,'FilopodiaClassificationProcess');
if needRun(clsProc, overwrite)
    pp4 = clsProc.funParams_;
    pp4.MinTipLifetime  = 5;
    pp4.MinLinearFrac   = 0.85;
    pp4.MinTipDist      = 6;
    pp4.ShaftBand       = 4;
    pp4.BodyMaxAngle    = 60;
    pp4.MinReachFrac    = 0.5;
    pp4.LenPenalty      = 0.6;
    pp4.OutputDirectory = fullfile(MD.outputDirectory_,'FilopodiaForcePackage','FilopodiaClassification');
    clsProc.setPara(pp4);
    tic; clsProc.run(); fprintf('  P4 classification %.0fs\n',toc);
else
    fprintf('  P4 already done, skipping.\n');
end

%% P5: FORCE / INTENSITY SAMPLING
smpProc = getOrAddProc(MD,'FilopodiaSamplingProcess');
if needRun(smpProc, overwrite)
    pp5 = smpProc.funParams_;
    pp5.ChannelIndex     = iChanTal;
    pp5.ShaftSampleStep  = 3;
    pp5.SampleRadius     = 1;
    pp5.OutputDirectory  = fullfile(MD.outputDirectory_,'FilopodiaForcePackage','FilopodiaSampling');
    smpProc.setPara(pp5);
    tic; smpProc.run(); fprintf('  P5 sampling %.0fs\n',toc);
else
    fprintf('  P5 already done, skipping.\n');
end

%% P6: STATISTICS
statProc = getOrAddProc(MD,'FilopodiaStatisticsProcess');
if needRun(statProc, overwrite)
    pp6 = statProc.funParams_;
    pp6.ChannelIndex        = iChanTal;
    pp6.MinLifetimeForStats = 3;
    pp6.MakeFigures         = true;
    pp6.OutputDirectory     = fullfile(MD.outputDirectory_,'FilopodiaForcePackage','FilopodiaStatistics');
    statProc.setPara(pp6);
    tic; statProc.run(); fprintf('  P6 statistics %.0fs\n',toc);
else
    fprintf('  P6 already done, skipping.\n');
end

end

% =====================================================================
function proc = getOrAddProc(MD, procName)
i = MD.getProcessIndex(procName,1,0);
if isempty(i)
    eval(sprintf('MD.addProcess(%s(MD));', procName));
    i = MD.getProcessIndex(procName,1,0);
end
proc = MD.processes_{i};
end

% =====================================================================
function tf = needRun(proc, overwrite)
if overwrite, tf = true; return; end
% run unless a non-empty output file already exists
tf = true;
try
    op = proc.outFilePaths_;
    for k = 1:numel(op)
        if ~isempty(op{k}) && exist(op{k},'file')==2, tf = false; return; end
    end
catch
end
end