function [] = filopodiaRun(MDPath, processesToRun, forceRun, iChanTal)
%FILOPODIARUN  Run the FilopodiaForcePackage (P1-P6) on one movie, mirroring
%tfmRun: each process runs only if it has not succeeded, is not updated, has
%changed, or forceRun is set. success_ is saved to disk after each process.
%
%   filopodiaRun(MDPath)                         % run whatever is not done
%   filopodiaRun(MDPath, [], false)
%   filopodiaRun(MDPath, 6, true)                % force-rerun only P6
%   filopodiaRun(MDPath, 4:6, false, 1)          % run P4-P6, channel 1
%
% processesToRun indexes the package slots 1..6:
%   1 Segmentation  2 Detection  3 Tracking
%   4 Classification 5 Sampling  6 Statistics
% Sangyoon J. Han / 2026

%% set up
if isa(MDPath,'MovieData')
    MD = MovieData.load(MDPath.getFullPath);
else
    MD = MovieData.load(MDPath,'askUser',false);
end

if nargin < 4 || isempty(iChanTal), iChanTal = 1; end
if nargin < 3 || isempty(forceRun),  forceRun  = false; end

% make sure the package + all six processes are registered with default
% (debug) settings; this is idempotent and does not re-run anything.
setupFilopodiaForceDefaults(MD, iChanTal);

iPack = MD.getPackageIndex('FilopodiaForcePackage');
pkg   = MD.getPackage(iPack);
status = pkg.sanityCheck;
%% which processes to run
if nargin < 2 || isempty(processesToRun)
    processesToRun = find(~status);
end

for ii = processesToRun
    curProcess = pkg.getProcess(ii);
    if isempty(curProcess), continue; end
    if ~curProcess.success_ || ~curProcess.updated_ || curProcess.procChanged_ || forceRun
        fprintf('  P%d %s: running...\n', ii, class(curProcess));
        curProcess.run
        MD.save
    else
        fprintf('  P%d %s: up to date, skip\n', ii, class(curProcess));
    end
end
end

% =====================================================================
function setupFilopodiaForceDefaults(MD, iChanTal)
% Register package + processes with the runFilopodiaForce_debug settings.
% Only sets parameters on processes that are not yet registered, so existing
% runs (and their success_ flags) are preserved. setupFilopodiaForcePackage
% is called LAST so every freshly added process gets wired into a slot.
tFrame = [];

% ---- P1 ----
if isempty(MD.getProcessIndex('FilopodiaSegmentationProcess',1,0))
    MD.addProcess(FilopodiaSegmentationProcess(MD));
    pr = MD.processes_{MD.getProcessIndex('FilopodiaSegmentationProcess',1,0)};
    pp = pr.funParams_;
    pp.ChannelIndex=iChanTal; pp.SteerableOrder=4; pp.SigmaArray=[1 2];
    pp.BodyThreshold='otsu'; pp.GaussianBlurSigma=2;
    pp.BodyOpenRadius=8; pp.BodyClosingRadius=8; pp.ProcessFrames=tFrame;
    pr.setPara(pp);
end
% ---- P2 ----
if isempty(MD.getProcessIndex('FilopodiaDetectionProcess',1,0))
    MD.addProcess(FilopodiaDetectionProcess(MD));
    pr = MD.processes_{MD.getProcessIndex('FilopodiaDetectionProcess',1,0)};
    defP2 = FilopodiaDetectionProcess.getDefaultParams(MD);
    pp = pr.funParams_; fn=fieldnames(defP2);
    for i=1:numel(fn), if ~isfield(pp,fn{i}), pp.(fn{i})=defP2.(fn{i}); end, end
    pp.ChannelIndex=iChanTal; pp.PSFsigma=2.1; pp.Alpha=0.05;
    pp.TipMaxDistFromBody=70; pp.MaxTipBaseDist=90; pp.BaseInsideBand=4;
    pp.DetectMode='auto'; pp.ProcessFrames=tFrame;
    pr.setPara(pp);
end
% ---- P3 ----
if isempty(MD.getProcessIndex('FilopodiaTrackingProcess',1,0))
    MD.addProcess(FilopodiaTrackingProcess(MD));
    pr = MD.processes_{MD.getProcessIndex('FilopodiaTrackingProcess',1,0)};
    pp = pr.funParams_;
    pp.MaxLinkDist=8; pp.MaxGapFrames=3; pp.MinTrackLength=3;
    pp.OutputDirectory=fullfile(MD.outputDirectory_,'FilopodiaForcePackage','FilopodiaTracking');
    pr.setPara(pp);
end
% ---- P4 ----
if isempty(MD.getProcessIndex('FilopodiaClassificationProcess',1,0))
    MD.addProcess(FilopodiaClassificationProcess(MD));
    pr = MD.processes_{MD.getProcessIndex('FilopodiaClassificationProcess',1,0)};
    pp = pr.funParams_;
    pp.MinTipLifetime=5; pp.MinLinearFrac=0.85; pp.MinTipDist=6;
    pp.ShaftBand=4; pp.BodyMaxAngle=60; pp.MinReachFrac=0.5; pp.LenPenalty=0.6;
    pp.OutputDirectory=fullfile(MD.outputDirectory_,'FilopodiaForcePackage','FilopodiaClassification');
    pr.setPara(pp);
end
% ---- P5 ----
if isempty(MD.getProcessIndex('FilopodiaSamplingProcess',1,0))
    MD.addProcess(FilopodiaSamplingProcess(MD));
    pr = MD.processes_{MD.getProcessIndex('FilopodiaSamplingProcess',1,0)};
    pp = pr.funParams_;
    pp.ChannelIndex=iChanTal; pp.ShaftSampleStep=3; pp.SampleRadius=1;
    pp.OutputDirectory=fullfile(MD.outputDirectory_,'FilopodiaForcePackage','FilopodiaSampling');
    pr.setPara(pp);
end
% ---- P6 ----
if isempty(MD.getProcessIndex('FilopodiaStatisticsProcess',1,0))
    MD.addProcess(FilopodiaStatisticsProcess(MD));
    pr = MD.processes_{MD.getProcessIndex('FilopodiaStatisticsProcess',1,0)};
    pp = pr.funParams_;
    pp.ChannelIndex=iChanTal; pp.MinLifetimeForStats=3; pp.MakeFigures=true;
    pp.OutputDirectory=fullfile(MD.outputDirectory_,'FilopodiaForcePackage','FilopodiaStatistics');
    pr.setPara(pp);
end

% wire everything (package registration + slot linking) now that all six
% processes exist on MD.
setupFilopodiaForcePackage(MD);
end
