function pkg = setupFilopodiaForcePackage(MD)
%SETUPFILOPODIAFORCEPACKAGE  Register FilopodiaForcePackage on a MovieData.
%
%   pkg = setupFilopodiaForcePackage(MD)
%
% Checks if FilopodiaForcePackage is already registered on MD. If not,
% creates and registers it, then saves MD. Either way, corrects the
% OutputDirectory of every Filopodia process to the canonical path under
% <movieDir>/FilopodiaForcePackage/. This is necessary when processes
% were first created before the folder convention was set.
%
% Sangyoon J. Han / 2026

assert(~isempty(MD.pixelSize_),    'MD.pixelSize_ is empty.');
assert(~isempty(MD.timeInterval_), 'MD.timeInterval_ is empty.');

% --- register package if not present ---
iP = cellfun(@(p) isa(p,'FilopodiaForcePackage'), MD.packages_);
if any(iP)
    pkg = MD.packages_{find(iP,1)};
    fprintf('FilopodiaForcePackage already registered (index %d).\n', find(iP,1));
else
    pkg = FilopodiaForcePackage(MD);
    MD.addPackage(pkg);
    fprintf('FilopodiaForcePackage registered.\n');
end

% --- canonical output root ---
root = fullfile(MD.outputDirectory_, 'FilopodiaForcePackage');

% map: process class -> subfolder name
procMap = {
    'FilopodiaSegmentationProcess',   'FilopodiaSegmentation';
    'FilopodiaDetectionProcess',      'FilopodiaDetection';
    'FilopodiaTrackingProcess',       'FilopodiaTracking';
    'FilopodiaClassificationProcess', 'FilopodiaClassification';
    'FilopodiaSamplingProcess',       'FilopodiaSampling';
    'FilopodiaStatisticsProcess',     'FilopodiaStatistics';
};

% correct OutputDirectory and fill any newly added default fields for every
% Filopodia process already on MD
corrected = 0;
for i = 1:numel(MD.processes_)
    proc = MD.processes_{i};
    if isempty(proc), continue; end
    row = find(strcmp(procMap(:,1), class(proc)));
    if isempty(row), continue; end

    pp = proc.funParams_;
    changed = false;

    % (a) fill missing fields from current defaults (handles params added later)
    def = proc.getDefaultParams(MD);
    fn = fieldnames(def);
    for k = 1:numel(fn)
        if ~isfield(pp, fn{k})
            pp.(fn{k}) = def.(fn{k});
            changed = true;
        end
    end

    % (b) force canonical OutputDirectory
    correctPath = fullfile(root, procMap{row,2});
    if ~strcmp(pp.OutputDirectory, correctPath)
        pp.OutputDirectory = correctPath;
        changed = true;
    end

    if changed
        proc.setPara(pp);
        corrected = corrected + 1;
        fprintf('  Updated params/path for %s\n', class(proc));
    end

    % (c) restore outFilePaths_ so packageGUI shows completion status
    outFileMap = struct( ...
        'FilopodiaSegmentationProcess',   '', ...
        'FilopodiaDetectionProcess',      'filoDetection.mat', ...
        'FilopodiaTrackingProcess',       'filoTracks.mat', ...
        'FilopodiaClassificationProcess', 'filoClassification.mat', ...
        'FilopodiaSamplingProcess',       'filoSamples.mat', ...
        'FilopodiaStatisticsProcess',     'filoStats.mat');
    cls = class(proc);
    if isfield(outFileMap, cls) && ~isempty(outFileMap.(cls))
        outFile = fullfile(proc.funParams_.OutputDirectory, outFileMap.(cls));
        if exist(outFile,'file')==2
            try
                nChan = numel(MD.channels_);
                outPaths = cell(1, nChan);
                outPaths{1,1} = outFile;
                proc.setOutFilePaths(outPaths);
                % mark process as successfully completed so packageGUI shows check
                % success_ is SetAccess=protected; use our subclass markSuccess()
                if ismethod(proc,'markSuccess')
                    proc.markSuccess();
                end
                fprintf('  Restored completion status for %s\n', cls);
            catch ME
                fprintf('  Warning: could not restore status for %s: %s\n', cls, ME.message);
            end
        end
    end
end
if corrected > 0
    MD.save();
    fprintf('Saved MD with %d updated process(es).\n', corrected);
end

% --- wire MD processes into package slots ---
% pkg.processes_ starts as {[] [] [] [] [] []}; fill slots from MD.processes_
classNames = pkg.getProcessClassNames();
for i = 1:numel(MD.processes_)
    proc = MD.processes_{i};
    if isempty(proc), continue; end
    slot = find(strcmp(classNames, class(proc)), 1);
    if ~isempty(slot) && isempty(pkg.processes_{slot})
        try
            pkg.setProcess(slot, proc);
            fprintf('  Linked %s to package slot %d\n', class(proc), slot);
        catch ME
            fprintf('  Warning: could not link %s: %s\n', class(proc), ME.message);
        end
    end
end
MD.save();

% --- report ---
fprintf('Output root: %s\n', root);
for i = 1:numel(MD.processes_)
    proc = MD.processes_{i};
    if isempty(proc), continue; end
    if any(strcmp(procMap(:,1), class(proc)))
        fprintf('  %s -> %s\n', class(proc), proc.funParams_.OutputDirectory);
    end
end
end
