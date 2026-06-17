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
    'FilopodiaSegmentationProcess', 'FilopodiaSegmentation';
    'FilopodiaDetectionProcess',    'FilopodiaDetection';
    'FilopodiaTrackingProcess',     'FilopodiaTracking';
    'FilopodiaWindowingProcess',    'FilopodiaWindowing';
    'FilopodiaSamplingProcess',     'FilopodiaSampling';
    'FilopodiaStatisticsProcess',   'FilopodiaStatistics';
};

% correct OutputDirectory for every Filopodia process already on MD
corrected = 0;
for i = 1:numel(MD.processes_)
    proc = MD.processes_{i};
    if isempty(proc), continue; end
    row = find(strcmp(procMap(:,1), class(proc)));
    if isempty(row), continue; end
    correctPath = fullfile(root, procMap{row,2});
    if ~strcmp(proc.funParams_.OutputDirectory, correctPath)
        pp = proc.funParams_;
        pp.OutputDirectory = correctPath;
        proc.setPara(pp);
        corrected = corrected + 1;
        fprintf('  Corrected OutputDirectory for %s\n', class(proc));
    end
end
if corrected > 0
    MD.save();
    fprintf('Saved MD with %d corrected path(s).\n', corrected);
end

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