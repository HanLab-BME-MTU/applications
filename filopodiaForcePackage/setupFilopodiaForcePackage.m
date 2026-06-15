function pkg = setupFilopodiaForcePackage(MD)
%SETUPFILOPODIAFORCEPACKAGE  Register FilopodiaForcePackage on a MovieData.
%
%   pkg = setupFilopodiaForcePackage(MD)
%
% Checks if FilopodiaForcePackage is already registered on MD. If not,
% creates and registers it, then saves MD. Returns the package object.
%
% Run once per MovieData before using the package. After this,
% MD.packages_ will include the package and it will appear in the GUI.
%
% Example:
%   MD = MovieData.load('/path/to/movieData.mat');
%   pkg = setupFilopodiaForcePackage(MD);
%   pkg.processes_{1}.run();   % run segmentation
%
% Sangyoon J. Han / 2026

% check if already registered
iP = cellfun(@(p) isa(p,'FilopodiaForcePackage'), MD.packages_);
if any(iP)
    pkg = MD.packages_{find(iP,1)};
    fprintf('FilopodiaForcePackage already registered (index %d).\n', find(iP,1));
    return;
end

% sanity: needs pixel size and time interval for physical units
assert(~isempty(MD.pixelSize_), ...
    'MD.pixelSize_ is empty. Fill in pixel size before setting up the package.');
assert(~isempty(MD.timeInterval_), ...
    'MD.timeInterval_ is empty. Fill in time interval before setting up the package.');

% create and register
pkg = FilopodiaForcePackage(MD);
MD.addPackage(pkg);
MD.save();

fprintf('FilopodiaForcePackage registered and MD saved.\n');
fprintf('Output root: %s\n', pkg.outputDirectory_);
fprintf('Processes:\n');
for i = 1:numel(pkg.processes_)
    if ~isempty(pkg.processes_{i})
        fprintf('  %d. %s\n', i, class(pkg.processes_{i}));
    else
        fprintf('  %d. (not yet created)\n', i);
    end
end
end
