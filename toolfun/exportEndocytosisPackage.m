%exportEndocytosisPackage(exportDir) exports the CME analysis software package to 'exportDir'
%
% Adapted from buildPackage.m

% Francois Aguet, 05/28/2013

function exportEndocytosisPackage(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('exportDir', '.', @ischar);
ip.addParamValue('IncludeSources', true, @islogical);
ip.parse(varargin{:});
exportDir = ip.Results.exportDir;

% list of variable names erroneously detected as functions
ignoreList = {'consist'; 'knn'};

% packages to include from extern
externRoot = [getParentDir(mfilename('fullpath'),3) 'extern' filesep];
externList = {'kdtree_andrea'};
externName = {'kdtree'};

%-----------------------------------------------------
% 1) Export pointSourceDetection.m separately
%-----------------------------------------------------
[fctList, toolboxList] = getFunDependencies('pointSourceDetection.m');
[fnames, fpaths, mexNames, mexPaths, sourceNames, sourcePaths, ignoreList] = parseFiles(fctList, ignoreList, externList);

% copy core functions
dest = [exportDir filesep 'PointSourceDetection' filesep];
[~,~] = mkdir(dest);
for i = 1:numel(fnames)
    copyfile([fpaths{i} filesep fnames{i}], [dest fnames{i}]);
end

% copy MEX functions
mdest = [dest filesep 'mex' filesep];
[~,~] = mkdir(mdest);
for i = 1:numel(mexNames)
    copyfile([mexPaths{i} filesep mexNames{i}], [mdest mexNames{i}]);
end

if ip.Results.IncludeSources
    % copy source files (into MEX directory)
    for i = 1:numel(sourceNames)
        copyfile([sourcePaths{i} filesep sourceNames{i}], [mdest sourceNames{i}]);
    end
    % for now, include 'stats.h' explicitly
    copyfile('/home/fa48/matlab/SVN/common/mex/include/stats.h', [mdest 'stats.h']);
end

% copy external packages
for i = 1:numel(externList)
    edest = [dest externName{i}];
    source = [externRoot externList{i}];
    %cmd = ['svn export ' source ' ' edest];
    %system(cmd);
    % -or-
    [~,~] = mkdir(edest);
    copyfile(source, edest);
end

% set permissions
cmd = ['chmod -R 777 ' dest];
system(cmd);

% Copy GPL and package
copyfile('/home/fa48/matlab/SVN/endocytosis/doc/GPL-License.txt', './PointSourceDetection/.');
system('zip -qr PointSourceDetection.zip PointSourceDetection');
system('mv PointSourceDetection.zip www/aguet/doc/.');


%-----------------------------------------------------
% 2) Export Endocytosis project
%-----------------------------------------------------
ts = loadTrackSettings();
masterList = {'cmeAnalysis.m', 'cmeDataViewer.m',...
    'plotPSNRDistribution.m', 'plotLifetimeComparison.m',...
    'stk2tiffDirs.m', 'selectConditionData.m',...
    'listSelectGUI.m', ...
    ts.costMatrices(1).funcName, ts.costMatrices(2).funcName,...
    ts.kalmanFunctions.reserveMem, ts.kalmanFunctions.initialize,...
    ts.kalmanFunctions.calcGain, ts.kalmanFunctions.timeReverse};
for i = 1:numel(masterList);
    [fctList{i}, toolboxList{i}] = getFunDependencies(masterList{i});
end
fctList = unique(vertcat(fctList{:}));
toolboxList = unique(vertcat(toolboxList{:}));
[fnames, fpaths, mexNames, mexPaths, ~, ~, ~] = parseFiles(fctList, ignoreList, externList);

% copy core functions
dest = [exportDir filesep 'cmeAnalysisPackage' filesep];
[~,~] = mkdir(dest);
for i = 1:numel(fnames)
    copyfile([fpaths{i} filesep fnames{i}], [dest fnames{i}]);
end

% copy GUIs (temp fix)
copyfile('/home/fa48/matlab/SVN/common/iofun/GUI/listSelectGUI.fig',...
    './cmeAnalysisPackage/listSelectGUI.fig');

% copy MEX functions
mdest = [dest filesep 'mex' filesep];
[~,~] = mkdir(mdest);
for i = 1:numel(mexNames)
    copyfile([mexPaths{i} filesep mexNames{i}], [mdest mexNames{i}]);
end

% set permissions
cmd = ['chmod -R 777 ' dest];
system(cmd);

disp('The package uses the following toolboxes:')
disp(toolboxList);

% Copy documentation and GPL
copyfile('/home/fa48/matlab/SVN/endocytosis/doc/GPL-License.txt',...
    './cmeAnalysisPackage/.');
copyfile('/home/fa48/matlab/SVN/endocytosis/doc/CME Analysis Package Manual.pdf',...
    '.');

system('mv PointSourceDetection cmeAnalysisPackage/.');
% system('zip -qr cmeAnalysisPackage.zip cmeAnalysisPackage');
system('zip -qr cmeAnalysisPackage.zip cmeAnalysisPackage "CME Analysis Package Manual.pdf"');

system('cp cmeAnalysisPackage.zip www/doc/.');
%unzip cmeAnalysisPackage.zip


function [fnames, fpaths, mexNames, mexPaths, sourceNames, sourcePaths, ignoreList] = parseFiles(fctList, ignoreList, externList)

[fpaths, fnames, fexts] = cellfun(@fileparts, fctList, 'unif', 0);
% remove unneeded functions
rmIdx = ismember(fnames, ignoreList);
tmp = cellfun(@(i) regexpi(fpaths, i, 'once'), externList, 'unif', 0);
tmp = horzcat(tmp{:});
rmIdx = rmIdx | any(~cellfun(@(i) isempty(i), tmp),2);
fpaths(rmIdx) = [];
fnames(rmIdx) = [];
fexts(rmIdx) = [];

% Update ignore list for main package
ignoreList = [ignoreList; fnames];

mexExts = {'.mexmaci64'; '.mexa64'; '.mexw64'};
mexIdx = find(ismember(fexts, mexExts));
% Expand MEX list to all platforms
mexNames = fnames(mexIdx);
mexNames = arrayfun(@(i) cellfun(@(e) [mexNames{i} e], mexExts, 'unif', 0), 1:numel(mexNames), 'unif', 0);
mexNames = vertcat(mexNames{:});
mexPaths = fpaths(mexIdx);
% retrieve source files (.c, .h, .cpp & .hpp, but does not include external headers!)
sourceNames = cellfun(@dir, mexPaths, 'unif', 0);
sourceNames = cellfun(@(i) {i(~cellfun(@isempty, regexpi({i.name}, '(\.c(pp)?|\.h(pp)?)$'))).name}, sourceNames, 'unif', 0);
sourcePaths = arrayfun(@(i) repmat(mexPaths(i), [numel(sourceNames{i}) 1]), 1:numel(mexPaths), 'unif', 0);
sourcePaths = vertcat(sourcePaths{:});
sourceNames = vertcat(sourceNames{:});

mexPaths = reshape(repmat(mexPaths', [numel(mexExts) 1]), [numel(mexExts)*numel(mexPaths) 1]);
fpaths(mexIdx) = [];
fnames(mexIdx) = [];
fexts(mexIdx) = [];
fnames = arrayfun(@(i) [fnames{i} fexts{i}], 1:numel(fnames), 'unif', 0)';
