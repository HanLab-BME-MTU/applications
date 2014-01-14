%exportEndocytosisPackage(exportDir) exports the CME analysis software package to 'exportDir'
%
% Adapted from buildPackage.m

% Francois Aguet, 05/28/2013

function exportEndocytosisPackage(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('exportDir', '.');
ip.addParameter('IncludeSources', true, @islogical);
ip.addParameter('Compress', false, @islogical);
ip.parse(varargin{:});
exportDir = ip.Results.exportDir;

% list of variable names erroneously detected as functions
ignoreList = {'consist'; 'knn'};

% packages to include from extern
% externRoot = [getParentDir(mfilename('fullpath'),3) 'extern' filesep];
% externList = {'kdtree_andrea'};
% externName = {'kdtree'};

%-----------------------------------------------------
% 1) Export pointSourceDetection.m to sub-directory
%-----------------------------------------------------
dest = [exportDir filesep 'cmeAnalysisPackage' filesep];

fnames = exportCode('PointSourceDetection.m', [dest 'PointSourceDetection' filesep],...
    'IncludeSources', ip.Results.IncludeSources, 'AddList', {'stats.h', 'KDTree.hpp'},...
    'IncludeGPL', false, 'Verbose', false);
ignoreList = [ignoreList; fnames];

% copy external packages
% for i = 1:numel(externList)
%     edest = [dest externName{i}];
%     source = [externRoot externList{i}];
%     %cmd = ['svn export ' source ' ' edest];
%     %system(cmd);
%     % -or-
%     [~,~] = mkdir(edest);
%     copyfile(source, edest);
% end

%-----------------------------------------------------
% 2) Export cmeAnalysis package
%-----------------------------------------------------
ts = loadTrackSettings();
masterList = {'cmeAnalysis.m', 'cmeDataViewer.m',...
    'plotPSNRDistribution.m', 'plotLifetimeComparison.m',...
    'stk2tiffDirs.m', 'selectConditionData.m',...
    'listSelectGUI.m', ...
    ts.costMatrices(1).funcName, ts.costMatrices(2).funcName,...
    ts.kalmanFunctions.reserveMem, ts.kalmanFunctions.initialize,...
    ts.kalmanFunctions.calcGain, ts.kalmanFunctions.timeReverse};

exportCode(masterList, dest,...
    'IncludeSources', false, 'AddList', [], 'IgnoreList', ignoreList);

% copy GUIs (temp fix)
copyfile(which('listSelectGUI.fig'), [dest 'listSelectGUI.fig']);

% set permissions
cmd = ['chmod -R 777 ' dest];
system(cmd);

% Copy documentation
docname = 'CME Analysis Package Manual.pdf';
copyfile(which(docname), [exportDir filesep docname]);


if ip.Results.Compress
    system(['zip -qr cmeAnalysisPackage.zip cmeAnalysisPackage "' docname '"']);
    
    % system('cp cmeAnalysisPackage.zip www/doc/.');
    %unzip cmeAnalysisPackage.zip
    
    % Remove temporary files
    system(['rm "' docname '"']);
    system('rm -r cmeAnalysisPackage');
end

disp('Export complete.');
