
% Francois Aguet, 08/2013

function [fnames, fpaths, mexNames, mexPaths, sourceNames, sourcePaths, ignoreList] = parseFilesForExport(fctList, ignoreList, externList)

[fpaths, fnames, fexts] = cellfun(@fileparts, fctList, 'unif', 0);

if ~isempty(ignoreList)
    % remove unneeded functions
    rmIdx = ismember(fnames, ignoreList);
    if ~isempty(externList)
        tmp = cellfun(@(i) regexpi(fpaths, i, 'once'), externList, 'unif', 0);
        tmp = horzcat(tmp{:});
        rmIdx = rmIdx | any(~cellfun(@(i) isempty(i), tmp),2);
    end
    fpaths(rmIdx) = [];
    fnames(rmIdx) = [];
    fexts(rmIdx) = [];
end

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
searchPaths = cellfun(@(i,j) [i filesep j '*'], mexPaths, fnames(mexIdx), 'unif', 0);
sourceNames = cellfun(@dir, searchPaths, 'unif', 0);
sourceNames = cellfun(@(i) {i(~cellfun(@isempty, regexpi({i.name}, '(\.c(pp)?|\.h(pp)?)$'))).name}, sourceNames, 'unif', 0);
sourcePaths = arrayfun(@(i) repmat(mexPaths(i), [numel(sourceNames{i}) 1]), 1:numel(mexPaths), 'unif', 0);
sourcePaths = vertcat(sourcePaths{:});
sourceNames = horzcat(sourceNames{:})';

mexPaths = reshape(repmat(mexPaths', [numel(mexExts) 1]), [numel(mexExts)*numel(mexPaths) 1]);
fpaths(mexIdx) = [];
fnames(mexIdx) = [];
fexts(mexIdx) = [];
fnames = arrayfun(@(i) [fnames{i} fexts{i}], 1:numel(fnames), 'unif', 0)';