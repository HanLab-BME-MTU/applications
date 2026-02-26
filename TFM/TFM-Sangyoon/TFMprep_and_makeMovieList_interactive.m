%% TFMprep_and_makeMovieList_interactive.m
% GUI: conditions table + beadChan + strict Pos matching toggle.
%
% Assumption:
%   - Each ref file is a separate OME-TIFF (one ref per position).
%   - One ref file should correspond to one movie file (matching by Pos# token).
%
% Sangyoon Han

clear; clc;

%% ==== USER SETTINGS ====
baseOutName = "TFMprep_bestZ";   % output folder name under each movDir
metric      = "tenengrad";
verbose     = true;

% MovieList output location (per condition, inside each movDir by default)
movieListsRoot = ""; % "" means use: fullfile(movDir, baseOutName, label, "MovieList")

%% ==== GET CONDITIONS + beadChan + strictPosMatching (GUI preferred) ====
[conds, beadChan, strictPosMatching] = getInteractiveInputs(2, true); % default beadChan=2, strict=true
assert(~isempty(conds) && size(conds,2)==3, 'No conditions provided.');
fprintf('Using beadChan = %d\n', beadChan);
fprintf('Strict Pos matching = %s\n', string(strictPosMatching));

%% ==== MAIN LOOP ====
for i = 1:size(conds,1)
    label  = string(conds{i,1});
    movDir = string(conds{i,2});
    refDir = string(conds{i,3});

    fprintf('\n====================================================\n');
    fprintf('Condition %d/%d: %s\n', i, size(conds,1), label);
    fprintf('  movDir: %s\n', movDir);
    fprintf('  refDir: %s\n', refDir);
    fprintf('====================================================\n');

    assert(exist(movDir,'dir')==7, 'movDir not found: %s', movDir);
    assert(exist(refDir,'dir')==7, 'refDir not found: %s', refDir);

    % ---- find movie and ref OME-TIFFs ----
    mFiles = dir(fullfile(movDir, "*.ome.tif"));
    % ---- filter Micro-Manager split OME-TIFF chunks (Pos0, Pos0_1, Pos0_2...) ----
    [mFiles, chunkLog] = filterMMChunkedOmeTiffs(mFiles);
    
    if ~isempty(chunkLog)
        fprintf('[ChunkFilter] Collapsed split movies:\n');
        fprintf('  %s\n', chunkLog{:});
    end

    mFiles = mFiles(~contains({mFiles.name}, "metadata", 'IgnoreCase', true));
    assert(~isempty(mFiles), 'No movie OME-TIFF found in %s', movDir);

    rFiles = dir(fullfile(refDir, "*.ome.tif"));
    rFiles = rFiles(~contains({rFiles.name}, "metadata", 'IgnoreCase', true));
    assert(~isempty(rFiles), 'No ref OME-TIFF found in %s', refDir);

    % ---- collect all MovieData for this condition ----
    allMD = MovieData.empty;

    % ---- log mapping for debugging ----
    mapLog = strings(0,1);

    for k = 1:numel(mFiles)
        movieFile = fullfile(mFiles(k).folder, mFiles(k).name);

        % ---- match ONE ref file per movie file (strict or lenient) ----
        refHit = pickRefFileForMovie_singleRefPerPos(mFiles(k).name, rFiles, strictPosMatching);
        refFile = fullfile(refHit.folder, refHit.name);

        % ---- output dir per movie file ----
        baseName = erase(string(mFiles(k).name), ".ome.tif");
        outDir = fullfile(movDir, baseOutName, label, baseName);

        fprintf('\n=== %s | movie=%s ===\n', label, mFiles(k).name);
        fprintf('    ref = %s\n', string(refFile));
        fprintf('    out = %s\n', string(outDir));

        % ---- series matched extraction -> returns MovieData array (one per series) ----
        mdList = make2DTFMFromRefBestZ_seriesMatched(movieFile, refFile, outDir, ...
            'beadChan', beadChan, ...
            'metric', metric, ...
            'verbose', verbose);

        if isempty(mdList)
            warning('No MovieData created for %s', mFiles(k).name);
            continue;
        end

        allMD = [allMD; mdList(:)];

        mapLog(end+1,1) = sprintf("movie=%s | ref=%s | out=%s | nSeriesMD=%d", ...
            string(mFiles(k).name), string(refHit.name), outDir, numel(mdList));
    end

    fprintf('\n[Condition %s] Total MovieData objects created: %d\n', label, numel(allMD));
    assert(~isempty(allMD), 'No MovieData objects created for condition %s', label);

    % ---- movieList output directory ----
    if strlength(movieListsRoot) == 0
        mlDir = fullfile(movDir, baseOutName, label, "MovieList");
    else
        mlDir = fullfile(movieListsRoot, label);
    end
    if ~exist(mlDir,'dir'); mkdir(mlDir); end
    mlDir = char(mlDir);
    
    % ---- Deduplicate MovieData to avoid analyzing identical movies twice ----
    [allMD, dedupReport] = deduplicateMDArray(allMD);
    
    fprintf('[Condition %s] Dedup: nMD=%d (unique=%d, removed=%d)\n', ...
        label, dedupReport.nOriginal, dedupReport.nUnique, dedupReport.nRemoved);
    
    % Save dedup report
    dedupFile = fullfile(mlDir, 'TFMprep_dedup_report.txt');
    fid = fopen(dedupFile,'w');
    fprintf(fid, "Condition: %s\nOriginal nMD: %d\nUnique nMD: %d\nRemoved: %d\n\n", ...
        label, dedupReport.nOriginal, dedupReport.nUnique, dedupReport.nRemoved);
    
    for gi = 1:numel(dedupReport.groups)
        g = dedupReport.groups{gi};
        fprintf(fid, "Group %d: keep=%d remove=%s\n", gi, g.keep, mat2str(g.remove));
        fprintf(fid, "  key: %s\n", dedupReport.keys(g.keep));
        for jj = 1:numel(g.allIdx)
            ii = g.allIdx(jj);
            fprintf(fid, "    %d) %s\n", ii, dedupReport.details(ii));
        end
        fprintf(fid, "\n");
    end
    fclose(fid);
    fprintf('[Condition %s] Saved dedup report: %s\n', label, dedupFile);
    ML = MovieList(allMD, mlDir);
    ML.setFilename('movieList.mat');
    ML.sanityCheck;
    ML.save();

    fprintf('[Condition %s] Saved MovieList: %s\n', fullfile(mlDir,'movieList.mat'));

    % ---- save mapping log ----
    logFile = fullfile(mlDir, 'TFMprep_seriesMatched_mapLog.txt');
    fid = fopen(logFile,'w');
    fprintf(fid, "Condition: %s\nmovDir: %s\nrefDir: %s\nbeadChan: %d\nstrictPosMatching: %d\n\n", ...
        label, movDir, refDir, beadChan, strictPosMatching);
    for j=1:numel(mapLog)
        fprintf(fid, "%s\n", mapLog(j));
    end
    fclose(fid);
    fprintf('[Condition %s] Saved map log: %s\n', label, logFile);
end

disp('Done: TFMprep + MovieList creation for all conditions.');

%% ===================== Interactive GUI/CLI =====================
function [conds, beadChan, strictPosMatching] = getInteractiveInputs(defaultBeadChan, defaultStrict)
if nargin < 1 || isempty(defaultBeadChan), defaultBeadChan = 2; end
if nargin < 2 || isempty(defaultStrict), defaultStrict = true; end

conds = {};
beadChan = defaultBeadChan;
strictPosMatching = logical(defaultStrict);

useGUI = usejava('desktop') && ~isempty(ver('matlab'));

if useGUI
    try
        [conds, beadChan, strictPosMatching] = conditionTableGUI_withBeadChanAndStrict(defaultBeadChan, defaultStrict);
        if isempty(conds)
            error('No conditions returned from GUI.');
        end
        return;
    catch ME
        warning(ME.identifier,'GUI entry failed (%s). Falling back to CLI.',...
            ME.message);
    end
end

% ---- CLI fallback ----
fprintf('\n[Interactive setup] CLI mode\n');

beadChan = input(sprintf('Bead channel index (default %d): ', defaultBeadChan));
if isempty(beadChan); beadChan = defaultBeadChan; end
validateBeadChan(beadChan);

strictPosMatching = input(sprintf('Strict Pos matching? 1=yes, 0=no (default %d): ', defaultStrict));
if isempty(strictPosMatching); strictPosMatching = defaultStrict; end
strictPosMatching = logical(strictPosMatching);

n = input('How many conditions? (e.g., 4): ');
if isempty(n) || ~isscalar(n) || n<1
    error('Invalid number of conditions.');
end

conds = cell(n,3);
for i = 1:n
    fprintf('\n--- Condition %d/%d ---\n', i, n);
    label = input('Condition label (e.g., TNFa_BBS): ', 's');
    if isempty(label), error('Label cannot be empty.'); end

    movDir = pickFolderInteractive(sprintf('Select MOVIE folder for "%s"', label));
    refDir = pickFolderInteractive(sprintf('Select REF folder for "%s"', label));

    conds{i,1} = string(label);
    conds{i,2} = string(movDir);
    conds{i,3} = string(refDir);
end

beadChan = double(beadChan);
end

function validateBeadChan(beadChan)
if ~isscalar(beadChan) || isnan(beadChan) || beadChan < 1 || mod(beadChan,1)~=0
    error('Invalid beadChan. Must be a positive integer.');
end
end

function folder = pickFolderInteractive(dialogTitle)
if usejava('desktop')
    folder = uigetdir(pwd, dialogTitle);
    if isequal(folder,0)
        error('Folder selection canceled.');
    end
else
    folder = input([dialogTitle ' (type full path): '], 's');
    if isempty(folder)
        error('Folder path cannot be empty.');
    end
end
end

function [conds, beadChan, strictPosMatching] = conditionTableGUI_withBeadChanAndStrict(defaultBeadChan, defaultStrict)
f = uifigure('Name','TFMprep: Conditions + Settings','Position',[100 100 940 520]);

% beadChan
uilabel(f,'Text','Bead channel index (within series):','Position',[20 480 220 22]);
beadField = uieditfield(f,'numeric', ...
    'Position',[250 480 80 24], ...
    'Value', defaultBeadChan, ...
    'Limits',[1 Inf], ...
    'RoundFractionalValues','on');

% strict toggle
strictBox = uicheckbox(f, ...
    'Text','Strict Pos matching (require Pos# token and unique ref match)', ...
    'Position',[360 480 520 22], ...
    'Value', logical(defaultStrict));

% table
defaultData = {
    "control",  "", "";
    "TNFa",     "", "";
};
tbl = uitable(f, ...
    'Data', defaultData, ...
    'ColumnName', {'Condition label','Movie folder','Ref folder'}, ...
    'ColumnEditable', [true true true], ...
    'Position', [20 110 900 350]);

% buttons
uibutton(f,'Text','Add row','Position',[20 60 90 30], ...
    'ButtonPushedFcn', @(~,~) addRow());
uibutton(f,'Text','Delete row','Position',[120 60 90 30], ...
    'ButtonPushedFcn', @(~,~) deleteRow());
uibutton(f,'Text','Pick MOV folder','Position',[240 60 120 30], ...
    'ButtonPushedFcn', @(~,~) pickFolderForSelectedCol(2,'Pick MOV folder'));
uibutton(f,'Text','Pick REF folder','Position',[370 60 120 30], ...
    'ButtonPushedFcn', @(~,~) pickFolderForSelectedCol(3,'Pick REF folder'));

uibutton(f,'Text','OK','Position',[790 60 50 30], ...
    'ButtonPushedFcn', @(~,~) uiresume(f));
uibutton(f,'Text','Cancel','Position',[850 60 70 30], ...
    'ButtonPushedFcn', @(~,~) cancel());

uiwait(f);

if ~isvalid(f)
    conds = {};
    beadChan = defaultBeadChan;
    strictPosMatching = logical(defaultStrict);
    return;
end

data = tbl.Data;
beadChan = beadField.Value;
strictPosMatching = strictBox.Value;

delete(f);

% validate beadChan
if ~isscalar(beadChan) || isnan(beadChan) || beadChan < 1 || mod(beadChan,1)~=0
    error('Invalid beadChan. Must be a positive integer.');
end
beadChan = double(beadChan);

% clean/validate conditions
if isempty(data)
    conds = {};
    return;
end

isEmptyLabel = cellfun(@(x) strlength(string(x))==0, data(:,1));
data = data(~isEmptyLabel,:);

for i = 1:size(data,1)
    lab = string(data{i,1});
    mov = string(data{i,2});
    ref = string(data{i,3});
    if strlength(mov)==0 || exist(mov,'dir')~=7
        error('Invalid/missing MOV folder for "%s" (row %d).', lab, i);
    end
    if strlength(ref)==0 || exist(ref,'dir')~=7
        error('Invalid/missing REF folder for "%s" (row %d).', lab, i);
    end
end

conds = data;

    function addRow()
        d = tbl.Data;
        d(end+1,:) = {"" ,"" ,""};
        tbl.Data = d;
    end
    function deleteRow()
        d = tbl.Data;
        if isempty(d), return; end
        sel = tbl.Selection;
        if isempty(sel)
            d(end,:) = [];
        else
            rows = unique(sel(:,1));
            d(rows,:) = [];
        end
        tbl.Data = d;
    end
    function pickFolderForSelectedCol(colIdx, titleTxt)
        d = tbl.Data;
        if isempty(d), return; end
        sel = tbl.Selection;
        if isempty(sel)
            row = size(d,1);
        else
            row = sel(1,1);
        end
        folder = uigetdir(pwd, titleTxt);
        if isequal(folder,0), return; end
        d{row,colIdx} = folder;
        tbl.Data = d;
    end
    function cancel()
        delete(f);
        conds = {};
    end
end

%% ===================== Helper: ref matching (strict/lenient) =====================
function refHit = pickRefFileForMovie_singleRefPerPos(movieName, rFiles, strictPosMatching)
if nargin < 3, strictPosMatching = true; end

tok = regexp(movieName, 'Pos\d+', 'match', 'once');

if strictPosMatching
    % Must have Pos token
    assert(~isempty(tok), 'Strict Pos matching ON: movie filename lacks Pos# token: %s', movieName);

    hit = contains({rFiles.name}, tok);
    nHit = nnz(hit);
    assert(nHit > 0, 'Strict Pos matching ON: no ref file matches token "%s" for movie=%s', tok, movieName);
    assert(nHit == 1, 'Strict Pos matching ON: %d ref files match token "%s" for movie=%s (ambiguous)', nHit, tok, movieName);

    refHit = rFiles(find(hit,1,'first'));
    return;
end

% ---- lenient mode (fallbacks) ----
if ~isempty(tok)
    hit = contains({rFiles.name}, tok);
    if any(hit)
        idx = find(hit);
        if numel(idx) > 1
            warning('Multiple ref files match token "%s" for movie=%s. Using first match: %s', ...
                tok, movieName, rFiles(idx(1)).name);
        end
        refHit = rFiles(idx(1));
        return;
    else
        warning('No ref file matches token "%s" for movie=%s. Falling back.', tok, movieName);
    end
end

if numel(rFiles) == 1
    refHit = rFiles(1);
else
    warning('Multiple ref files exist but cannot uniquely match movie=%s. Using first ref file: %s', ...
        movieName, rFiles(1).name);
    refHit = rFiles(1);
end
end

%% deduplicateMDArray(MDarr)
function [MDu, report] = deduplicateMDArray(MDarr)
    n = numel(MDarr);
    keys = strings(n,1);
    details = strings(n,1);
    
    for i = 1:n
        [k, d] = getMDRawSignature(MDarr(i));
        keys(i) = k;
        details(i) = d;
    end
    
    [ukeys, ~, g] = unique(keys, 'stable');
    keepIdx = zeros(numel(ukeys),1);
    groups = cell(0,1);
    
    for ui = 1:numel(ukeys)
        idx = find(g==ui);
        keepIdx(ui) = idx(1);
        if numel(idx) > 1
            grp.keep = idx(1);
            grp.remove = idx(2:end);
            grp.allIdx = idx(:)';
            groups{end+1} = grp; %#ok<AGROW>
        end
    end
    
    keepIdx = sort(keepIdx);
    MDu = MDarr(keepIdx);
    
    report = struct();
    report.nOriginal = n;
    report.nUnique = numel(keepIdx);
    report.nRemoved = n - numel(keepIdx);
    report.keys = keys;
    report.details = details;
    report.groups = groups;
    end
    
    function [key, detail] = getMDRawSignature(MD)
    % Robust ?same movie? signature:
    % prefer channelPath + first filename + bytes + datenum.
    key = "";
    detail = "";
    
    try
        if ~isempty(MD.channels_)
            ch = MD.channels_(1);
    
            if ismethod(ch,'getImageFileNames')
                fns = ch.getImageFileNames();
                if ~isempty(fns)
                    if isprop(ch,'channelPath_') && ~isempty(ch.channelPath_)
                        fullp = fullfile(ch.channelPath_, fns{1});
                    else
                        fullp = fns{1};
                    end
    
                    fullp = char(fullp);
                    d = dir(fullp);
                    if ~isempty(d)
                        stat = sprintf('bytes=%d|datenum=%.0f', d.bytes, d.datenum);
                    else
                        stat = 'nostat';
                    end
    
                    key = string(fullp) + "::" + string(stat);
                    detail = "raw=" + string(fullp) + " | " + string(stat) + " | out=" + string(MD.outputDirectory_);
                    return;
                end
            end
        end
    catch
        % fall through
    end
    
    key = "OUTDIR:" + string(MD.outputDirectory_);
    detail = "fallback_out=" + string(MD.outputDirectory_);
end

function [mFilesOut, logLines] = filterMMChunkedOmeTiffs(mFilesIn)
% Collapses Micro-Manager split OME-TIFF chunk files:
%   MMStack_Pos0.ome.tif, MMStack_Pos0_1.ome.tif, MMStack_Pos0_2.ome.tif ...
% Keeps only one file per logical dataset.
%
% Preference:
%   - keep the base (no _\d+) if present
%   - else keep the smallest suffix

    names = string({mFilesIn.name});
    logLines = {};
    
    % Extract "base key" by removing the optional _<digits> right before .ome.tif
    % Example: MMStack_Pos0_3.ome.tif -> MMStack_Pos0.ome.tif (key)
    keys = regexprep(names, '(_\d+)?\.ome\.tif$', '.ome.tif', 'ignorecase');
    
    [ukeys, ~, g] = unique(keys, 'stable');
    
    keep = false(numel(mFilesIn),1);
    
    for ui = 1:numel(ukeys)
        idx = find(g==ui);
        if numel(idx)==1
            keep(idx) = true;
            continue;
        end
    
        % Prefer the true base file (no _\d+ before .ome.tif)
        isBase = ~contains(names(idx), regexpPattern('_\d+\.ome\.tif$'), 'IgnoreCase', true);
    
        if any(isBase)
            % If multiple base candidates (unlikely), keep the first
            pick = idx(find(isBase,1,'first'));
        else
            % Otherwise pick smallest numeric suffix
            % Parse suffix; missing -> Inf (but here missing won't happen)
            suf = inf(numel(idx),1);
            for k = 1:numel(idx)
                tok = regexp(names(idx(k)), '_([0-9]+)\.ome\.tif$', 'tokens', 'once', 'ignorecase');
                if ~isempty(tok)
                    suf(k) = str2double(tok{1});
                end
            end
            [~,minI] = min(suf);
            pick = idx(minI);
        end
    
        keep(pick) = true;
    
        removed = setdiff(idx, pick);
        logLines{end+1} = sprintf('%s  | kept: %s  | removed: %s', ...
            ukeys(ui), names(pick), strjoin(names(removed), ", "));
    end
    
    mFilesOut = mFilesIn(keep);
end