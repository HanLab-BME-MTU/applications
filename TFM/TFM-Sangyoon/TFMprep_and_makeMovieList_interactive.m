%% TFMprep_and_makeMovieList_interactive.m
%
% Improvements over previous version:
%   1. GUI dialogs removed - pure CLI input only
%   2. Full BioFormats metadata extraction into MD/Channel:
%      - emissionWavelength_, excitationWavelength_, name_ (from setupMovieDataFromND pattern)
%      - numericalAperture_, camBitdepth_ (from BioFormats store)
%      - timeInterval_ (interactive prompt if missing, as in setupMovieDataFromND)
%   3. Multi-format support: .ome.tif, .nd, .nd2, .czi
%
% Sangyoon Han

clear; clc;

%% ==== USER SETTINGS ====
baseOutName        = "TFMprep_bestZ";
metric             = "tenengrad";
verbose            = true;
movieListsRoot     = "";   % "" = inside each movDir/baseOutName/label/MovieList
supportedExts      = {'.ome.tif', '.nd', '.nd2', '.czi'};

%% ==== CLI INPUT: conditions + beadChan + strictPosMatching ====
[conds, beadChan, strictPosMatching] = getInputsCLI(2, true);
assert(~isempty(conds) && size(conds,2)==3, 'No conditions provided.');
fprintf('beadChan = %d | strictPosMatching = %s\n', beadChan, string(strictPosMatching));

%% ==== MAIN LOOP ====
for i = 1:size(conds,1)
    label  = string(conds{i,1});
    movDir = string(conds{i,2});
    refDir = string(conds{i,3});

    fprintf('\n====================================================\n');
    fprintf('Condition %d/%d: %s\n', i, size(conds,1), label);
    fprintf('  movDir: %s\n  refDir: %s\n', movDir, refDir);
    fprintf('====================================================\n');

    assert(exist(movDir,'dir')==7, 'movDir not found: %s', movDir);
    assert(exist(refDir,'dir')==7, 'refDir not found: %s', refDir);

    % ---- find movie files (all supported formats) ----
    mFiles = findMovieFiles(movDir, supportedExts);
    [mFiles, chunkLog] = filterMMChunkedOmeTiffs(mFiles);
    if ~isempty(chunkLog)
        fprintf('[ChunkFilter] Collapsed split movies:\n');
        fprintf('  %s\n', chunkLog{:});
    end
    mFiles = mFiles(~contains({mFiles.name}, "metadata", 'IgnoreCase', true));
    assert(~isempty(mFiles), 'No movie files found in %s', movDir);

    % ---- find ref files ----
    rFiles = findMovieFiles(refDir, supportedExts);
    rFiles = rFiles(~contains({rFiles.name}, "metadata", 'IgnoreCase', true));
    assert(~isempty(rFiles), 'No ref files found in %s', refDir);

    allMD  = MovieData.empty;
    mapLog = strings(0,1);

    % ---- prompt time interval once per condition (as in setupMovieDataFromND) ----
    timeLapseAsked = false;
    timeLapse = [];

    for k = 1:numel(mFiles)
        movieFile = fullfile(mFiles(k).folder, mFiles(k).name);
        refHit    = pickRefFileForMovie(mFiles(k).name, rFiles, strictPosMatching);
        refFile   = fullfile(refHit.folder, refHit.name);

        [~, baseName] = fileparts(mFiles(k).name);
        baseName = regexprep(baseName, '\.ome$', '');   % strip .ome if present
        outDir   = fullfile(movDir, baseOutName, label, baseName);

        fprintf('\n=== %s | movie=%s ===\n', label, mFiles(k).name);
        fprintf('    ref = %s\n    out = %s\n', string(refFile), string(outDir));

        % ---- extract MovieData via make2DTFMFromRefBestZ_seriesMatched ----
        mdList = make2DTFMFromRefBestZ_seriesMatched(movieFile, refFile, outDir, ...
            'beadChan', beadChan, ...
            'metric',   metric, ...
            'verbose',  verbose);

        if isempty(mdList)
            warning('No MovieData created for %s', mFiles(k).name);
            continue;
        end

        % ---- enrich each MD with full BioFormats metadata ----
        for mi = 1:numel(mdList)
            md = mdList(mi);

            % time interval: ask once if missing and nFrames>1
            if ~timeLapseAsked && (isempty(md.timeInterval_) || md.timeInterval_<=0) && md.nFrames_>1
                timeLapse = input('What is the time interval in seconds? ');
                if isempty(timeLapse), timeLapse = 1; end
                timeLapseAsked = true;
            end
            if ~isempty(timeLapse)
                md.timeInterval_ = timeLapse;
            end

            % per-channel metadata from BioFormats store
            md = enrichMDFromBioFormats(md, movieFile, mi-1, verbose);

            md.sanityCheck();
            md.save();
            mdList(mi) = md;
        end

        allMD = [allMD; mdList(:)]; %#ok<AGROW>
        mapLog(end+1,1) = sprintf("movie=%s | ref=%s | out=%s | nMD=%d", ...
            string(mFiles(k).name), string(refHit.name), outDir, numel(mdList));
    end

    fprintf('\n[%s] Total MD: %d\n', label, numel(allMD));
    assert(~isempty(allMD), 'No MovieData created for condition %s', label);

    % ---- MovieList output dir ----
    if strlength(movieListsRoot) == 0
        mlDir = fullfile(movDir, baseOutName, label, "MovieList");
    else
        mlDir = fullfile(movieListsRoot, label);
    end
    if ~exist(mlDir,'dir'); mkdir(mlDir); end
    mlDir = char(mlDir);

    % ---- deduplicate ----
    [allMD, dedupReport] = deduplicateMDArray(allMD);
    fprintf('[%s] Dedup: %d -> %d unique (removed %d)\n', ...
        label, dedupReport.nOriginal, dedupReport.nUnique, dedupReport.nRemoved);

    saveDedupReport(mlDir, label, dedupReport);

    ML = MovieList(allMD, mlDir);
    ML.setFilename('movieList.mat');
    ML.sanityCheck;
    ML.save();
    fprintf('[%s] Saved MovieList: %s\n', label, fullfile(mlDir,'movieList.mat'));

    saveMapLog(mlDir, label, movDir, refDir, beadChan, strictPosMatching, mapLog);
end

disp('Done.');

%% =====================================================================
%% LOCAL FUNCTIONS
%% =====================================================================

% ------------------------------------------------------------------
% 1. CLI input (bead channel + strictPosMatching via CLI,
%    folder selection via uigetdir dialog)
% ------------------------------------------------------------------
function [conds, beadChan, strictPosMatching] = getInputsCLI(defaultBeadChan, defaultStrict)
if nargin < 1, defaultBeadChan = 2; end
if nargin < 2, defaultStrict = true; end

fprintf('\n======= TFMprep Interactive Setup =======\n');

beadChan = input(sprintf('Bead channel index [default %d]: ', defaultBeadChan));
if isempty(beadChan), beadChan = defaultBeadChan; end
assert(isscalar(beadChan) && beadChan>=1 && mod(beadChan,1)==0, 'Invalid beadChan.');

ans_ = input(sprintf('Strict Pos matching? 1=yes 0=no [default %d]: ', defaultStrict));
if isempty(ans_), strictPosMatching = logical(defaultStrict);
else,             strictPosMatching = logical(ans_); end

n = input('Number of conditions: ');
assert(isscalar(n) && n>=1, 'Invalid number.');

conds = cell(n, 3);
for i = 1:n
    fprintf('\n--- Condition %d/%d ---\n', i, n);
    label = strtrim(input(sprintf('Label (e.g. WT_40nm): '), 's'));
    assert(~isempty(label), 'Label cannot be empty.');

    movDir = uigetdir(pwd, sprintf('Select MOVIE folder for "%s"', label));
    assert(~isequal(movDir, 0), 'Movie folder selection cancelled.');

    refDir_ = uigetdir(movDir, sprintf('Select REF folder for "%s"', label));
    assert(~isequal(refDir_, 0), 'Ref folder selection cancelled.');

    conds{i,1} = label;
    conds{i,2} = movDir;
    conds{i,3} = refDir_;
end
end

% ------------------------------------------------------------------
% 2. Find movie files across supported formats
% ------------------------------------------------------------------
function files = findMovieFiles(folder, exts)
files = struct('name',{},'folder',{},'bytes',{},'datenum',{});
for ei = 1:numel(exts)
    ext = exts{ei};
    if strcmp(ext, '.ome.tif')
        pat = '*.ome.tif';
    else
        pat = ['*' ext];
    end
    d = dir(fullfile(folder, pat));
    for di = 1:numel(d)
        files(end+1) = struct('name', d(di).name, 'folder', d(di).folder, ...
            'bytes', d(di).bytes, 'datenum', d(di).datenum); %#ok<AGROW>
    end
end
% Remove duplicates by name
if ~isempty(files)
    [~, ui] = unique({files.name}, 'stable');
    files = files(ui);
end
end

% ------------------------------------------------------------------
% 3. Enrich MD channels with BioFormats metadata
%    (emission, excitation wavelength, NA, bit depth, channel name)
%    Follows setupMovieDataFromND pattern: name2wavelength for emission,
%    plus direct OME store queries for NA and camBitdepth_.
% ------------------------------------------------------------------
function md = enrichMDFromBioFormats(md, movieFile, seriesIdx, verbose)
if nargin < 4, verbose = false; end

try
    reader = bfGetReader(char(movieFile));
    reader.setSeries(seriesIdx);
    store = reader.getMetadataStore();
    nChan = reader.getSizeC();

    for c0 = 0:nChan-1
        ci = c0 + 1;
        if ci > numel(md.channels_), break; end
        ch = md.channels_(ci);

        % ---- channel name ----
        try
            cname = char(store.getChannelName(seriesIdx, c0));
            if ~isempty(cname) && isprop(ch,'name_')
                ch.name_ = cname;
            end
        catch; end

        % ---- emission wavelength (nm) ----
        try
            w = store.getChannelEmissionWavelength(seriesIdx, c0);
            if ~isempty(w) && isprop(ch,'emissionWavelength_')
                ch.emissionWavelength_ = w.value().doubleValue();
            elseif isprop(ch,'emissionWavelength_') && ~isempty(ch.name_)
                % fallback: name2wavelength (from setupMovieDataFromND)
                wNm = name2wavelength(ch.name_);
                if ~isnan(wNm)
                    ch.emissionWavelength_ = wNm * 1e9;
                end
            end
        catch; end

        % ---- excitation wavelength (nm) ----
        try
            wx = store.getChannelExcitationWavelength(seriesIdx, c0);
            if ~isempty(wx) && isprop(ch,'excitationWavelength_')
                ch.excitationWavelength_ = wx.value().doubleValue();
            end
        catch; end

        % ---- camera bit depth ----
        try
            bd = store.getPixelsSignificantBits(seriesIdx);
            if ~isempty(bd) && isprop(ch,'camBitdepth_')
                ch.camBitdepth_ = double(bd.getValue());
            end
        catch; end

        % ---- numerical aperture ----
        try
            % OME: ObjectiveSettings -> Objective NA
            objSettingsID = store.getObjectiveSettingsID(seriesIdx);
            if ~isempty(objSettingsID)
                % find the objective in the list
                nInst = store.getInstrumentCount();
                for inst = 0:nInst-1
                    nObj = store.getObjectiveCount(inst);
                    for oi = 0:nObj-1
                        objID = store.getObjectiveID(inst, oi);
                        if strcmp(char(objID), char(objSettingsID))
                            na = store.getObjectiveLensNA(inst, oi);
                            if ~isempty(na) && isprop(ch,'numericalAperture_')
                                ch.numericalAperture_ = na.doubleValue();
                            end
                        end
                    end
                end
            end
        catch; end

        if verbose
            fprintf('  Ch%d: name=%s | em=%.0f nm | ex=%.0f nm | NA=%.2f | bits=%d\n', ci, ...
                getpropOrDash(ch,'name_'), ...
                getpropOrDash(ch,'emissionWavelength_'), ...
                getpropOrDash(ch,'excitationWavelength_'), ...
                getpropOrDash(ch,'numericalAperture_'), ...
                getpropOrDash(ch,'camBitdepth_'));
        end

        md.channels_(ci) = ch;
    end

    try; reader.close(); catch; end
catch ME
    if verbose, warning(ME.identifier, 'enrichMDFromBioFormats: %s', ME.message); end
end
end

function val = getpropOrDash(obj, prop)
try
    v = obj.(prop);
    if isempty(v) || (isnumeric(v) && isnan(v)), val = NaN;
    else, val = v; end
catch
    val = NaN;
end
end

% ------------------------------------------------------------------
% 4. Ref file matching
% ------------------------------------------------------------------
function refHit = pickRefFileForMovie(movieName, rFiles, strictPosMatching)
if nargin < 3, strictPosMatching = true; end

tok = regexp(movieName, 'Pos\d+', 'match', 'once');

if strictPosMatching
    assert(~isempty(tok), 'Strict Pos matching: movie has no Pos# token: %s', movieName);
    hit = contains({rFiles.name}, tok);
    assert(nnz(hit)==1, 'Strict Pos matching: %d ref files match "%s" for %s', nnz(hit), tok, movieName);
    refHit = rFiles(find(hit,1));
    return;
end

if ~isempty(tok)
    hit = contains({rFiles.name}, tok);
    if any(hit)
        idx = find(hit);
        if numel(idx)>1
            warning('Multiple refs match "%s" for %s. Using first.', tok, movieName);
        end
        refHit = rFiles(idx(1)); return;
    end
end

refHit = rFiles(1);
warning('Could not uniquely match ref for %s. Using: %s', movieName, rFiles(1).name);
end

% ------------------------------------------------------------------
% 5. Deduplication (unchanged from original)
% ------------------------------------------------------------------
function [MDu, report] = deduplicateMDArray(MDarr)
n = numel(MDarr);
keys = strings(n,1); details = strings(n,1);
for i = 1:n
    [keys(i), details(i)] = getMDRawSignature(MDarr(i));
end
[ukeys, ~, g] = unique(keys, 'stable');
keepIdx = zeros(numel(ukeys),1);
groups = cell(0,1);
for ui = 1:numel(ukeys)
    idx = find(g==ui);
    keepIdx(ui) = idx(1);
    if numel(idx)>1
        grp.keep=idx(1); grp.remove=idx(2:end); grp.allIdx=idx(:)';
        groups{end+1}=grp; %#ok<AGROW>
    end
end
MDu = MDarr(sort(keepIdx));
report = struct('nOriginal',n,'nUnique',numel(keepIdx),'nRemoved',n-numel(keepIdx),...
    'keys',keys,'details',details,'groups',{groups});
end

function [key, detail] = getMDRawSignature(MD)
key = ""; detail = "";
try
    if ~isempty(MD.channels_)
        ch = MD.channels_(1);
        if ismethod(ch,'getImageFileNames')
            fns = ch.getImageFileNames();
            if ~isempty(fns)
                fullp = char(fullfile(ch.channelPath_, fns{1}));
                d = dir(fullp);
                stat = '';
                if ~isempty(d), stat = sprintf('bytes=%d|datenum=%.0f',d.bytes,d.datenum); end
                key = string(fullp) + "::" + string(stat);
                detail = "raw=" + string(fullp) + " | " + stat + " | out=" + MD.outputDirectory_;
                return;
            end
        end
    end
catch; end
key = "OUTDIR:" + string(MD.outputDirectory_);
detail = "fallback_out=" + string(MD.outputDirectory_);
end

% ------------------------------------------------------------------
% 6. Chunk filter (unchanged)
% ------------------------------------------------------------------
function [mFilesOut, logLines] = filterMMChunkedOmeTiffs(mFilesIn)
names = string({mFilesIn.name});
logLines = {};
keys = regexprep(names, '(_\d+)?\.ome\.tif$', '.ome.tif', 'ignorecase');
[ukeys, ~, g] = unique(keys, 'stable');
keep = false(numel(mFilesIn),1);
for ui = 1:numel(ukeys)
    idx = find(g==ui);
    if numel(idx)==1, keep(idx)=true; continue; end
    isBase = ~contains(names(idx), regexpPattern('_\d+\.ome\.tif$'), 'IgnoreCase', true);
    if any(isBase)
        pick = idx(find(isBase,1));
    else
        suf = inf(numel(idx),1);
        for k=1:numel(idx)
            tok=regexp(names(idx(k)),'_([0-9]+)\.ome\.tif$','tokens','once','ignorecase');
            if ~isempty(tok), suf(k)=str2double(tok{1}); end
        end
        [~,minI]=min(suf); pick=idx(minI);
    end
    keep(pick)=true;
    removed=setdiff(idx,pick);
    logLines{end+1}=sprintf('%s | kept: %s | removed: %s', ukeys(ui), names(pick), strjoin(names(removed),", ")); %#ok<AGROW>
end
mFilesOut = mFilesIn(keep);
end

% ------------------------------------------------------------------
% 7. Report saving helpers
% ------------------------------------------------------------------
function saveDedupReport(mlDir, label, r)
fid = fopen(fullfile(mlDir,'TFMprep_dedup_report.txt'),'w');
fprintf(fid,"Condition: %s\nOriginal: %d\nUnique: %d\nRemoved: %d\n\n", ...
    label, r.nOriginal, r.nUnique, r.nRemoved);
for gi=1:numel(r.groups)
    g=r.groups{gi};
    fprintf(fid,"Group %d: keep=%d remove=%s\n  key: %s\n",gi,g.keep,mat2str(g.remove),r.keys(g.keep));
    for jj=1:numel(g.allIdx)
        fprintf(fid,"    %d) %s\n",g.allIdx(jj),r.details(g.allIdx(jj)));
    end
    fprintf(fid,"\n");
end
fclose(fid);
end

function saveMapLog(mlDir, label, movDir, refDir, beadChan, strict, mapLog)
fid = fopen(fullfile(mlDir,'TFMprep_mapLog.txt'),'w');
fprintf(fid,"Condition: %s\nmovDir: %s\nrefDir: %s\nbeadChan: %d\nstrictPos: %d\n\n",...
    label, movDir, refDir, beadChan, strict);
for j=1:numel(mapLog), fprintf(fid,"%s\n",mapLog(j)); end
fclose(fid);
end