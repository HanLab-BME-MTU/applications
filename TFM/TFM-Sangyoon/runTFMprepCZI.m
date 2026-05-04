%% runTFMprepCZI.m

% General TFM preprocessing script.
% Asks for data folder and output folder, then automatically:
%   1) Finds all CZI files
%   2) Separates movies from reference z-stacks using flexible heuristics:
%      - Refs contain '_Z-' or 'Z stack' or 'z-stack' in the filename
%      - Movies are the rest
%   3) Pairs by matching position number (P-02, P02, P04, etc.)
%      and longest common prefix
%   4) Runs make2DTFMFromRefMatched for each pair

%% ---- Ask for paths ----
dataDir = input('Enter data folder path (containing CZI files): ', 's');
dataDir = strtrim(dataDir);
if ~isfolder(dataDir)
    error('Folder not found: %s', dataDir);
end

outRoot = input('Enter output root folder path: ', 's');
outRoot = strtrim(outRoot);
if ~isfolder(outRoot)
    mkdir(outRoot);
    fprintf('Created output folder: %s\n', outRoot);
end

%% ---- Common options ----
beadChan       = 2;
refBeadChan    = 2;
templateFrame  = 1;
refFrame       = 1;
cropROI        = [];
downsample     = 1;
maxShiftFrac   = 0.9;
scoreThreshold = 0.3;
verbose        = true;
writeRef       = true;
overwrite      = true;

%% ---- Find CZI files and classify ----
cziFiles = dir(fullfile(dataDir, '*.czi'));
if isempty(cziFiles)
    error('No CZI files found in: %s', dataDir);
end

% Ref heuristic: filename contains '_Z-' or 'Z stack' or 'z stack' (case-insensitive)
% Anything else is treated as a movie.
refKeywords = {'_Z-', 'Z stack', 'z-stack', 'Zstack'};

movieFiles = {};
refFiles   = {};
for k = 1:numel(cziFiles)
    name = cziFiles(k).name;
    isRef = false;
    for kw = refKeywords
        if contains(name, kw{1}, 'IgnoreCase', true)
            isRef = true;
            break;
        end
    end
    if isRef
        refFiles{end+1}   = name;
    else
        movieFiles{end+1} = name;
    end
end

fprintf('\nFound %d movie file(s) and %d reference z-stack file(s).\n', ...
    numel(movieFiles), numel(refFiles));

if isempty(refFiles)
    error(['No reference z-stack files found.\n' ...
           'Refs must contain one of: _Z-, Z stack, z-stack, Zstack\n' ...
           'Files found:\n  %s'], strjoin({cziFiles.name}, '\n  '));
end

%% ---- Extract position number from filename ----
% Handles formats: _P-02-, _P02_, _P-02_, P-02-, P02-, etc.
function pos = extractPos(name)
    % Handle formats: -P01-, _P-02-, _P02_, P04_, -P04-
    % Step 1: try with word boundary around Pdigits
    tok = regexp(name, '[-_\s]P[-_]?0*(\d+)[-_\s]', 'tokens', 'once');
    if ~isempty(tok), pos = tok{1}; return; end
    % Step 2: P at end of token (e.g. P01-Airyscan)
    tok = regexp(name, '[-_\s]P[-_]?0*(\d+)\b', 'tokens', 'once');
    if ~isempty(tok), pos = tok{1}; return; end
    % Step 3: any P followed by digits
    tok = regexp(name, '\bP0*(\d+)\b', 'tokens', 'once');
    if ~isempty(tok), pos = tok{1}; return; end
    pos = '';
end

%% ---- Auto-pair movies with refs ----
% Strategy:
%   1. Extract position number from each file
%   2. For each movie, find refs with same position number
%   3. Among those, pick the ref with the longest common prefix

pairs = {};

for m = 1:numel(movieFiles)
    mName = movieFiles{m};
    mPos  = extractPos(mName);

    if isempty(mPos)
        warning('Cannot extract position from movie: %s', mName);
        continue;
    end

    % Find refs with same position number (compare as integers to handle leading zeros)
    mPosInt = str2double(mPos);
    candidates = {};
    for r = 1:numel(refFiles)
        rPos = extractPos(refFiles{r});
        if ~isempty(rPos) && str2double(rPos) == mPosInt
            candidates{end+1} = refFiles{r};
        end
    end

    if isempty(candidates)
        warning('No ref found for movie (pos=%s): %s', mPos, mName);
        continue;
    end

    % Pick candidate with longest common prefix with the movie name
    bestRef   = candidates{1};
    bestLen   = commonPrefixLen(mName, candidates{1});
    for c = 2:numel(candidates)
        L = commonPrefixLen(mName, candidates{c});
        if L > bestLen
            bestLen = L;
            bestRef = candidates{c};
        end
    end

    % Build label
    label = buildLabel(mName, mPos);
    pairs{end+1} = {label, mName, bestRef};
end

if isempty(pairs)
    fprintf('\nMovie files:\n');
    for k = 1:numel(movieFiles), fprintf('  %s\n', movieFiles{k}); end
    fprintf('Ref files:\n');
    for k = 1:numel(refFiles),   fprintf('  %s\n', refFiles{k});   end
    error('No movie-reference pairs could be matched.');
end

%% ---- Show pairs and confirm ----
fprintf('\n======== Matched pairs ========\n');
for i = 1:numel(pairs)
    fprintf('[%d] %s\n',       i, pairs{i}{1});
    fprintf('    movie : %s\n',   pairs{i}{2});
    fprintf('    ref   : %s\n\n', pairs{i}{3});
end
fprintf('================================\n');

confirm = input('Proceed with these pairs? (y/n): ', 's');
if ~strcmpi(strtrim(confirm), 'y')
    disp('Cancelled.');
    return;
end

%% ---- Main loop ----
mdAll = {};

for i = 1:numel(pairs)
    label     = pairs{i}{1};
    movieFile = fullfile(dataDir, pairs{i}{2});
    refFile   = fullfile(dataDir, pairs{i}{3});
    outDir    = fullfile(outRoot, label);

    fprintf('\n========================================\n');
    fprintf('[%d/%d] %s\n', i, numel(pairs), label);
    fprintf('  movie : %s\n', pairs{i}{2});
    fprintf('  ref   : %s\n', pairs{i}{3});
    fprintf('========================================\n');

    try
        mdList = make2DTFMFromRefMatched(movieFile, refFile, outDir, ...
            'beadChan',       beadChan,       ...
            'refBeadChan',    refBeadChan,    ...
            'templateFrame',  templateFrame,  ...
            'refFrame',       refFrame,       ...
            'cropROI',        cropROI,        ...
            'downsample',     downsample,     ...
            'maxShiftFrac',   maxShiftFrac,   ...
            'scoreThreshold', scoreThreshold, ...
            'verbose',        verbose,        ...
            'writeRefSlice',  writeRef,       ...
            'overwrite',      overwrite);

        mdAll{end+1} = struct('label', label, 'mdList', mdList);
        fprintf('[%s] Done. %d MovieData object(s) created.\n', label, numel(mdList));

    catch ME
        warning('[%s] FAILED: %s\n%s', label, ME.message, getReport(ME));
    end
end

%% ---- Summary ----
fprintf('\n\n======== Summary ========\n');
for i = 1:numel(mdAll)
    lb  = mdAll{i}.label;
    mds = mdAll{i}.mdList;
    fprintf('%s : %d MovieData\n', lb, numel(mds));
    for j = 1:numel(mds)
        fprintf('   [%d] %s\n', j, mds(j).getFullPath());
    end
end
fprintf('=========================\n');

%% ---- Local helpers ----
function n = commonPrefixLen(a, b)
    n = 0;
    for k = 1:min(numel(a), numel(b))
        if a(k) == b(k), n = n + 1; else, break; end
    end
end

function label = buildLabel(movieName, posNum)
    % Remove extension and sanitise
    [~, base, ~] = fileparts(movieName);
    % Strip everything from position token onwards
    base = regexprep(base, '[_\s]P[-_]?' + string(posNum) + '.*', '', 'ignorecase');
    base = regexprep(base, '[^\w]+', '_');
    base = regexprep(base, '^_|_$', '');
    label = [base '_P' posNum];
end