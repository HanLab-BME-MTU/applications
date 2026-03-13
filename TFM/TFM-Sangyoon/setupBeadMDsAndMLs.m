function setupBeadMDsAndMLs()
% VERSION 4 ? 2026-03-12
fprintf('=== setupBeadMDsAndMLs v4 starting ===\n');
fprintf('Running from: %s\n', mfilename('fullpath'));

% Walks a root directory containing per-bead-image subfolders (each with
% one *.ome.tif), uses bfImport to create a MovieData per file, then
% groups MovieDatas into MovieLists based on shared folder-name prefix
% (everything before the trailing integer, e.g. "BBSG1after" from "BBSG1after3").
%
% Folder naming convention:
%   <GroupName><Index>   e.g. BBSG1after1, BBSG1after2, BBSG2P3 ...
%
% Output:
%   movieData.mat per subfolder (output dir = subfolder itself)
%   movieList.mat per group in rootDir/<GroupName>/
%   allMovieLists.mat in rootDir
%
% Sangyoon Han, 2026.

%% ---- User settings -------------------------------------------------------
rootDir = uigetdir(pwd, 'Select root folder containing bead subfolders');
if rootDir == 0, error('No folder selected.'); end

% Emission wavelength (nm) ? used to fill missing channel metadata and to
% compute psfSigma for display. Adjust to your bead fluorophore.
emissionWL   = 560;    % nm
defaultNA    = 1.45;   % fallback if NA absent from ome metadata
defaultPixNm = 108;    % nm ? fallback if pixel size absent or clearly wrong

%% ---- Discover subfolders -------------------------------------------------
allEntries = dir(rootDir);
allEntries = allEntries([allEntries.isdir]);
allEntries = allEntries(~ismember({allEntries.name}, {'.','..'}));
if isempty(allEntries)
    error('No subfolders found in %s', rootDir)
end
fprintf('Found %d subfolders in %s\n', numel(allEntries), rootDir)

%% ---- Per-subfolder: bfImport + fix metadata ------------------------------
validFolders = {};
validMDPaths = {};

for fi = 1:numel(allEntries)
    folderName = allEntries(fi).name;
    folderPath = fullfile(rootDir, folderName);

    % Find ome.tif (prefer *.ome.tif, fall back to any *.tif)
    tifFiles = dir(fullfile(folderPath, '*.ome.tif'));
    if isempty(tifFiles)
        tifFiles = dir(fullfile(folderPath, '*.tif'));
    end
    if isempty(tifFiles)
        fprintf('  [SKIP] %s ? no .tif found\n', folderName)
        continue
    end

    tifFile = fullfile(folderPath, tifFiles(1).name);
    fprintf('  [%d/%d] %s  ->  %s\n', fi, numel(allEntries), folderName, tifFiles(1).name)

    try
        %% bfImport handles all metadata, Channel/MD construction, sanityCheck.
        MD = bfImport(tifFile, 'outputDirectory', folderPath);
        MD = MD(1);

        %% Fix pixel size if metadata value is wrong (e.g. the 650 nm bug)
        if isempty(MD.pixelSize_) || MD.pixelSize_ > 300
            warning('setupBeadMDs:badPixelSize', ...
                '%s: metadata pixel size = %g nm ? overriding to %g nm', ...
                folderName, MD.pixelSize_, defaultPixNm)
            MD.pixelSize_ = defaultPixNm;
        end

        %% Fix NA if missing
        if isempty(MD.numAperture_) || MD.numAperture_ == 0
            MD.numAperture_ = defaultNA;
        end

        %% Ensure emission wavelength is set on all channels.
        % psfSigma_ is derived by Channel.calculatePSFSigma() from
        % emissionWavelength_, NA, and pixelSize_ ? do not set it directly.
        for iChan = 1:numel(MD.channels_)
            if isempty(MD.channels_(iChan).emissionWavelength_)
                MD.channels_(iChan).emissionWavelength_ = emissionWL;
            end
        end

        %% Save MD.
        % Do NOT call MD.setFilename() ? movieDataFileName_ is locked
        % (read-only via checkProperty) once bfImport has already set it.
        % bfImport names the file after the ome.tif base name; just use
        % whatever path bfImport chose via MD.getFullPath().
        MD.save;
        mdFullPath = MD.getFullPath();

        % Compute sigma for logging (psfSigma_ is derived/read-only)
        psfSigma = 0.21 * emissionWL / MD.numAperture_ / MD.pixelSize_;
        fprintf('    -> MD saved: %s\n', mdFullPath)
        fprintf('       pixSize=%g nm  NA=%.2f  sigma~%.3f px\n', ...
                MD.pixelSize_, MD.numAperture_, psfSigma)

        validFolders{end+1} = folderName;   %#ok<SAGROW>
        validMDPaths{end+1} = mdFullPath;   %#ok<SAGROW>

    catch ME
        fprintf('    [FAILED] %s\n', folderName)
        fprintf('      %s\n', ME.message)
        for si = 1:numel(ME.stack)
            fprintf('        at %s (line %d)\n', ME.stack(si).name, ME.stack(si).line)
        end
    end
end

fprintf('\n%d / %d MovieDatas created successfully.\n', ...
        numel(validFolders), numel(allEntries))
if isempty(validFolders)
    error('No MovieDatas were created ? check Bio-Formats and folder contents.')
end

%% ---- Group by condition (strip embedded series number) ------------------
% Folder names follow <CoatingName><SeriesNum><Condition><RepNum>
% e.g. "BBSG1after3" -> coating="BBSG", series=1, condition="after", rep=3
% We combine all series numbers into one ML per coating+condition:
%   BBSG1P*, BBSG2P*, BBSG3P*, BBSG4P*  -> "BBSGp"
%   BBSG1after*, BBSG2after*, ...        -> "BBSGafter"
%   BSG1P*, BSG2P*, BSG3P*              -> "BSGp"
%   CG1P*, CG2P*, CG3P*                 -> "CGp"
%   CG1after*, CG2after*, CG3after*     -> "CGafter"
%
% Rule: strip the digits that appear between the coating letters and the
% condition suffix, AND strip the trailing replicate number.
% Regex: match letters, then digits (series), then letters+optional digits (condition+rep)
% We keep only the letter-only prefix + letter-only condition.
prefixes = regexprep(validFolders, '^([A-Za-z]+)\d+([A-Za-z]+)\d+$', '$1$2');
uniquePrefixes = unique(prefixes, 'stable');

fprintf('\nFound %d groups:\n', numel(uniquePrefixes))
for gi = 1:numel(uniquePrefixes)
    members = validFolders(strcmp(prefixes, uniquePrefixes{gi}));
    fprintf('  %-20s (%d movies): %s\n', uniquePrefixes{gi}, numel(members), ...
            strjoin(members, ', '))
end

%% ---- Build one MovieList per group --------------------------------------
% Each ML is saved directly into rootDir/MLs_for_beadAnalysis/<GroupPrefix>/
% as 'movieList.mat' ? this IS the canonical save location so no copying
% is needed and paths remain consistent when beadAnalysisBatch loads the ML.
repRootDir = fullfile(rootDir, 'MLs_for_beadAnalysis');
if ~exist(repRootDir, 'dir'), mkdir(repRootDir); end

MLPaths = cell(numel(uniquePrefixes), 1);

for gi = 1:numel(uniquePrefixes)
    grpPrefix  = uniquePrefixes{gi};
    grpMDPaths = validMDPaths(strcmp(prefixes, grpPrefix));

    grpMDs = cellfun(@(p) MovieData.load(p), grpMDPaths, 'UniformOutput', false);

    % Save ML directly into its representative subfolder
    grpMLDir = fullfile(repRootDir, grpPrefix);
    if ~exist(grpMLDir, 'dir'), mkdir(grpMLDir); end

    ML = MovieList(grpMDs, grpMLDir);
    ML.setFilename('movieList.mat');
    ML.setPath(grpMLDir);
    ML.sanityCheck;
    ML.save;

    MLPaths{gi} = fullfile(grpMLDir, 'movieList.mat');
    fprintf('  ML saved: %s  (%d movies)\n', MLPaths{gi}, numel(grpMDs))
end

fprintf('\nAll ML folders ready in:\n  %s\n', repRootDir)

%% ---- Summary ------------------------------------------------------------
save(fullfile(rootDir, 'allMovieLists.mat'), 'MLPaths', 'uniquePrefixes', '-v7.3')
fprintf('\n=== Setup complete ===\n')
fprintf('  %d MovieDatas,  %d MovieLists\n', numel(validFolders), numel(uniquePrefixes))
fprintf('Point beadAnalysisBatch at:\n  %s\n', repRootDir)

%% ---- Note on psfSigma_ --------------------------------------------------
% psfSigma_ is a computed/read-only property of Channel in this u-Track
% version. When estimateBeadDistance calls beadChan.psfSigma_, Channel
% will derive it from emissionWavelength_ and the MovieData's numAperture_
% and pixelSize_. This is why we ensure emissionWavelength_ is set above
% rather than trying to assign psfSigma_ directly.
% If your Channel class does NOT auto-derive psfSigma_, you will need to
% subclass Channel or compute sigma manually in estimateBeadDistance using:
%   sigma = 0.21 * beadChan.emissionWavelength_ / MD.numAperture_ / MD.pixelSize_;

end  % setupBeadMDsAndMLs