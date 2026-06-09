function setupTFMPackageForMovieList(movieListInput, tfmParams)
% setupTFMPackageForMovieList
%
% Sets up TFMPackage for every MD in a MovieList.
%
% KEY BEHAVIORS:
%   1. NO metadata overwrite: timeInterval_, pixelSize_, channel properties
%      already in MD (set by bfImport/TFMprep) are preserved as-is.
%   2. beadChan: resolves from .beadChan or .iBeadChan (either field name
%      accepted). Used for proc1 shift estimation AND proc2 displacement
%      field tracking.
%   3. referenceFramePath is OPTIONAL. Only set when the file is actually
%      found on disk. If not found, the field is left at its default to
%      avoid the "Image#1" display bug in the u-track GUI.
%
% INPUTS:
%   movieListPath : full path to movieList.mat
%   tfmParams     : struct with any of the fields in defaultTFMParams().
%                   Accepts both .beadChan and .iBeadChan for bead channel.

% Accept either:
%   1) movieList.mat path
%   2) already-loaded MovieList object
if isa(movieListInput, 'MovieList')
    ML = movieListInput;
    fprintf('[setupTFMPackage] Using pre-loaded MovieList object.\n');
else
    movieListPath = char(movieListInput);
    assert(exist(movieListPath,'file')==2, ...
        'movieListPath not found: %s', movieListPath);
    fprintf('[setupTFMPackage] Loading MovieList from disk...\n');
    ML = MovieList.load(movieListPath);
end

nMovies = numel(ML.movieDataFile_);

% ---- Merge user params over defaults ----
p = defaultTFMParams();
if nargin >= 2 && ~isempty(tfmParams)
    fn = fieldnames(tfmParams);
    for i = 1:numel(fn)
        p.(fn{i}) = tfmParams.(fn{i});
    end
end

% ---- Resolve beadChan: accept both .beadChan and .iBeadChan ----
% run_TFMpipeline_interactive uses iBeadChan; this ensures both work.
if isfield(p,'iBeadChan') && (~isfield(p,'beadChan') || p.beadChan == defaultTFMParams().beadChan)
    p.beadChan = p.iBeadChan;
end
fprintf('[setupTFMPackage] beadChan = %d\n', p.beadChan);

for i = 1:nMovies
    MD = ML.getMovie(i);
    fprintf('[setupTFMPackage] %d/%d: %s\n', i, nMovies, MD.outputDirectory_);

    % ---- NO metadata overwrite ----
    % timeInterval_, pixelSize_, channel emissionWavelength_ etc. are
    % already populated by TFMprep / bfImport. Do not touch them here.

    % ---- Resolve reference image path (OPTIONAL) ----
    % Returns '' if not found. We only set proc params when non-empty.
    % Leaving referenceFramePath unset avoids the u-track "Image#1" bug.
    refImgPath = resolveRefBestZPath(MD.outputDirectory_, false);
    if ~isempty(refImgPath)
        fprintf('  ref image: %s\n', refImgPath);
    else
        fprintf('  ref image: not found (referenceFramePath will not be set)\n');
    end

    % ---- Get or create TFMPackage ----
    iPack = MD.getPackageIndex('TFMPackage');
    if isempty(iPack)
        MD.addPackage(TFMPackage(MD));
        iPack = MD.getPackageIndex('TFMPackage');
    end
    tfmPack = MD.getPackage(iPack);

    % =================================================================
    % Process 1: EfficientSubpixelRegistration (stage drift correction)
    % ChannelIndex  = ALL channels (shift is applied to every channel)
    % iBeadChannel  = p.beadChan  (shift is ESTIMATED from bead channel)
    % =================================================================
    if isempty(tfmPack.processes_{1})
        tfmPack.createDefaultProcess(1);
    end
    proc1 = tfmPack.getProcess(1);
    fp1   = proc1.funParams_;

    % All channels get the shift applied
    if isfield(fp1,'ChannelIndex')
        fp1.ChannelIndex = 1:numel(MD.channels_);
    end

    % Bead channel used to ESTIMATE the shift
    fp1 = setField(fp1, {'iBeadChannel','BeadChannelIndex','beadChannel','refChannel'}, p.beadChan);

    % Reference image: only set if the file actually exists
    if ~isempty(refImgPath)
        fp1 = setField(fp1, {'referenceFramePath','referenceFrame'}, refImgPath);
    end

    if isfield(fp1,'usfac'), fp1.usfac = p.usfac; end

    proc1.setPara(fp1);

    % =================================================================
    % Process 2: Displacement field tracking
    % ChannelIndex  = p.beadChan ONLY (tracking runs on bead channel)
    % =================================================================
    if isempty(tfmPack.processes_{2})
        tfmPack.createDefaultProcess(2);
    end
    proc2 = tfmPack.getProcess(2);
    fp2   = proc2.funParams_;

    % Displacement tracking: bead channel only
    if isfield(fp2,'ChannelIndex')
        fp2.ChannelIndex = p.beadChan;
    end

    % Reference image: only set if the file actually exists
    if ~isempty(refImgPath)
        fp2 = setField(fp2, {'referenceFramePath','referenceFrame'}, refImgPath);
    end

    if isfield(fp2,'alpha'),                  fp2.alpha                  = p.alpha;                  end
    if isfield(fp2,'minCorLength'),           fp2.minCorLength           = p.minCorLength;           end
    if isfield(fp2,'addNonLocMaxBeads'),      fp2.addNonLocMaxBeads      = p.addNonLocMaxBeads;      end
    if isfield(fp2,'maxFlowSpeed'),           fp2.maxFlowSpeed           = p.maxFlowSpeed;           end
    if isfield(fp2,'highRes'),                fp2.highRes                = p.highRes;                end
    if isfield(fp2,'useGrid'),                fp2.useGrid                = p.useGrid;                end
    if isfield(fp2,'mode'),                   fp2.mode                   = p.mode;                  end
    if isfield(fp2,'noFlowOutwardOnBorder'),  fp2.noFlowOutwardOnBorder  = p.noFlowOutwardOnBorder; end
    if isfield(fp2,'trackSuccessively'),      fp2.trackSuccessively      = p.trackSuccessively;     end

    proc2.setPara(fp2);

    % =================================================================
    % Process 3: Displacement post-processing (delete if exists)
    % =================================================================
    if numel(tfmPack.processes_) >= 3 && ~isempty(tfmPack.processes_{3})
        MD.deleteProcess(tfmPack.processes_{3});
        MD.save();
    end

    % =================================================================
    % Process 4: Force reconstruction (FTTC/BEM)
    % =================================================================
    if numel(tfmPack.processes_) >= 4
        if isempty(tfmPack.processes_{4})
            tfmPack.createDefaultProcess(4);
        end
        proc4 = tfmPack.getProcess(4);
        fp4   = proc4.funParams_;

        if isfield(fp4,'YoungModulus'),   fp4.YoungModulus   = p.YoungModulus;   end
        if isfield(fp4,'regParam'),       fp4.regParam       = p.regParam;       end
        if isfield(fp4,'method'),         fp4.method         = p.method;         end
        if isfield(fp4,'solMethodBEM'),   fp4.solMethodBEM   = p.solMethodBEM;   end
        if isfield(fp4,'useLcurve'),      fp4.useLcurve      = p.useLcurve;      end

        if isfield(fp4,'basisClassTblPath') && ~isempty(p.basisClassTblPath)
            fp4.basisClassTblPath = p.basisClassTblPath;
        end

        proc4.setPara(fp4);
    end

    % =================================================================
    % Process 5: Strain energy (if present)
    % =================================================================
    if numel(tfmPack.processes_) >= 5 && ~isempty(tfmPack.processes_{5})
        proc5 = tfmPack.getProcess(5);
        fp5   = proc5.funParams_;
        proc5.setPara(fp5);
    end

    MD.save();
    fprintf('[setupTFMPackage] Done: %s\n', MD.outputDirectory_);
end
end

%% ===================== Helpers =====================

function fp = setField(fp, fieldNames, value)
% Set the first field in fieldNames that already exists in fp.
% Never creates a new field - only updates existing ones.
for k = 1:numel(fieldNames)
    if isfield(fp, fieldNames{k})
        fp.(fieldNames{k}) = value;
        return;
    end
end
end

function refImgPath = resolveRefBestZPath(mdOutDir, strict)
% Resolve the reference bestZ image under <mdOutDir>/reference/.
%
% Returns '' when not found and strict=false.
% Returns error when not found and strict=true.
%
% Candidate priority:
%   1) ref_beads_matched.tif   (make2DTFMFromRefMatched output)
%   2) ref_beads_bestZ.tif     (make2DTFMFromRefBestZ output)
%   3) ref_beads_bestZ_series*.tif
%   4) *bestZ*.tif or *matched*.tif
%   5) *ref*.tif

if nargin < 2, strict = false; end

refImgPath = '';  % default: not found

refDir = fullfile(mdOutDir, 'reference');
if exist(refDir,'dir') ~= 7
    msg = sprintf('Reference folder not found: %s', refDir);
    if strict, error(msg); %#ok<SPERR>
    else, warning('setupTFMPackage:noRefDir', '%s', msg); return; end
end

% Priority 1 & 2: exact names
for name = {'ref_beads_matched.tif', 'ref_beads_bestZ.tif'}
    p = fullfile(refDir, name{1});
    if exist(p,'file') == 2
        refImgPath = p;
        return;
    end
end

% Priority 3-5: wildcard patterns
patterns = {'ref_beads_bestZ_series*.tif', '*bestZ*.tif', '*matched*.tif', '*ref*.tif'};
for pi = 1:numel(patterns)
    d = dir(fullfile(refDir, patterns{pi}));
    d = d(~[d.isdir]);
    if isempty(d), continue; end

    if numel(d) == 1
        refImgPath = fullfile(d(1).folder, d(1).name);
        return;
    end

    % Multiple candidates: pick by name score
    names  = string({d.name});
    score  = contains(lower(names),'bestz')   * 3 ...
           + contains(lower(names),'matched') * 3 ...
           + contains(lower(names),'beads')   * 2;
    [~,idx] = max(score);
    refImgPath = fullfile(d(idx).folder, d(idx).name);
    if numel(d) > 1
        warning('setupTFMPackage:multipleRefs', ...
            'Multiple ref candidates in %s. Using: %s', refDir, d(idx).name);
    end
    return;
end

% Not found
msg = sprintf('No reference image found in: %s (referenceFramePath will not be set)', refDir);
if strict, error(msg); %#ok<SPERR>
else, warning('setupTFMPackage:noRefImg', '%s', msg); end
end

function p = defaultTFMParams()
p = struct();
p.beadChan             = 2;
p.usfac                = 20;
p.alpha                = 0.05;
p.minCorLength         = 17;
p.addNonLocMaxBeads    = false;   % default off
p.maxFlowSpeed         = 20;
p.highRes              = true;
p.useGrid              = true;
p.mode                 = 'accurate';
p.noFlowOutwardOnBorder= 1;
p.trackSuccessively    = false;
p.outlierThreshold     = 2;
p.fillVectors          = false;
p.YoungModulus         = 2700;
p.regParam             = 0.005;
p.method               = 'FTTC';
p.solMethodBEM         = 'QR';
p.useLcurve            = false;
p.basisClassTblPath    = '';
end