function setupTFMPackageForMovieList(movieListPath, tfmParams)
% setupTFMPackageForMovieList
% MovieList ? ?? MD? TFMPackage? ??.
%
% movieListPath: ".../movieList.mat"
% tfmParams: struct (??)

movieListPath = char(movieListPath);
assert(exist(movieListPath,'file')==2, 'movieListPath not found: %s', movieListPath);

ML = MovieList.load(movieListPath);
nMovies = numel(ML.movieDataFile_);

% ?? ????
p = defaultTFMParams();
if nargin >= 2 && ~isempty(tfmParams)
    fn = fieldnames(tfmParams);
    for i=1:numel(fn), p.(fn{i}) = tfmParams.(fn{i}); end
end

for i = 1:nMovies
    MD = ML.getMovie(i);

    % ---- Ensure frame rate exists (required by TFMPackage sanityCheck) ----
    if isempty(MD.timeInterval_) || ~(MD.timeInterval_ > 0)
        MD.timeInterval_ = p.timeInterval;   % seconds per frame
    end
    if isempty(MD.pixelSize_) || ~(MD.pixelSize_ > 0)
        MD.pixelSize_ = p.pixelSize;         % nm/pixel
    end
    MD.save;

    % ---- Flexible reference bestZ path resolution ----
    refImgPath = resolveRefBestZPath(MD.outputDirectory_, true); % strict=true

    % --- TFMPackage ?? ---
    iPack = MD.getPackageIndex('TFMPackage');
    if isempty(iPack)
        MD.addPackage(TFMPackage(MD));
        iPack = MD.getPackageIndex('TFMPackage');
    end
    tfmPack = MD.getPackage(iPack);

    % --- Process 1: EfficientSubpixelRegistration ---
    if isempty(tfmPack.processes_{1})
        tfmPack.createDefaultProcess(1);
    end
    proc1 = tfmPack.getProcess(1);
    fp = proc1.funParams_;

    % Apply shift to ALL channels
    if isfield(fp,'ChannelIndex')
        fp.ChannelIndex = 1:numel(MD.channels_);
    end

    % Use bead channel to compute shift
    if isfield(fp,'iBeadChannel')
        fp.iBeadChannel = p.beadChan;
    elseif isfield(fp,'BeadChannelIndex')
        fp.BeadChannelIndex = p.beadChan;
    elseif isfield(fp,'refChannel')
        fp.refChannel = p.beadChan;
    end

    % Reference image (bestZ slice)
    if isfield(fp,'referenceFramePath')
        fp.referenceFramePath = refImgPath;
    elseif isfield(fp,'referenceFrame')
        fp.referenceFrame = refImgPath;
    end

    if isfield(fp,'usfac')
        fp.usfac = p.usfac;
    end

    proc1.setPara(fp);

    % --- Process 2: Displacement field ---
    if isempty(tfmPack.processes_{2}), tfmPack.createDefaultProcess(2); end
    proc2 = tfmPack.getProcess(2);
    fp = proc2.funParams_;

    if isfield(fp,'referenceFramePath'); fp.referenceFramePath = refImgPath; end
    if isfield(fp,'referenceFrame');     fp.referenceFrame     = refImgPath; end
    if isfield(fp,'ChannelIndex');       fp.ChannelIndex       = p.beadChan; end

    if isfield(fp,'alpha'); fp.alpha = p.alpha; end
    if isfield(fp,'minCorLength'); fp.minCorLength = p.minCorLength; end
    if isfield(fp,'addNonLocMaxBeads'); fp.addNonLocMaxBeads = p.addNonLocMaxBeads; end
    if isfield(fp,'maxFlowSpeed'); fp.maxFlowSpeed = p.maxFlowSpeed; end
    if isfield(fp,'highRes'); fp.highRes = p.highRes; end
    if isfield(fp,'useGrid'); fp.useGrid = p.useGrid; end
    if isfield(fp,'mode'); fp.mode = p.mode; end
    if isfield(fp,'noFlowOutwardOnBorder'); fp.noFlowOutwardOnBorder = p.noFlowOutwardOnBorder; end
    if isfield(fp,'trackSuccessively'); fp.trackSuccessively = p.trackSuccessively; end

    proc2.setPara(fp);

    % --- Process 3: Displacement post-processing (delete if exists) ---
    if ~isempty(tfmPack.processes_{3})
        MD.deleteProcess(tfmPack.processes_{3})
        MD.save
    end

    % --- Process 4: Force reconstruction ---
    if isempty(tfmPack.processes_{4}), tfmPack.createDefaultProcess(4); end
    proc4 = tfmPack.getProcess(4);
    fp = proc4.funParams_;

    if isfield(fp,'YoungModulus'); fp.YoungModulus = p.YoungModulus; end
    if isfield(fp,'regParam');     fp.regParam     = p.regParam; end
    if isfield(fp,'method');       fp.method       = p.method; end
    if isfield(fp,'solMethodBEM'); fp.solMethodBEM = p.solMethodBEM; end
    if isfield(fp,'useLcurve');    fp.useLcurve    = p.useLcurve; end

    if isfield(fp,'basisClassTblPath') && ~isempty(p.basisClassTblPath)
        fp.basisClassTblPath = p.basisClassTblPath;
    end

    proc4.setPara(fp);

    % --- Process 5: Strain energy (if present) ---
    if numel(tfmPack.processes_) >= 5
        if isempty(tfmPack.processes_{5}), tfmPack.createDefaultProcess(5); end
        proc5 = tfmPack.getProcess(5);
        fp = proc5.funParams_;
        proc5.setPara(fp);
    end

    MD.save;
    fprintf('[setupTFMPackage] %d/%d done: %s\n', i, nMovies, MD.outputDirectory_);
end

end

%% ===================== Helper: resolve ref bestZ path =====================
function refImgPath = resolveRefBestZPath(mdOutDir, strict)
% Looks for reference bestZ image under:
%   <mdOutDir>/reference/
% Priority:
%   1) ref_beads_bestZ.tif
%   2) ref_beads_bestZ_series*.tif
%   3) *bestZ*.tif
% strict=true => error on missing or ambiguous matches

if nargin < 2, strict = true; end

refDir = fullfile(mdOutDir, 'reference');
assert(exist(refDir,'dir')==7, 'Reference folder missing: %s', refDir);

p1 = fullfile(refDir, 'ref_beads_bestZ.tif');
if exist(p1,'file')==2
    refImgPath = p1;
    return;
end

cands = [];
d = dir(fullfile(refDir, 'ref_beads_bestZ_series*.tif'));
d = d(~[d.isdir]);
if ~isempty(d), cands = d; end

if isempty(cands)
    d = dir(fullfile(refDir, '*bestZ*.tif'));
    d = d(~[d.isdir]);
    cands = d;
end

if isempty(cands)
    msg = sprintf('Reference image missing in: %s (expected ref_beads_bestZ.tif or ref_beads_bestZ_series*.tif)', refDir);
    if strict, error(msg); else, warning(msg); end
    refImgPath = p1; % placeholder
    return;
end

if numel(cands) > 1
    % If multiple, try to pick the best candidate by name score; still error if strict.
    names = string({cands.name});
    score = zeros(numel(names),1);
    score = score + 5*contains(lower(names),'ref_beads_bestz_series');
    score = score + 3*contains(lower(names),'beads');
    score = score + 3*contains(lower(names),'bestz');

    [~,idx] = max(score);
    picked = cands(idx);

    if strict
        error('Ambiguous reference bestZ files in %s. Found: %s. Please keep only one, or rename to ref_beads_bestZ.tif', ...
            refDir, strjoin(names, ", "));
    else
        warning('Multiple candidates in %s. Using: %s', refDir, picked.name);
        refImgPath = fullfile(picked.folder, picked.name);
        return;
    end
else
    picked = cands(1);
    refImgPath = fullfile(picked.folder, picked.name);
end

end

function p = defaultTFMParams()
% tfmBatch.m ??? ??? + ??

p = struct();
p.beadChan = 2;
p.usfac = 20;

% Displacement params
p.alpha = 0.05;
p.minCorLength = 17;
p.addNonLocMaxBeads = true;
p.maxFlowSpeed = 20;
p.highRes = true;
p.useGrid = true;
p.mode = 'accurate';
p.noFlowOutwardOnBorder = 1;
p.trackSuccessively = false;

% Post-proc
p.outlierThreshold = 2;
p.fillVectors = false;

% Force recon
p.YoungModulus = 2700;   % Pa
p.regParam = 0.005;      % <<< enforce your desired value here
p.method = 'FTTC';
p.solMethodBEM = 'QR';
p.useLcurve = false;
p.basisClassTblPath = '';

% Metadata defaults
p.timeInterval = 120;  % seconds per frame
p.pixelSize = 325;     % nm/pixel
end