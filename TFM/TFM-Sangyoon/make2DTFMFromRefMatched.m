function mdList = make2DTFMFromRefMatched(movieFile, refFile, outRoot, varargin)
% make2DTFMFromRefMatched
%
% Finds the ref z-stack layer most similar to the movie bead channel image
% (drift-aware: phase correlation + overlap NCC), saves the matched ref image,
% then creates MovieData directly from the CZI via bfImport (no TIFF extraction).
%
% bfImport handles:
%   - Channel creation from CZI path directly (BioFormatsReader)
%   - All metadata: pixelSize, timeInterval, camBitdepth, NA, emission wavelength
%
% Steps per series:
%   1) Load movie bead channel template (templateFrame, z=1)
%   2) Phase correlation + overlap NCC across all ref z-layers -> bestZref
%   3) Save ref image at bestZref + z_score_profile.mat for QC
%   4) bfImport(movieFile, outputDirectory) -> MovieData with full metadata
%
% Inputs:
%   movieFile : cell-present movie CZI
%   refFile   : reference bead z-stack CZI
%   outRoot   : root output folder (each series gets a subfolder)
%
% Name-Value options:
%   'beadChan'      : bead channel index, 1-based (default: 2)
%   'templateFrame' : frame in movie used as matching template (default: 1)
%   'refFrame'      : time frame in ref z-stack to scan (default: 1)
%   'cropROI'       : [x y w h] crop before scoring, [] = full image (default: [])
%   'downsample'    : downsampling factor for scoring speed, >=1 (default: 1)
%   'maxShiftFrac'  : max allowed XY shift as fraction of image size (default: 0.9)
%   'verbose'       : true/false (default: true)
%   'writeRefSlice' : save bestZref image to reference/ subfolder (default: true)
%   'overwrite'     : reprocess even if movieData.mat already exists (default: false)
%
% Output:
%   mdList : MovieData array (one per series)

ip = inputParser;
ip.addParameter('beadChan',      2,     @(x) isnumeric(x) && isscalar(x) && x>=1);
ip.addParameter('refBeadChan',   [],    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>=1));
% refBeadChan: bead channel index in the ref z-stack (1-based).
%   Default [] -> uses the same value as beadChan.
%   Set explicitly when the ref z-stack has a different channel order than the movie.
%   e.g. movie: ch1=green(talin), ch2=red(bead) -> beadChan=2
%        ref  : ch1=red(bead) only              -> refBeadChan=1
ip.addParameter('templateFrame', 1,     @(x) isnumeric(x) && isscalar(x) && x>=1);
ip.addParameter('refFrame',      1,     @(x) isnumeric(x) && isscalar(x) && x>=1);
ip.addParameter('cropROI',       [],    @(x) isempty(x) || (isnumeric(x) && numel(x)==4));
ip.addParameter('downsample',    1,     @(x) isnumeric(x) && isscalar(x) && x>=1);
ip.addParameter('maxShiftFrac',  0.9,   @(x) isnumeric(x) && isscalar(x) && x>0 && x<=1);
ip.addParameter('scoreThreshold', 0.3,  @(x) isnumeric(x) && isscalar(x));
% scoreThreshold: minimum corr2 score to accept bestZref automatically.
%   If best score < threshold, the function shows the top candidates and
%   prompts the user to enter the correct z manually.
%   Set to -Inf to always accept without prompting.
ip.addParameter('verbose',       true,  @(x) islogical(x) && isscalar(x));
ip.addParameter('writeRefSlice', true,  @(x) islogical(x) && isscalar(x));
ip.addParameter('overwrite',     false, @(x) islogical(x) && isscalar(x));
ip.parse(varargin{:});

beadChan      = ip.Results.beadChan;
refBeadChan   = ip.Results.refBeadChan;
if isempty(refBeadChan), refBeadChan = beadChan; end  % default: same as movie
templateFrame = ip.Results.templateFrame;
refFrame      = ip.Results.refFrame;
cropROI       = ip.Results.cropROI;
ds            = ip.Results.downsample;
maxShiftFrac   = ip.Results.maxShiftFrac;
scoreThreshold = ip.Results.scoreThreshold;
verbose       = ip.Results.verbose;
writeRef      = ip.Results.writeRefSlice;
overwrite     = ip.Results.overwrite;

movieFile = char(movieFile);
refFile   = char(refFile);
outRoot   = char(outRoot);
if ~exist(outRoot,'dir'); mkdir(outRoot); end

if verbose
    fprintf('[TFMprep] movieFile   : %s\n', movieFile);
    fprintf('[TFMprep] refFile     : %s\n', refFile);
    fprintf('[TFMprep] outRoot     : %s\n', outRoot);
    fprintf('[TFMprep] beadChan=%d (movie) | refBeadChan=%d (ref z-stack) | templateFrame=%d | refFrame=%d\n', ...
        beadChan, refBeadChan, templateFrame, refFrame);
    fprintf('[TFMprep] maxShiftFrac=%.2f | downsample=%d | overwrite=%d\n', ...
        maxShiftFrac, ds, overwrite);
end

% Open readers (ref for scoring; movie only for loading the template image)
refR = bfGetReader(refFile);
movR = bfGetReader(movieFile);

nSref = refR.getSeriesCount();
nSmov = movR.getSeriesCount();
nS    = min(nSmov, nSref);

if verbose
    fprintf('[TFMprep] movie series=%d | ref series=%d -> processing %d\n', nSmov, nSref, nS);
end

mdCell = cell(nS, 1);

for s = 0:nS-1
    refR.setSeries(s);
    movR.setSeries(s);

    % Series label for output subfolder
    posName = sprintf('series_%03d', s);
    try
        nm = movR.getMetadataStore().getImageName(s);
        if ~isempty(nm), posName = char(nm); end
    catch
    end
    posNameClean = regexprep(posName, '[^\w\-\.]+', '_');
    outDir = fullfile(outRoot, posNameClean);

    % Overwrite check
    [~, movName, movExt] = fileparts(movieFile);
    mdFile = fullfile(outDir, [movName movExt '.mat']);
    if nS > 1
        sStr   = num2str(s+1, ['_s%0' num2str(floor(log10(nS))+1) '.f']);
        mdFile = fullfile(outDir, [movName sStr '.mat']);
    end

    if ~overwrite && exist(mdFile, 'file')
        if verbose
            fprintf('[TFMprep] series=%d: movieData.mat exists, loading (overwrite=false).\n', s);
        end
        tmp = load(mdFile);
        fn = fieldnames(tmp);
        mdCell{s+1} = tmp.(fn{1});
        continue;
    end

    Zref = refR.getSizeZ();  Cref = refR.getSizeC();
    Zmov = movR.getSizeZ();  C    = movR.getSizeC();  T = movR.getSizeT();

    if verbose
        fprintf('\n[TFMprep] series=%d (%s)\n', s, posName);
        fprintf('  movie: Z=%d C=%d T=%d | ref: Z=%d C=%d\n', Zmov, C, T, Zref, Cref);
    end

    if Cref < refBeadChan || C < beadChan || Zmov<1 || T<1 || Zref<1
        warning('series %d: invalid dimensions or beadChan out of range. Skipping.', s);
        continue;
    end

    tplFrame = min(templateFrame, T);
    rFrame   = min(refFrame, refR.getSizeT());

    % ==============================================================
    % Step 1: Load movie bead template (z=1 for single-plane movie)
    % ==============================================================
    movBeadPlane = movR.getIndex(0, beadChan-1, tplFrame-1) + 1;
    if verbose
        fprintf('  [Step 1] Movie template: ch%d (0-based idx=%d), frame=%d, plane=%d\n', ...
            beadChan, beadChan-1, tplFrame, movBeadPlane);
    end
    I0 = double(bfGetPlane(movR, movBeadPlane));
    if ~isempty(cropROI), I0 = imcrop(I0, cropROI); end
    if ds > 1, I0 = imresize(I0, 1/ds); end
    I0 = normalizeForScoring(I0);

    % ==============================================================
    % Step 2: Drift-aware scoring across all ref z-layers
    %
    % KEY INSIGHT: XY stage drift is fixed for all z-layers of the same
    % position. Estimating shift per z-layer is unreliable because
    % out-of-focus bead images produce noisy phase correlation results.
    %
    % Strategy:
    %   2a) Compute tenengrad focus score for all ref z-layers
    %   2b) Estimate XY shift from top-N sharpest z-layers -> median shift
    %   2c) Apply ONE fixed shift to all z-layers
    %   2d) Score all z-layers with corr2 over the fixed overlap region
    % ==============================================================
    if verbose
        fprintf('  [Step 2] Ref z-stack scan: ch%d (0-based idx=%d), refFrame=%d, Zref=%d\n', ...
            refBeadChan, refBeadChan-1, rFrame, Zref);
    end

    % --- 2a) Load all ref z-layers and compute focus score ---
    refImgs     = cell(Zref, 1);
    focusScores = zeros(Zref, 1);
    for z = 1:Zref
        refPlane = refR.getIndex(z-1, refBeadChan-1, rFrame-1) + 1;
        J = double(bfGetPlane(refR, refPlane));
        if ~isempty(cropROI), J = imcrop(J, cropROI); end
        if ds > 1, J = imresize(J, 1/ds); end
        refImgs{z}     = normalizeForScoring(J);
        focusScores(z) = tenengradScore(J);
    end

    % --- 2b) Estimate global XY shift from top-N sharpest z-layers ---
    nShift = min(5, Zref);
    [~, sharpIdx] = sort(focusScores, 'descend');
    topIdx = sharpIdx(1:nShift);

    dyVec = zeros(nShift, 1);
    dxVec = zeros(nShift, 1);
    for k = 1:nShift
        [I0c, Jc] = matchSize(I0, refImgs{topIdx(k)});
        [dyVec(k), dxVec(k)] = phaseCorrelationShift(I0c, Jc);
    end
    dy_global = round(median(dyVec));
    dx_global = round(median(dxVec));

    if verbose
        fprintf('  Global XY shift (from top-%d sharpest z): dy=%d dx=%d\n', ...
            nShift, dy_global, dx_global);
        fprintf('  Sharpest z-layers used for shift: ');
        fprintf('z=%d ', topIdx); fprintf('\n');
    end

    % --- 2c-d) Score all z-layers using fixed shift ---
    scores = nan(Zref, 1);
    for z = 1:Zref
        [I0c, Jc] = matchSize(I0, refImgs{z});
        [Iov, Jov] = overlapRegion(I0c, Jc, dy_global, dx_global);

        if isempty(Iov) || numel(Iov) < 100
            scores(z) = -Inf;
            continue;
        end
        if std(Iov(:)) < eps || std(Jov(:)) < eps
            scores(z) = -Inf;
        else
            scores(z) = corr2(Iov, Jov);
        end
    end

    [bestScore, bestZref] = max(scores);

    if verbose
        fprintf('  bestZref=%d  score=%.4f\n', bestZref, bestScore);
        validMask = isfinite(scores);
        if any(validMask)
            tmpScores = scores;
            tmpScores(~validMask) = NaN;
            [sS, sI] = sort(tmpScores, 'descend', 'MissingPlacement', 'last');
            fprintf('  Top-5: ');
            for k = 1:min(5, sum(validMask))
                fprintf('z=%d(%.4f) ', sI(k), sS(k));
            end
            fprintf('\n');
        else
            fprintf('  WARNING: overlap region empty. Check maxShiftFrac (now=%.2f) or cropROI.\n', maxShiftFrac);
        end
    end

    % ==============================================================
    % Confidence check: if best score is below threshold, prompt user
    % ==============================================================
    if bestScore < scoreThreshold
        fprintf('\n  *** LOW CONFIDENCE WARNING ***\n');
        fprintf('  series=%d (%s)\n', s, posName);
        fprintf('  Best score = %.4f (threshold = %.2f)\n', bestScore, scoreThreshold);
        fprintf('  Top-5 candidates:\n');
        validMask = isfinite(scores);
        tmpScores = scores; tmpScores(~validMask) = NaN;
        [sS, sI] = sort(tmpScores, 'descend', 'MissingPlacement', 'last');
        for k = 1:min(5, sum(validMask))
            fprintf('    z=%d  score=%.4f\n', sI(k), sS(k));
        end
        fprintf('  Ref z-stack has %d layers total.\n', Zref);
        userZ = input(sprintf('  Enter correct z manually (1-%d), or press Enter to use z=%d: ', Zref, bestZref));
        if ~isempty(userZ) && isnumeric(userZ) && userZ >= 1 && userZ <= Zref
            bestZref = round(userZ);
            fprintf('  Using user-specified bestZref = %d\n', bestZref);
        else
            fprintf('  Keeping auto-selected bestZref = %d\n', bestZref);
        end
    end

    % ==============================================================
    % Step 3: Save bestZref ref image + QC score profile
    %         Saved directly under outRoot (the label folder, e.g. WT_P04/)
    %         so it is easy to find alongside the MovieData folder.
    % ==============================================================
    if writeRef
        refDir = fullfile(outRoot, 'reference');
        if ~exist(refDir,'dir'); mkdir(refDir); end

        refBestPlane = refR.getIndex(bestZref-1, refBeadChan-1, rFrame-1) + 1;
        Ibest = bfGetPlane(refR, refBestPlane);
        if ~isempty(cropROI), Ibest = imcrop(Ibest, cropROI); end
        imwrite(castToUInt16(Ibest), fullfile(refDir, 'ref_beads_matched.tif'));
        save(fullfile(refDir, 'z_score_profile.mat'), 'scores', 'focusScores', 'dy_global', 'dx_global', 'bestZref');

        if verbose
            fprintf('  ref image saved (bestZref=%d): %s\n', bestZref, refDir);
        end
    end

    % ==============================================================
    % Step 4: Create MovieData via bfImport (CZI read directly)
    %         bfImport sets BioFormatsReader and all metadata automatically.
    % ==============================================================
    if verbose
        fprintf('  Running bfImport for MovieData creation...\n');
    end

    mdArr = bfImport(movieFile, 'outputDirectory', outRoot);

    % bfImport returns one MD per series; pick the right one
    if numel(mdArr) >= s+1
        md = mdArr(s+1);
    else
        md = mdArr(1);
    end

    md.save();

    if verbose
        fprintf('  MovieData saved: %s\n', md.getFullPath());
    end

    mdCell{s+1} = md;
end

% Close readers
try refR.close(); catch, end
try movR.close(); catch, end

mdList = [mdCell{~cellfun(@isempty, mdCell)}];

end

%% ======== Local helper functions ========

function s = tenengradScore(I)
% Tenengrad focus score: mean squared gradient magnitude
[Gx, Gy] = imgradientxy(double(I));
s = mean(Gx(:).^2 + Gy(:).^2);
end

function A = normalizeForScoring(A)
A = double(A);
A = A - mean(A(:));
sd = std(A(:));
if sd > 0, A = A / sd; end
end

function [A2, B2] = matchSize(A, B)
h = min(size(A,1), size(B,1));
w = min(size(A,2), size(B,2));
A2 = A(1:h, 1:w);
B2 = B(1:h, 1:w);
end

function [dy, dx] = phaseCorrelationShift(A, B)
FA = fft2(A);
FB = fft2(B);
R  = FA .* conj(FB);
denom = abs(R);
denom(denom < eps) = eps;
R  = R ./ denom;
pc = real(ifft2(R));

[~, idx] = max(pc(:));
[r, c]   = ind2sub(size(pc), idx);
[H, W]   = size(pc);
if r > H/2, r = r - H; end
if c > W/2, c = c - W; end
dy = r - 1;
dx = c - 1;
end

function [Iov, Jov] = overlapRegion(I, J, dy, dx)
[H, W] = size(I);
if dy >= 0
    Ir = (dy+1):H;    Jr = 1:(H-dy);
else
    Ir = 1:(H+dy);    Jr = (-dy+1):H;
end
if dx >= 0
    Ic = (dx+1):W;    Jc = 1:(W-dx);
else
    Ic = 1:(W+dx);    Jc = (-dx+1):W;
end
if isempty(Ir) || isempty(Ic) || isempty(Jr) || isempty(Jc)
    Iov = []; Jov = []; return;
end
Iov = I(Ir, Ic);
Jov = J(Jr, Jc);
end

function U = castToUInt16(I)
if isinteger(I)
    U = I;
else
    I = double(I);
    I = I - min(I(:));
    mx = max(I(:));
    if mx > 0, I = I / mx; end
    U = uint16(round(I * 65535));
end
end