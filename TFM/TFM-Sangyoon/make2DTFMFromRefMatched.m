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
    % Step 1: Load movie bead template + identify well-focused tiles
    %
    % When only part of the movie is in-focus (e.g. left side in P07),
    % whole-image NCC is diluted by the blurry region, leading to wrong
    % bestZ. Solution: tile the template, score each tile's focus, and
    % use only the top-K% tiles for NCC against the ref z-stack.
    %
    % tileSize   : tile size in pixels (after downsampling)
    % focusTilePct: keep tiles in the top this% of focus scores (0-100)
    % ==============================================================
    tileSize     = 64;    % pixels per tile side
    focusTilePct = 30;    % use top 30% sharpest tiles

    movBeadPlane = movR.getIndex(0, beadChan-1, tplFrame-1) + 1;
    if verbose
        fprintf('  [Step 1] Movie template: ch%d (0-based idx=%d), frame=%d, plane=%d\n', ...
            beadChan, beadChan-1, tplFrame, movBeadPlane);
    end
    I0_raw = double(bfGetPlane(movR, movBeadPlane));
    if ~isempty(cropROI), I0_raw = imcrop(I0_raw, cropROI); end
    if ds > 1, I0_raw = imresize(I0_raw, 1/ds); end
    I0 = normalizeForScoring(I0_raw);

    % Tile the movie template and score each tile
    [H0, W0] = size(I0);
    [tileCoords, tileFocusScores] = computeTileFocusScores(I0_raw, tileSize);
    nTiles = size(tileCoords, 1);

    % Select top-K% tiles by focus score
    focusThresh   = prctile(tileFocusScores, 100 - focusTilePct);
    selectedTiles = tileFocusScores >= focusThresh;
    nSelected     = sum(selectedTiles);

    if verbose
        fprintf('  [Step 1] Tile analysis: %d tiles total, %d selected (top %.0f%%)\n', ...
            nTiles, nSelected, focusTilePct);
        % Show where well-focused tiles are (column distribution)
        selCoords = tileCoords(selectedTiles, :);
        fprintf('  [Step 1] Well-focused tile column range: %d-%d px (of %d)\n', ...
            min(selCoords(:,2)), max(selCoords(:,2)+tileSize-1), W0);
    end
    % Note: tileCoords and selectedTiles are in movie-template pixel space.
    % They will be recomputed in ref pixel space after size check in Step 2.

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

    % If movie and ref have different sizes (e.g. 881 vs 428),
    % downsample I0 to ref size for shift estimation in Step 2.
    % (NCC and phase correlation require same-size images)
    [Href_native, Wref_native] = size(refImgs{1});
    if size(I0,1) ~= Href_native || size(I0,2) ~= Wref_native
        I0_forShift = imresize(I0, [Href_native Wref_native], 'bicubic');
        if verbose
            fprintf('  [Step 2] Downsampled movie template %dx%d -> %dx%d for shift estimation\\n', ...
                size(I0,1), size(I0,2), Href_native, Wref_native);
        end
    else
        I0_forShift = I0;
    end

    % Recompute tile selection in ref pixel space (I0_forShift)
    [tileCoords, tileFocusScores] = computeTileFocusScores(I0_forShift, tileSize);
    nTiles        = size(tileCoords, 1);
    focusThresh   = prctile(tileFocusScores, 100 - focusTilePct);
    selectedTiles = tileFocusScores >= focusThresh;
    nSelected     = sum(selectedTiles);
    if verbose
        fprintf('  [Step 2] Tile recomputed in ref space: %d tiles, %d selected\n', nTiles, nSelected);
    end

    % --- 2b) Estimate global XY shift using ALL z-layers -> median ---
    %
    % PROBLEM with top-N tenengrad approach:
    %   Near-coverslip z-layers (z=1-5) have bright reflections that give
    %   artificially high Tenengrad scores. These are NOT bead focal planes.
    %   Using them for shift estimation gives wrong/noisy shifts, making
    %   NCC scores uniformly low (~0.035) for all z-layers.
    %
    % FIX: Compute phase correlation for ALL z-layers, take median shift.
    %   - In-focus bead layers give consistent, correct shifts
    %   - Out-of-focus/artifact layers give noisy shifts (outliers)
    %   - Median robustly rejects outliers -> correct global shift
    %
    % Optimization: skip layers with very low intensity (empty/out-of-range)
    %   to avoid spending time on pure-noise phase correlations.
    dyVec = nan(Zref, 1);
    dxVec = nan(Zref, 1);
    meanIntRef = cellfun(@(J) mean(J(:)), refImgs);
    intensityThresh = prctile(meanIntRef, 10);  % skip bottom 10% lowest-intensity layers

    for z = 1:Zref
        if meanIntRef(z) < intensityThresh, continue; end
        [I0c, Jc] = matchSize(I0_forShift, refImgs{z});
        [dyVec(z), dxVec(z)] = phaseCorrelationShift(I0c, Jc);
    end
    dy_global = round(median(dyVec(~isnan(dyVec))));
    dx_global = round(median(dxVec(~isnan(dxVec))));

    validDy  = dyVec(~isnan(dyVec));
    validDx  = dxVec(~isnan(dxVec));
    nUsed    = numel(validDy);
    dyStd    = std(validDy);
    dxStd    = std(validDx);

    if verbose
        fprintf('  Global XY shift (median over %d/%d z-layers): dy=%d dx=%d\n', ...
            nUsed, Zref, dy_global, dx_global);
        fprintf('  Shift consistency (std): dy_std=%.1f dx_std=%.1f (lower=more consistent)\n', ...
            dyStd, dxStd);
    end

    % ---- Alias correction when shift std is large ----
    % Large dx_std or dy_std means phase correlation is giving inconsistent
    % results across z-layers, typically because the true shift exceeds N/2
    % (image half-size) and circular aliasing causes some layers to report
    % shift and others to report shift +/- N.
    %
    % Strategy: test the original median shift AND its aliased variants
    % (shift +/- image_width or +/- image_height).
    % For each candidate, compute mean NCC over a sample of z-layers.
    % Pick the candidate that maximises mean NCC = most physically correct.
    [H_img, W_img] = size(I0);
    shiftStdThresh = 30;  % px: above this, suspect aliasing

    if dyStd > shiftStdThresh || dxStd > shiftStdThresh
        if verbose
            fprintf('  [AliasCheck] High shift std - testing aliased variants...\n');
        end

        % Candidate shifts: original + aliases in each axis
        dy_cands = unique([dy_global, dy_global + H_img, dy_global - H_img]);
        dx_cands = unique([dx_global, dx_global + W_img, dx_global - W_img]);

        % Sample z-layers evenly for fast NCC testing
        nTest    = min(15, Zref);
        testZs   = round(linspace(1, Zref, nTest));

        best_dy  = dy_global;
        best_dx  = dx_global;
        best_ncc = -Inf;

        for cdy = dy_cands
            for cdx = dx_cands
                nccs = nan(nTest, 1);
                for ki = 1:nTest
                    z_t = testZs(ki);
                    Jref_t = refImgs{z_t};
                    tileNccs_a = nan(nTiles,1);
                    for ti2 = 1:nTiles
                        if ~selectedTiles(ti2), continue; end
                        r0=tileCoords(ti2,1); c0=tileCoords(ti2,2);
                        r1=min(r0+tileSize-1,Href_native); c1=min(c0+tileSize-1,Wref_native);
                        rr0=r0-cdy; rr1=r1-cdy; cc0=c0-cdx; cc1=c1-cdx;
                        if rr0<1||cc0<1||rr1>size(Jref_t,1)||cc1>size(Jref_t,2), continue; end
                        Ip=I0_forShift(r0:r1,c0:c1); Jp=Jref_t(rr0:rr1,cc0:cc1);
                        if numel(Ip)<9||std(Ip(:))<eps||std(Jp(:))<eps, continue; end
                        tileNccs_a(ti2) = corr2(Ip,Jp);
                    end
                    vn = tileNccs_a(~isnan(tileNccs_a));
                    if numel(vn)>=3, nccs(ki) = mean(vn); end
                end
                meanNcc = nanmean(nccs);
                if meanNcc > best_ncc
                    best_ncc = meanNcc;
                    best_dy  = cdy;
                    best_dx  = cdx;
                end
            end
        end

        if best_dy ~= dy_global || best_dx ~= dx_global
            if verbose
                fprintf('  [AliasCheck] Corrected: dy %d->%d  dx %d->%d  (mean_ncc=%.3f)\n', ...
                    dy_global, best_dy, dx_global, best_dx, best_ncc);
            end
            dy_global = best_dy;
            dx_global = best_dx;
        else
            if verbose
                fprintf('  [AliasCheck] Original shift confirmed (mean_ncc=%.3f)\n', best_ncc);
            end
        end
    end

    % --- 2c-d) Score all z-layers using tile-based partial NCC ---
    %
    % For each z-layer, compute NCC only in the well-focused tiles of I0
    % (identified in Step 1). This ignores blurry regions and focuses on
    % the area where beads are sharp (e.g. left side in P07).
    %
    % For each selected tile at (r0, c0) in the movie:
    %   - Find corresponding region in ref after applying global shift
    %   - Compute corr2 for that patch pair
    % Final score for z = mean corr2 across all selected tiles.
    scores = nan(Zref, 1);
    for z = 1:Zref
        Jref = refImgs{z};

        tileNccs = nan(nTiles, 1);
        for ti = 1:nTiles
            if ~selectedTiles(ti), continue; end

            r0 = tileCoords(ti, 1);
            c0 = tileCoords(ti, 2);
            r1 = min(r0 + tileSize - 1, Href_native);
            c1 = min(c0 + tileSize - 1, Wref_native);

            % Corresponding ref tile with global shift applied
            rr0 = r0 - dy_global;  rr1 = r1 - dy_global;
            cc0 = c0 - dx_global;  cc1 = c1 - dx_global;

            % Check bounds in ref
            if rr0 < 1 || cc0 < 1 || rr1 > size(Jref,1) || cc1 > size(Jref,2)
                continue;
            end

            Ipatch = I0_forShift(r0:r1, c0:c1);
            Jpatch = Jref(rr0:rr1, cc0:cc1);

            if numel(Ipatch) < 9, continue; end
            if std(Ipatch(:)) < eps || std(Jpatch(:)) < eps, continue; end

            tileNccs(ti) = corr2(Ipatch, Jpatch);
        end

        validNccs = tileNccs(~isnan(tileNccs));
        if numel(validNccs) >= 3
            scores(z) = mean(validNccs);
        else
            scores(z) = -Inf;
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
    % Step 3: Save ref image + QC score profile
    %
    % Size mismatch handling (e.g. ref 428x428, movie 881x881):
    %   If the movie template is larger than the ref slice, apply
    %   weighted focal stack fusion using bestZ ± fusionHalfWidth slices,
    %   weighted by per-tile hybrid focus scores, then bicubic upsample
    %   to match the movie template size.
    %
    %   Fusion is skipped when sizes match (no mismatch), saving
    %   the single bestZ slice directly.
    %
    % fusionHalfWidth: number of z-slices on each side of bestZ to include
    % ==============================================================
    if writeRef
        refDir = fullfile(outRoot, 'reference');
        if ~exist(refDir,'dir'); mkdir(refDir); end

        % Load single bestZ ref slice (native resolution)
        refBestPlane = refR.getIndex(bestZref-1, refBeadChan-1, rFrame-1) + 1;
        Ibest        = double(bfGetPlane(refR, refBestPlane));
        if ~isempty(cropROI), Ibest = imcrop(Ibest, cropROI); end

        % Check size mismatch with movie template
        [Href, Wref] = size(Ibest);
        [Hmov, Wmov] = size(I0_raw);   % original (non-normalized) movie template
        sizeMismatch = (Href ~= Hmov) || (Wref ~= Wmov);

        if sizeMismatch
            if verbose
                fprintf('  [RefUpscale] Size mismatch: ref=%dx%d, movie=%dx%d\n', ...
                    Href, Wref, Hmov, Wmov);
            end

            % -------------------------------------------------------
            % Exact 2x Fourier-domain upsampling + RL deconvolution
            %
            % RATIONALE:
            %   Bicubic upscale 428->881 changes the FOV (beads appear
            %   larger). The correct approach:
            %     1) Fourier-domain 2x upsample: 428->856 (exact sinc
            %        interpolation, preserves all spatial frequencies,
            %        same FOV at half pixel size)
            %     2) Mild Richardson-Lucy deconvolution to sharpen beads
            %        (Airyscan PSF ~0.5px at 85nm/px, so gentle)
            %     3) Zero-pad symmetrically 856->Hmov (e.g. 881):
            %        adds (Hmov-856)/2 px border of zeros
            %        ? FOV matches movie exactly
            %
            % Fourier upsample is lossless (no blurring, no stretching).
            % -------------------------------------------------------

            % Step A: Fourier-domain exact 2x upsampling (428 -> 856)
            H2x = Href * 2;
            W2x = Wref * 2;
            IrefFourier = fourierUpsample2x(Ibest, H2x, W2x);

            % Step B: Richardson-Lucy deconvolution (sharpen beads)
            % PSF for Airyscan-processed image at 85nm/px:
            %   effective FWHM ~ 140nm -> sigma ~ 0.83px at 85nm/px
            %   After 2x upsample, pixel is 42.5nm -> sigma ~ 1.65px
            rlIterations = 8;
            psfSigma     = 1.65;   % px in 2x-upsampled space
            IrefDeconv   = richardsonLucy(IrefFourier, psfSigma, rlIterations);

            % Step C: Symmetric zero-pad to match movie size (856 -> Hmov)
            padH = Hmov - H2x;   % total padding rows (can be 0 or positive)
            padW = Wmov - W2x;
            if padH >= 0 && padW >= 0
                padT = floor(padH/2);  padB = ceil(padH/2);
                padL = floor(padW/2);  padR = ceil(padW/2);
                IrefOut = padarray(IrefDeconv, [padT padL], 0, 'pre');
                IrefOut = padarray(IrefOut,    [padB padR], 0, 'post');
            else
                % Movie smaller than 2x ref (unusual): center-crop
                r0 = floor(-padH/2)+1; r1 = r0 + Hmov - 1;
                c0 = floor(-padW/2)+1; c1 = c0 + Wmov - 1;
                IrefOut = IrefDeconv(r0:r1, c0:c1);
            end

            if verbose
                fprintf('  [RefUpscale] Fourier 2x upsample: %dx%d -> %dx%d\n', ...
                    Href, Wref, H2x, W2x);
                fprintf('  [RefUpscale] RL deconvolution: %d iter, PSF sigma=%.2f px\n', ...
                    rlIterations, psfSigma);
                fprintf('  [RefUpscale] Zero-pad: %dx%d -> %dx%d (pad T/B=%d/%d L/R=%d/%d)\n', ...
                    H2x, W2x, Hmov, Wmov, padT, padB, padL, padR);
            end
        else
            % No size mismatch: save bestZ slice directly
            IrefOut = Ibest;
            if verbose
                fprintf('  [RefFusion] No size mismatch, saving bestZ slice directly.\n');
            end
        end

        imwrite(castToUInt16(IrefOut), fullfile(refDir, 'ref_beads_matched.tif'));
        save(fullfile(refDir, 'z_score_profile.mat'), ...
            'scores', 'focusScores', 'dy_global', 'dx_global', 'bestZref', ...
            'sizeMismatch');

        if verbose
            fprintf('  ref image saved (bestZref=%d): %s\n', bestZref, refDir);
        end
    end

    % ==============================================================
    % Step 4: Create MovieData via bfImport, then enrich channel
    %         metadata from BioFormats store.
    %
    % bfImport creates channels with BioFormatsReader but often misses
    % emission/excitation wavelengths (the "Fluorophore not supported"
    % warning). We re-read the store and set them explicitly.
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

    % Enrich channel metadata from BioFormats store
    try
        md = enrichChannelMetadata(md, movieFile, s, beadChan, verbose);
    catch ME
        warning(ME.identifier, 'enrichChannelMetadata failed: %s', ME.message);
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

function [tileCoords, scores] = computeTileFocusScores(I, tileSize)
% Divide image I into non-overlapping tiles and compute a hybrid focus
% score for each tile using geometric mean of Tenengrad and FFT high-freq
% energy ratio.
%
% For ~1px beads (40nm bead at 85nm/px):
%   Tenengrad : detects sharp edges -> high when bead is in-focus
%   FFT score : high-freq energy ratio -> high when PSF is narrow (in-focus)
%               out-of-focus PSF spreads energy to low-freq -> score drops
%   Hybrid    : geometric mean -> BOTH conditions must be met
%               (avoids high score from background texture alone)
[H, W] = size(I);
rStarts = 1:tileSize:H;
cStarts = 1:tileSize:W;
nTiles  = numel(rStarts) * numel(cStarts);
tileCoords = zeros(nTiles, 2);
scores     = zeros(nTiles, 1);
ti = 0;
for ri = 1:numel(rStarts)
    for ci = 1:numel(cStarts)
        ti = ti + 1;
        r0 = rStarts(ri); r1 = min(r0+tileSize-1, H);
        c0 = cStarts(ci); c1 = min(c0+tileSize-1, W);
        tileCoords(ti,:) = [r0, c0];
        patch = double(I(r0:r1, c0:c1));

        % --- Tenengrad ---
        [Gx, Gy] = imgradientxy(patch);
        tScore = mean(Gx(:).^2 + Gy(:).^2);

        % --- FFT high-frequency energy ratio ---
        % For ~1px point-source beads, in-focus PSF is narrow -> energy
        % spread across all frequencies including high-freq.
        % Out-of-focus PSF expands to 3-5px Gaussian -> low-pass -> high-
        % freq energy drops. Cutoff at Nyquist/4 (tileSize/8).
        fScore = fftHighFreqScore(patch, tileSize);

        % --- Hybrid: geometric mean ---
        % Requires BOTH sharp edges AND high-freq energy to score high.
        % Prevents background texture (high Tenengrad but low FFT) or
        % uniform bright regions (high FFT but low Tenengrad) from
        % being selected as well-focused bead tiles.
        scores(ti) = sqrt(tScore * fScore);
    end
end
end


function score = fftHighFreqScore(patch, tileSize)
% High-frequency energy ratio in the FFT power spectrum.
% Cutoff: r > tileSize/8 (Nyquist/4).
F  = fft2(patch - mean(patch(:)));
P  = abs(fftshift(F)).^2;
[H, W] = size(P);
[cols, rows] = meshgrid(1:W, 1:H);
cx = W/2 + 1;  cy = H/2 + 1;
r_px     = sqrt((rows-cy).^2 + (cols-cx).^2);
r_cutoff = tileSize / 8;   % ~8px cutoff for 64px tile
highMask = r_px > r_cutoff;
score    = sum(P(highMask)) / (sum(P(:)) + eps);
end

function emNm = guessEmissionFromFilename(fname, chanIdx, beadChan)
% Guess emission wavelength from movie filename and channel index.
%
% Strategy:
%   1) Look for explicit wavelength numbers: "647", "561", "488" etc.
%      near bead-related keywords
%   2) Look for color keywords: R/Red -> 625nm, G/Green -> 515nm,
%      B/Blue -> 460nm, Far-red/Cy5/647 -> 647nm
%   3) Only apply color guess to the bead channel; leave other channels
%      as NaN (they likely have their own metadata)

emNm = NaN;
isBead = (chanIdx == beadChan);

% Step 1: explicit wavelength numbers in filename (e.g. "647nm", "561nm")
% Search near bead-related context
toks = regexp(fname, '(\d{3})(?:\s*nm)?', 'tokens');
for k = 1:numel(toks)
    v = str2double(toks{k}{1});
    % plausible emission range: 400-800nm
    if v >= 400 && v <= 800
        emNm = v;
        return;  % use first plausible wavelength found
    end
end

% Step 2: color keywords (apply only to bead channel)
if ~isBead, return; end

fname_lower = lower(fname);

% Far-red / Cy5 / Alexa647
if contains(fname_lower, {'cy5','alexa647','alexa 647','a647','far.red','farred'})
    emNm = 668; return;
end
% Red / Cy3 / Alexa555/568/594
if contains(fname_lower, {'bead-r','bead_r','_r-','_r_','-r-', ...
        'cy3','alexa555','alexa 555','alexa568','alexa594','red','rhodamine'})
    emNm = 625; return;
end
% Green
if contains(fname_lower, {'bead-g','bead_g','_g-','_g_','-g-', ...
        'gfp','fitc','alexa488','alexa 488','green'})
    emNm = 515; return;
end
% Blue
if contains(fname_lower, {'bead-b','bead_b','_b-','_b_','-b-', ...
        'dapi','hoechst','blue','cf405','alexa405'})
    emNm = 460; return;
end
end


function md = enrichChannelMetadata(md, movieFile, seriesIdx, beadChan, verbose)
% enrichChannelMetadata
% Reads channel properties from BioFormats metadata store and sets them
% in the MD Channel objects. Fills gaps using channel name heuristics.
%
% Priority for emission wavelength:
%   1) OME store: getChannelEmissionWavelength()
%   2) Channel name -> name2wavelength() heuristic
%   3) Leave as-is (do not overwrite if already set)
%
% Sets per-channel: emissionWavelength_, excitationWavelength_, name_,
%                   camBitdepth_, numericalAperture_

if nargin < 5, verbose = false; end

try
    reader = bfGetReader(char(movieFile));
    reader.setSeries(seriesIdx);
    store  = reader.getMetadataStore();
    nChan  = reader.getSizeC();

    for c0 = 0:nChan-1
        ci = c0 + 1;
        if ci > numel(md.channels_), break; end
        ch = md.channels_(ci);

        % ---- Channel name ----
        try
            cname = char(store.getChannelName(seriesIdx, c0));
            if ~isempty(cname) && ~isempty(strtrim(cname))
                if isprop(ch,'name_') || isfield(ch,'name_')
                    ch.name_ = cname;
                end
            end
        catch, end

        % ---- Emission wavelength (nm) ----
        emSet = false;
        try
            w = store.getChannelEmissionWavelength(seriesIdx, c0);
            if ~isempty(w)
                emNm = w.value().doubleValue();
                if emNm > 0 && isprop(ch,'emissionWavelength_')
                    ch.emissionWavelength_ = emNm;
                    emSet = true;
                end
            end
        catch, end

        % Fallback: name2wavelength from channel name
        if ~emSet
            try
                cname = '';
                if isprop(ch,'name_'), cname = ch.name_; end
                if isempty(cname)
                    cname = char(store.getChannelName(seriesIdx, c0));
                end
                if ~isempty(cname)
                    wNm = name2wavelength(cname);
                    if ~isnan(wNm) && wNm > 0 && isprop(ch,'emissionWavelength_')
                        ch.emissionWavelength_ = wNm * 1e9;  % name2wavelength returns meters
                        emSet = true;
                    end
                end
            catch, end
        end

        % Fallback 2: filename-based color heuristic
        % Looks for color indicators in the movie filename:
        %   -R / BEAD-R / Red  -> ~625nm (generic red bead)
        %   -G / BEAD-G / Green -> ~515nm
        %   -B / BEAD-B / Blue  -> ~460nm
        %   explicit nm value (e.g. "647nm", "561") -> use directly
        if ~emSet
            try
                [~, fname] = fileparts(movieFile);
                emGuess = guessEmissionFromFilename(fname, ci, beadChan);
                if ~isnan(emGuess) && isprop(ch,'emissionWavelength_')
                    ch.emissionWavelength_ = emGuess;
                    emSet = true;
                    if verbose
                        fprintf('  Ch%d: emission guessed from filename: %.0fnm\\n', ci, emGuess);
                    end
                end
            catch, end
        end
        try
            wx = store.getChannelExcitationWavelength(seriesIdx, c0);
            if ~isempty(wx)
                exNm = wx.value().doubleValue();
                if exNm > 0 && isprop(ch,'excitationWavelength_')
                    ch.excitationWavelength_ = exNm;
                end
            end
        catch, end

        % ---- Camera bit depth ----
        try
            bd = store.getPixelsSignificantBits(seriesIdx);
            if ~isempty(bd) && isprop(ch,'camBitdepth_')
                ch.camBitdepth_ = double(bd.getValue());
            end
        catch, end

        % ---- Numerical aperture ----
        try
            objSettingsID = store.getObjectiveSettingsID(seriesIdx);
            if ~isempty(objSettingsID)
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
        catch, end

        md.channels_(ci) = ch;

        if verbose
            emVal = nan; exVal = nan; bdVal = nan; naVal = nan;
            try, emVal = ch.emissionWavelength_;  catch, end
            try, exVal = ch.excitationWavelength_; catch, end
            try, bdVal = ch.camBitdepth_;          catch, end
            try, naVal = ch.numericalAperture_;    catch, end
            isBead = (ci == beadChan);
            fprintf('  Ch%d%s: em=%.0fnm ex=%.0fnm bits=%g NA=%.2f\n', ...
                ci, iif(isBead,' [BEAD]',''), emVal, exVal, bdVal, naVal);
        end
    end

    try; reader.close(); catch, end
catch ME
    warning(ME.identifier, 'enrichChannelMetadata: %s', ME.message);
end
end


function v = iif(cond, a, b)
if cond, v = a; else, v = b; end
end


function Iout = fourierUpsample2x(I, H2x, W2x)
% Exact 2x upsampling via Fourier zero-padding (sinc interpolation).
% Places the original spectrum in the center of a larger array,
% zero-pads the high-frequency bands, then inverse FFT.
% This is equivalent to ideal sinc interpolation: lossless, no blurring.
[H, W] = size(I);
F  = fftshift(fft2(double(I)));

% Place original spectrum in center of zero-padded array
Fpad = zeros(H2x, W2x);
rOff = floor((H2x - H) / 2);
cOff = floor((W2x - W) / 2);
Fpad(rOff+1:rOff+H, cOff+1:cOff+W) = F;

% Scale to preserve mean intensity (energy conservation for 2x upsample)
Iout = real(ifft2(ifftshift(Fpad))) * (H2x * W2x) / (H * W);
Iout = max(Iout, 0);  % clip negative values (ringing artefact suppression)
end


function Iout = richardsonLucy(I, psfSigma, nIter)
% Richardson-Lucy deconvolution with Gaussian PSF.
% psfSigma: PSF standard deviation in pixels
% nIter   : number of RL iterations (5-10 typically sufficient)
%
% RL update: I_{k+1} = I_k * (conv((Y / conv(I_k, PSF)), PSF_flip))
% For symmetric PSF, PSF_flip = PSF.
I    = double(I);
I    = max(I, eps);   % avoid division by zero

% Build Gaussian PSF kernel (7-sigma wide to avoid edge artefacts)
kRadius = ceil(3.5 * psfSigma);
[kx, ky] = meshgrid(-kRadius:kRadius, -kRadius:kRadius);
psf = exp(-(kx.^2 + ky.^2) / (2 * psfSigma^2));
psf = psf / sum(psf(:));

Iest = I;  % initial estimate = observed image
for k = 1:nIter
    Iblur  = imfilter(Iest, psf, 'symmetric', 'same');
    Iblur  = max(Iblur, eps);
    ratio  = I ./ Iblur;
    Iest   = Iest .* imfilter(ratio, psf, 'symmetric', 'same');
    Iest   = max(Iest, 0);
end
Iout = Iest;
end


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