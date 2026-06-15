function segmentMovieFilopodia(movieData)
%SEGMENTMOVIEFILOPODIA  Process 1 wrapper. Cell body mask + steerable maps.
%
% Body = intensity threshold + cleanup. Shaft enhancement = multiscale
% steerable ridge response with hysteresis. Saves per-frame:
%   bodyMask, res, theta, nms, scaleMap.
% Sangyoon J. Han / 2026

%% process & params
iProc = movieData.getProcessIndex('FilopodiaSegmentationProcess', 1, 0);
assert(~isempty(iProc), ['No FilopodiaSegmentationProcess found. Create it ' ...
    'and call process.run.']);
proc = movieData.processes_{iProc};
p = parseProcessParams(proc);

iChan = p.ChannelIndex;
frames = p.ProcessFrames; if isempty(frames), frames = 1:movieData.nFrames_; end

%% I/O
inFilePaths = cell(1, numel(movieData.channels_));
inFilePaths{1, iChan} = movieData.getChannelPaths{iChan};
proc.setInFilePaths(inFilePaths);

outDir = p.OutputDirectory; mkClrDir(outDir);
outFilePaths = cell(1, numel(movieData.channels_));
outFilePaths{1, iChan} = outDir;            % loadChannelOutput builds filenames
proc.setOutFilePaths(outFilePaths);

%% run
for t = frames
    img = double(movieData.channels_(iChan).loadImage(t));

    % --- cell body: blurred threshold, then despike + smooth ---
    % Filopodia are recovered separately by the steerable filter, so the body
    % must exclude them and end smooth/rounded. Opening removes the thin
    % filopodia roots; closing rounds the boundary and fills the notches
    % between them. Blur first so noise does not roughen the edge.
    imgB = imgaussfilt(img, p.GaussianBlurSigma);
    switch lower(num2str(p.BodyThreshold))
        case 'otsu',  level = thresholdOtsu(imgB);
        case 'rosin', level = thresholdRosin(imgB);
        otherwise,    level = p.BodyThreshold;   % numeric
    end
    bodyMask = imgB > level;
    bodyMask = imfill(bodyMask, 'holes');
    if any(bodyMask(:))
        bodyMask = bwareafilt(bodyMask, 1);                 % keep largest object
    end
    if p.BodyOpenRadius > 0
        bodyMask = imopen(bodyMask, strel('disk', p.BodyOpenRadius));    % despike
    end
    if p.BodyClosingRadius > 0
        bodyMask = imclose(bodyMask, strel('disk', p.BodyClosingRadius)); % round/smooth
    end
    bodyMask = imfill(bodyMask, 'holes');
    bodyMask = bwareaopen(bodyMask, p.BodyMinArea);

    % --- steerable ridge maps (on the original, unblurred image) ---
    [res, theta, nms, scaleMap] = multiscaleSteerableDetector(img, ...
        p.SteerableOrder, p.SigmaArray);

    % --- hysteresis-enhanced shaft mask (stored in nms-derived field) ---
    bg = res(~bodyMask);
    hi = p.HysteresisHigh; lo = p.HysteresisLow;
    if isempty(hi), hi = median(bg) + 3 * mad(bg, 1) * 1.4826; end
    if isempty(lo), lo = median(bg) + 1 * mad(bg, 1) * 1.4826; end
    seed = res > hi; cand = res > lo;
    shaftMask = imreconstruct(seed & cand, cand);   %#ok<NASGU> % saved with nms

    fname = fullfile(outDir, sprintf('filoSeg_frame_%04d.mat', t));
    save(fname, 'bodyMask', 'res', 'theta', 'nms', 'scaleMap', 'shaftMask');
end
disp('Filopodia segmentation done.');
end