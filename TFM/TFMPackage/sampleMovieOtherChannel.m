
function sampleMovieOtherChannel(movieData, varargin)
% sampleMovieOtherChannel
% Sample intensity statistics from a specified channel inside a binary mask
% for each frame, optionally computing dF/F0 and per-object stats.
%
% INPUT
%   movieData : MovieData
%   paramsIn  : struct (optional)
%
% OUTPUT (saved to funParams.OutputDirectory)
%   otherChannelSampling.mat with struct 'S'
%
% Notes
% - Mask is expected to be a binary image per frame (logical or 0/1).
% - If StageDriftCorrectionProcess exists in TFMPackage and
%   UseStageDriftCorrection=true, masks are warped per-frame using T shifts.
%
% 2026 - Sangyoon Han lab

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn', [], @(x) isempty(x) || isstruct(x));
ip.parse(movieData, varargin{:});
paramsIn = ip.Results.paramsIn;

% Ensure process exists
iProc = movieData.getProcessIndex('OtherChannelSamplingProcess', 1, 0);
if isempty(iProc)
    iProc = numel(movieData.processes_) + 1;
    movieData.addProcess(OtherChannelSamplingProcess(movieData, ...
        fullfile(movieData.outputDirectory_, 'otherChannelSampling')));
end
proc = movieData.processes_{iProc};

% Parse params
p = parseProcessParams(proc, paramsIn);

% Output dir
if ~exist(p.OutputDirectory,'dir'); mkdir(p.OutputDirectory); end
outMat = fullfile(p.OutputDirectory, 'otherChannelSampling.mat');
proc.setOutFilePaths({outMat});

nFrames = movieData.nFrames_;
ch = p.ChannelIndex;

% --- Resolve mask process ---
maskProc = [];
if ~isempty(p.MaskProcessName)
    iMask = movieData.getProcessIndex(p.MaskProcessName, 1, 0);
    if ~isempty(iMask)
        maskProc = movieData.getProcess(iMask);
    end
else
    % common candidates (in priority)
    cand = {'MaskIntersectionProcess','MaskRefinementProcess','ThresholdProcess'};
    for i = 1:numel(cand)
        iMask = movieData.getProcessIndex(cand{i}, 1, 0);
        if ~isempty(iMask)
            maskProc = movieData.getProcess(iMask);
            p.MaskProcessName = cand{i};
            break;
        end
    end
end

if isempty(maskProc)
    error('OtherChannelSamplingProcess:NoMask', ...
        'Could not find a mask process. Set params.MaskProcessName or run a mask/threshold process first.');
end

maskChan = p.MaskChannelIndex;
if isempty(maskChan) || ~isscalar(maskChan) || maskChan < 1
    maskChan = 1;
end

% --- Optional: Stage drift correction shifts (T) ---
useSDC = isfield(p,'UseStageDriftCorrection') && p.UseStageDriftCorrection;
T = [];
existSDC = false;
if useSDC
    try
        iPack = movieData.getPackageIndex('TFMPackage');
        if ~isempty(iPack)
            pack = movieData.getPackage(iPack);
            sdcProc = pack.processes_{1};
            if ~isempty(sdcProc) && isprop(sdcProc,'outFilePaths_')
                % bead channel for SDC output:
                try
                    iBeadChan = sdcProc.funParams_.iBeadChannel;
                catch
                    iBeadChan = 1;
                end
                if size(sdcProc.outFilePaths_,1) >= 3 && size(sdcProc.outFilePaths_,2) >= iBeadChan
                    s = load(sdcProc.outFilePaths_{3,iBeadChan}, 'T');
                    T = s.T;
                    existSDC = true;
                end
            end
        end
    catch
        existSDC = false;
    end
end

ref_obj = [];
if existSDC
    ref_obj = imref2d(movieData.imSize_);
end

% --- Preallocate output ---
S = struct();
S.params = p;
S.frame = (1:nFrames)';
S.mean = nan(nFrames,1);
S.median = nan(nFrames,1);
S.areaPix = nan(nFrames,1);
S.nObjects = nan(nFrames,1);
S.objectMeans = cell(nFrames,1);
S.objectAreas = cell(nFrames,1);
S.objectCentroids = cell(nFrames,1);

if p.SavePerFrameTifPreview
    previewDir = fullfile(p.OutputDirectory, 'preview');
    if ~exist(previewDir,'dir'); mkdir(previewDir); end
end

% --- Main loop ---
for ii = 1:nFrames

    % Load image frame from desired channel
    I = movieData.channels_(ch).loadImage(ii);
    I = double(I);

    % Load mask (binary)
    mask = [];
    try
        % Many mask processes use: loadChannelOutput(iChan, iFrame)
        mask = maskProc.loadChannelOutput(maskChan, ii);
    catch
        try
            mask = maskProc.loadChannelOutput(ii);
        catch
            % last resort: try reading from masks folder if it exists
            mask = [];
        end
    end
    if isempty(mask)
        error('OtherChannelSamplingProcess:MaskLoadFailed', ...
            'Failed to load mask for frame %d from %s.', ii, class(maskProc));
    end
    mask = logical(mask);

    % Apply stage drift correction to mask if requested and available
    if existSDC && size(T,1) >= ii
        Tr = affine2d([1 0 0; 0 1 0; fliplr(T(ii,:)) 1]);
        mask = imwarp(mask, Tr, 'OutputView', ref_obj);
        mask = mask > 0;
    end

    % Summary stats within mask
    pix = I(mask);
    if isempty(pix)
        S.mean(ii) = NaN;
        S.median(ii) = NaN;
        S.areaPix(ii) = 0;
        S.nObjects(ii) = 0;
        continue;
    end

    S.mean(ii) = mean(pix, 'omitnan');
    S.median(ii) = median(pix, 'omitnan');
    S.areaPix(ii) = nnz(mask);

    % Per-object stats (connected components on mask)
    if isfield(p,'UseLabeling') && p.UseLabeling
        CC = bwconncomp(mask);
        if CC.NumObjects == 0
            S.nObjects(ii) = 0;
        else
            % Filter by area
            areas = cellfun(@numel, CC.PixelIdxList);
            keep = areas >= p.MinAreaPix;
            CC.PixelIdxList = CC.PixelIdxList(keep);
            areas = areas(keep);
            CC.NumObjects = numel(CC.PixelIdxList);

            S.nObjects(ii) = CC.NumObjects;
            if CC.NumObjects > 0
                objMeans = nan(CC.NumObjects,1);
                cent = nan(CC.NumObjects,2);
                for k = 1:CC.NumObjects
                    idx = CC.PixelIdxList{k};
                    objMeans(k) = mean(I(idx), 'omitnan');
                end
                % centroids using regionprops on label matrix
                L = labelmatrix(CC);
                rp = regionprops(L, 'Centroid');
                for k = 1:numel(rp)
                    cent(k,:) = rp(k).Centroid;
                end

                S.objectMeans{ii} = objMeans;
                S.objectAreas{ii} = areas(:);
                S.objectCentroids{ii} = cent;
            end
        end
    end

    % Optional preview
    if p.SavePerFrameTifPreview
        rgb = repmat(mat2gray(I),1,1,3);
        b = bwperim(mask);
        rgb(:,:,1) = max(rgb(:,:,1), b);
        imwrite(rgb, fullfile(previewDir, sprintf('frame_%03d.png', ii)));
    end
end

% --- dF/F0 ---
if isfield(p,'ComputeDFF0') && p.ComputeDFF0
    baseFrames = p.BaselineFrames;
    baseFrames = baseFrames(baseFrames>=1 & baseFrames<=nFrames);
    if isempty(baseFrames)
        baseFrames = 1;
    end
    F0 = mean(S.mean(baseFrames), 'omitnan');
    if isempty(F0) || ~isfinite(F0) || F0==0
        S.F0 = NaN;
        S.dFF0 = nan(nFrames,1);
    else
        S.F0 = F0;
        S.dFF0 = (S.mean - F0) ./ F0;
    end
end

save(outMat, 'S');

movieData.save();
end
