function mdList = make2DTFMFromRefBestZ_seriesMatched(movieFile, refFile, outDir, varargin)
% make2DTFMFromRefBestZ_seriesMatched
% - Handles multi-series Micro-Manager OME-TIFF (per Pos file)
% - Matches movie series <-> ref series (default: same series index)
% - For each series:
%     1) find bestZ from ref series using bead channel
%     2) extract that bestZ from movie series for all C,T
%     3) write per-channel TIFF series
%     4) create MovieData with metadata (pixelSize_, timeInterval_, emission wavelengths)
%
% Output:
%   mdList : MovieData array (one per series). If none, empty.

ip = inputParser;
ip.addParameter('beadChan', 2, @(x) isnumeric(x) && isscalar(x) && x>=1);
ip.addParameter('metric', 'tenengrad', @(s) ischar(s) || isstring(s));
ip.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
ip.addParameter('seriesMap', [], @(x) isempty(x) || isnumeric(x)); 
% seriesMap: optional mapping from movie series (0-based) to ref series (0-based).
% If empty: assumes same series index for both.

ip.addParameter('writeRefSlice', true, @(x) islogical(x) && isscalar(x));
ip.addParameter('frameDigits', 4, @(x) isnumeric(x) && isscalar(x) && x>=3);

ip.parse(varargin{:});
beadChan   = ip.Results.beadChan;
metric     = lower(string(ip.Results.metric));
verbose    = ip.Results.verbose;
seriesMap  = ip.Results.seriesMap;
writeRef   = ip.Results.writeRefSlice;
frameDigits= ip.Results.frameDigits;

movieFile = char(movieFile);
refFile   = char(refFile);
outDir    = char(outDir);
if ~exist(outDir,'dir'); mkdir(outDir); end

if verbose
    fprintf('[TFMprep] movieFile: %s\n', movieFile);
    fprintf('[TFMprep] refFile  : %s\n', refFile);
    fprintf('[TFMprep] outDir   : %s\n', outDir);
end

refReader = bfGetReader(refFile);
movReader = bfGetReader(movieFile);

nSeriesRef = refReader.getSeriesCount();
nSeriesMov = movReader.getSeriesCount();

if verbose
    fprintf('[TFMprep] SeriesCount: movie=%d, ref=%d\n', nSeriesMov, nSeriesRef);
end

% Build default series mapping (0-based)
if isempty(seriesMap)
    % Map movie series s -> ref series s (clipped to available ref series)
    seriesMap = (0:nSeriesMov-1)';
    seriesMap(seriesMap > (nSeriesRef-1)) = nSeriesRef-1;
else
    % user supplied: assume it is movieSeries->refSeries mapping (0-based)
    assert(numel(seriesMap) == nSeriesMov, ...
        'seriesMap must have length = movie series count (%d).', nSeriesMov);
end

mdCell = cell(nSeriesMov,1);

for sMov = 0:nSeriesMov-1
    sRef = seriesMap(sMov+1);

    % ---- set series for both readers ----
    refReader.setSeries(sRef);
    movReader.setSeries(sMov);

    % Use reader sizes (more reliable than metadata store in multi-series files)
    Zref = refReader.getSizeZ();
    Cref = refReader.getSizeC();
    Tref = refReader.getSizeT();

    Z = movReader.getSizeZ();
    C = movReader.getSizeC();
    T = movReader.getSizeT();

    if verbose
        fprintf('\n[TFMprep] series movie=%d (Z=%d C=%d T=%d) | ref=%d (Z=%d C=%d T=%d)\n', ...
            sMov, Z, C, T, sRef, Zref, Cref, Tref);
    end

    if Cref < beadChan
        warning('series %d: beadChan=%d exceeds ref C=%d. Skipping this series.', sMov, beadChan, Cref);
        continue;
    end
    if Zref < 1 || Z < 1 || T < 1 || C < 1
        warning('series %d: invalid dimensions. Skipping.', sMov);
        continue;
    end

    % ---- 1) Find bestZ from ref bead channel (t=0) ----
    scores = zeros(Zref,1);
    for z = 1:Zref
        plane = refReader.getIndex(z-1, beadChan-1, 0) + 1;
        I = bfGetPlane(refReader, plane);
        scores(z) = focusScore(I, metric);
    end
    [~, bestZref] = max(scores);

    % Map bestZ into movie Z if sizes differ (rare). Here: clamp.
    bestZ = min(bestZref, Z);

    if verbose
        fprintf('[TFMprep] series %d: bestZref=%d -> bestZ(movie)=%d\n', sMov, bestZref, bestZ);
    end

    % ---- Output dirs per series ----
    seriesDir = fullfile(outDir, sprintf('series_%03d', sMov));
    if ~exist(seriesDir,'dir'); mkdir(seriesDir); end

    if writeRef
        refDir = fullfile(seriesDir,'reference');
        if ~exist(refDir,'dir'); mkdir(refDir); end
        planeBest = refReader.getIndex(bestZref-1, beadChan-1, 0) + 1;
        Ibest = bfGetPlane(refReader, planeBest);
        imwrite(Ibest, fullfile(refDir, sprintf('ref_beads_bestZ_series%03d.tif', sMov)));
    end

    % ---- 2) Extract bestZ from movie series for all C,T ----
    chDirs = cell(C,1);
    for c = 1:C
        chDirs{c} = fullfile(seriesDir, sprintf('ch%02d', c));
        mkClrDir(chDirs{c});
    end

    ffmt = ['%0' num2str(frameDigits) 'd'];

    for t = 1:T
        for c = 1:C
            plane = movReader.getIndex(bestZ-1, c-1, t-1) + 1;
            I = bfGetPlane(movReader, plane);

            % enforce uint16 if needed (BioFormats sometimes returns uint8/uint16)
            % keep original type
            imwrite(I, fullfile(chDirs{c}, sprintf(['frame_' ffmt '.tif'], t)));
        end
        if verbose && (mod(t, max(1,round(T/10)))==0)
            fprintf('  [series %d] wrote frame %d/%d\n', sMov, t, T);
        end
    end

    % ---- 3) Create MovieData with metadata ----
    channels = Channel.empty(0,C);
    for c = 1:C
        channels(c) = Channel(chDirs{c});  % folder with frame_####.tif
        % (optional) add channel name if you want: channels(c).setChannelName(...)
    end

    md = MovieData(channels, seriesDir, ...
        'movieDataPath_', seriesDir, ...
        'movieDataFileName_', 'movieData.mat');

    % Populate metadata from OME store if possible
    try
        % pixel size (µm/pixel) -> MovieData.pixelSize_ wants nm/pixel in u-track
        pixSizeNm = getPixelSizeNmFromStore(movReader.getMetadataStore(), sMov);
        if ~isempty(pixSizeNm) && isprop(md,'pixelSize_')
            md.pixelSize_ = pixSizeNm;
        end
    catch ME
        if verbose, warning('series %d: pixel size read failed: %s', sMov, ME.message); end
    end

    try
        dtSec = getTimeIntervalSecFromStore(movReader.getMetadataStore(), sMov);
        if ~isempty(dtSec) && isprop(md,'timeInterval_')
            md.timeInterval_ = dtSec;
        end
    catch ME
        if verbose, warning('series %d: time interval read failed: %s', sMov, ME.message); end
    end

    try
        emNm = getEmissionNmFromStore(movReader.getMetadataStore(), sMov, C);
        if ~isempty(emNm)
            for c = 1:min(C, numel(emNm))
                if isprop(md.channels_{c}, 'emissionWavelength_')
                    md.channels_{c}.emissionWavelength_ = emNm(c);
                end
            end
        end
    catch ME
        if verbose, warning('series %d: emission wavelength read failed: %s', sMov, ME.message); end
    end

    % Make sure MD infers nFrames_ etc from files
    md.sanityCheck();
    md.save();

    if verbose
        fprintf('[TFMprep] Saved MovieData: %s\n', fullfile(seriesDir,'movieData.mat'));
    end

    mdCell{sMov+1} = md;
end

% close readers
try refReader.close(); catch, end
try movReader.close(); catch, end

% return as array
mdList = [mdCell{~cellfun(@isempty, mdCell)}];

end

%% -------- Helpers --------

function s = focusScore(I, metric)
I = double(I);
switch metric
    case "variance"
        s = var(I(:));
    case "tenengrad"
        [Gx,Gy] = imgradientxy(I);
        s = mean(Gx(:).^2 + Gy(:).^2);
    otherwise
        error('Unknown metric: %s', metric);
end
end

function pixNm = getPixelSizeNmFromStore(store, series0)
% Tries PixelsPhysicalSizeX in µm -> convert to nm
pixNm = [];
try
    % OME store uses image index; often equals series index, but not guaranteed.
    % We attempt series0; if missing, try 0.
    idx = series0;
    v = store.getPixelsPhysicalSizeX(idx);
    if isempty(v); v = store.getPixelsPhysicalSizeX(0); idx = 0; end
    if ~isempty(v)
        % v is Length object, typically in micrometers
        pixUm = v.value().doubleValue(); % µm
        pixNm = pixUm * 1000;            % nm
    end
catch
    % leave empty
end
end

function dtSec = getTimeIntervalSecFromStore(store, series0)
dtSec = [];
try
    idx = series0;
    % Try time increment (global)
    v = store.getPixelsTimeIncrement(idx);
    if isempty(v); v = store.getPixelsTimeIncrement(0); idx = 0; end
    if ~isempty(v)
        % v is Time object, usually seconds
        dtSec = v.value().doubleValue();
        return;
    end
catch
end

% Fallback: attempt delta between PlaneDeltaT(0) and PlaneDeltaT(1) if present
try
    idx = series0;
    t0 = store.getPlaneDeltaT(idx, 0);
    t1 = store.getPlaneDeltaT(idx, 1);
    if ~isempty(t0) && ~isempty(t1)
        dtSec = t1.value().doubleValue() - t0.value().doubleValue();
    end
catch
end
end

function emNm = getEmissionNmFromStore(store, series0, C)
emNm = nan(1,C);
try
    idx = series0;
    % OME: Channel index per Pixels; sometimes store.getChannelEmissionWavelength(image, channel)
    for c0 = 0:C-1
        w = store.getChannelEmissionWavelength(idx, c0);
        if isempty(w)
            w = store.getChannelEmissionWavelength(0, c0);
        end
        if ~isempty(w)
            emNm(c0+1) = w.value().doubleValue(); % nm
        end
    end
catch
    emNm = [];
end
if all(isnan(emNm)), emNm = []; end
end