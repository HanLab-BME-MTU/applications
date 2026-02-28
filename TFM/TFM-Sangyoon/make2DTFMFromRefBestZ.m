function mdList = make2DTFMFromRefBestZ(movieFile, refFile, outRoot, varargin)
% make2DTFMFromRefBestZ
% - Use ONE movie OME-TIFF that contains ALL positions as series
% - Use ONE ref OME-TIFF that contains ALL positions as series
% - For each series (position):
%    1) setSeries(s) on both readers
%    2) find bestZ from ref bead channel (t=0)
%    3) extract bestZ from movie for all C,T to per-channel TIFF folders
%    4) create MovieData, populate pixelSize_/timeInterval_/emission wavelengths
%
% Output:
%   mdList: MovieData array (one per position/series)

ip = inputParser;
ip.addParameter('beadChan', 2, @(x) isnumeric(x) && isscalar(x) && x>=1);
ip.addParameter('metric', 'tenengrad', @(s) ischar(s) || isstring(s));
ip.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
ip.addParameter('frameDigits', 4, @(x) isnumeric(x) && isscalar(x) && x>=3);
ip.addParameter('writeRefSlice', true, @(x) islogical(x) && isscalar(x));
ip.parse(varargin{:});

beadChan   = ip.Results.beadChan;
metric     = lower(string(ip.Results.metric));
verbose    = ip.Results.verbose;
frameDigits= ip.Results.frameDigits;
writeRef   = ip.Results.writeRefSlice;

movieFile = char(movieFile);
refFile   = char(refFile);
outRoot   = char(outRoot);
if ~exist(outRoot,'dir'); mkdir(outRoot); end

refR = bfGetReader(refFile);
movR = bfGetReader(movieFile);

nSref = refR.getSeriesCount();
nSmov = movR.getSeriesCount();

if verbose
    fprintf('[TFMprep] movie series=%d | ref series=%d\n', nSmov, nSref);
end

nS = min(nSmov, nSref);
stMov = movR.getMetadataStore();

mdCell = cell(nS,1);
ffmt = ['%0' num2str(frameDigits) 'd'];

for s = 0:nS-1
    refR.setSeries(s);
    movR.setSeries(s);

    % ImageName used as position label
    posName = sprintf('series_%03d', s);
    try
        nm = stMov.getImageName(s);
        if ~isempty(nm)
            posName = char(nm);
        end
    catch
    end
    posNameClean = regexprep(posName,'[^\w\-\.\%]+','_'); % safe folder name

    Zref = refR.getSizeZ(); Cref = refR.getSizeC(); Tref = refR.getSizeT();
    Z    = movR.getSizeZ(); C    = movR.getSizeC(); T    = movR.getSizeT();

    if verbose
        fprintf('\n[TFMprep] %s | series=%d | movie(Z=%d C=%d T=%d) ref(Z=%d C=%d T=%d)\n', ...
            posName, s, Z, C, T, Zref, Cref, Tref);
    end

    if Cref < beadChan
        warning('%s: beadChan=%d exceeds ref C=%d. Skipping.', posName, beadChan, Cref);
        continue;
    end

    % ---- 1) bestZ from ref bead channel ----
    scores = zeros(Zref,1);
    for z = 1:Zref
        plane = refR.getIndex(z-1, beadChan-1, 0) + 1;
        I = bfGetPlane(refR, plane);
        scores(z) = focusScore(I, metric);
    end
    [~, bestZref] = max(scores);
    bestZ = min(bestZref, Z);

    if verbose
        fprintf('[TFMprep] %s: bestZref=%d -> bestZ(movie)=%d\n', posName, bestZref, bestZ);
    end

    % ---- 2) output dirs ----
    outDir = fullfile(outRoot, posNameClean);
    if ~exist(outDir,'dir'); mkdir(outDir); end

    if writeRef
        refDir = fullfile(outDir,'reference');
        if ~exist(refDir,'dir'); mkdir(refDir); end
        planeBest = refR.getIndex(bestZref-1, beadChan-1, 0) + 1;
        Ibest = bfGetPlane(refR, planeBest);
        imwrite(Ibest, fullfile(refDir, 'ref_beads_bestZ.tif'));
    end

    % ---- 3) extract bestZ for all C,T ----
    chDirs = cell(C,1);
    for c = 1:C
        chDirs{c} = fullfile(outDir, sprintf('ch%02d', c));
        mkClrDir(chDirs{c});
    end

    for t = 1:T
        for c = 1:C
            plane = movR.getIndex(bestZ-1, c-1, t-1) + 1;
            I = bfGetPlane(movR, plane);
            imwrite(I, fullfile(chDirs{c}, sprintf(['frame_' ffmt '.tif'], t)));
        end
    end

    % ---- 4) create MovieData + metadata ----
    channels = Channel.empty(0,C);
    for c = 1:C
        channels(c) = Channel(chDirs{c});
    end

    md = MovieData(channels, outDir, 'movieDataPath_', outDir, 'movieDataFileName_', 'movieData.mat');

    % pixelSize_, timeInterval_, emission
    try
        pixNm = getPixelSizeNmFromStore(stMov, s);
        if ~isempty(pixNm) && isprop(md,'pixelSize_'), md.pixelSize_ = pixNm; end
    catch, end

    try
        dtSec = getTimeIntervalSecFromStore(stMov, s);
        if ~isempty(dtSec) && isprop(md,'timeInterval_'), md.timeInterval_ = dtSec; end
    catch, end

    try
        emNm = getEmissionNmFromStore(stMov, s, C);
        if ~isempty(emNm)
            for c = 1:min(C,numel(emNm))
                if isprop(md.channels_{c}, 'emissionWavelength_')
                    md.channels_{c}.emissionWavelength_ = emNm(c);
                end
            end
        end
    catch, end

    md.sanityCheck();
    md.save();

    mdCell{s+1} = md;
end

try refR.close(); catch, end
try movR.close(); catch, end

mdList = [mdCell{~cellfun(@isempty, mdCell)}];
end

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

function pixNm = getPixelSizeNmFromStore(store, imgIdx)
pixNm = [];
try
    v = store.getPixelsPhysicalSizeX(imgIdx);
    if isempty(v), v = store.getPixelsPhysicalSizeX(0); end
    if ~isempty(v)
        pixUm = v.value().doubleValue();
        pixNm = pixUm * 1000;
    end
catch
end
end

function dtSec = getTimeIntervalSecFromStore(store, imgIdx)
dtSec = [];
try
    v = store.getPixelsTimeIncrement(imgIdx);
    if isempty(v), v = store.getPixelsTimeIncrement(0); end
    if ~isempty(v)
        dtSec = v.value().doubleValue();
        return;
    end
catch
end
try
    t0 = store.getPlaneDeltaT(imgIdx, 0);
    t1 = store.getPlaneDeltaT(imgIdx, 1);
    if ~isempty(t0) && ~isempty(t1)
        dtSec = t1.value().doubleValue() - t0.value().doubleValue();
    end
catch
end
end

function emNm = getEmissionNmFromStore(store, imgIdx, C)
emNm = nan(1,C);
try
    for c0 = 0:C-1
        w = store.getChannelEmissionWavelength(imgIdx, c0);
        if isempty(w), w = store.getChannelEmissionWavelength(0, c0); end
        if ~isempty(w)
            emNm(c0+1) = w.value().doubleValue();
        end
    end
catch
    emNm = [];
end
if all(isnan(emNm)), emNm = []; end
end