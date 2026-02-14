function md2d = make2DTFMFromRefBestZ(movieFile, refFile, outDir, varargin)
% make2DTFMFromRefBestZ (bfGetReader-only, robust for Micro-Manager OME-TIFF)
% - Find best-focus Z from reference using bead channel
% - Extract that Z from the movie for all frames & channels
% - Save per-channel TIFF series and create MovieData

ip = inputParser;
ip.addParameter('beadChan', 2, @(x) isnumeric(x) && isscalar(x) && x>=1);
ip.addParameter('metric', 'tenengrad', @(s) ischar(s) || isstring(s));
ip.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
ip.parse(varargin{:});
beadChan = ip.Results.beadChan;
metric   = lower(string(ip.Results.metric));
verbose  = ip.Results.verbose;

movieFile = char(movieFile);
refFile   = char(refFile);
outDir    = char(outDir);
if ~exist(outDir,'dir'); mkdir(outDir); end

%% ---- 1) Reference: best-focus Z from bead channel ----
refReader = bfGetReader(refFile);
refMeta = refReader.getMetadataStore();

Zref = refMeta.getPixelsSizeZ(0).getValue();
Cref = refMeta.getPixelsSizeC(0).getValue();
Tref = refMeta.getPixelsSizeT(0).getValue();

if verbose
    fprintf('[TFMprep] ref dims: Z=%d, C=%d, T=%d\n', Zref, Cref, Tref);
end

assert(beadChan <= Cref, 'beadChan=%d exceeds reference channels C=%d.', beadChan, Cref);

scores = zeros(Zref,1);
for z = 1:Zref
    plane = refReader.getIndex(z-1, beadChan-1, 0) + 1; % z,c,t (0-based) -> 1-based plane
    I = bfGetPlane(refReader, plane);
    scores(z) = focusScore(I, metric);
end
[~, bestZ] = max(scores);

refDir = fullfile(outDir,'reference');
if ~exist(refDir,'dir'); mkdir(refDir); end
planeBest = refReader.getIndex(bestZ-1, beadChan-1, 0) + 1;
Ibest = bfGetPlane(refReader, planeBest);
imwrite(Ibest, fullfile(refDir,'ref_beads_bestZ.tif'));

refReader.close();

if verbose
    fprintf('[TFMprep] bestZ=%d (metric=%s)\n', bestZ, metric);
end

%% ---- 2) Movie: extract bestZ for all C,T ----
movReader = bfGetReader(movieFile);
movMeta = movReader.getMetadataStore();

Z = movMeta.getPixelsSizeZ(0).getValue();
C = movMeta.getPixelsSizeC(0).getValue();
T = movMeta.getPixelsSizeT(0).getValue();

if verbose
    fprintf('[TFMprep] movie dims: Z=%d, C=%d, T=%d\n', Z, C, T);
end

assert(bestZ <= Z, 'bestZ=%d exceeds movie SizeZ=%d.', bestZ, Z);

% Make channel dirs
chDirs = cell(C,1);
for c = 1:C
    chDirs{c} = fullfile(outDir, sprintf('ch%02d', c));
    mkClrDir(chDirs{c});
end

for t = 1:T
    for c = 1:C
        plane = movReader.getIndex(bestZ-1, c-1, t-1) + 1;
        I = bfGetPlane(movReader, plane);
        imwrite(I, fullfile(chDirs{c}, sprintf('frame_%03d.tif', t)));
    end
    if verbose && mod(t,10)==0
        fprintf('  wrote frame %d/%d\n', t, T);
    end
end
movReader.close();

%% ---- 3) Create MovieData ----
channels(C) = Channel();
for c = 1:C
    channels(c) = Channel(chDirs{c});
end
md2d = MovieData(channels, outDir, 'movieDataPath_', outDir, 'movieDataFileName_', 'movieData.mat');

% Do NOT set read-only properties (nFrames_, zSize_) manually.
% Let MovieData infer them from the channel image files.
md2d.reader = [];   % ok if this one is settable; if not, remove this line too

md2d.sanityCheck;   % this should populate nFrames_ based on frame_###.tif count
md2d.save;

if verbose
    fprintf('[TFMprep] Created md2d: %s\n', fullfile(outDir,'movieData.mat'));
end

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