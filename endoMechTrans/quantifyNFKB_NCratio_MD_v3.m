function out = quantifyNFKB_NCratio_MD_v3(md, varargin)
% quantifyNFKB_NCratio_MD_v3
% Robust nuclei detection for z-stack DAPI in IF: top-K projection + adaptive + LoG seed rescue.
% Cytoplasm remains nucleus-centered annulus (ring) to avoid boundary issues at 10x.
%
% OUTPUT same structure as v2.

% -----------------------
% Parameters
% -----------------------
p = inputParser;

% Channels
p.addParameter('ChNFKB', 1, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('ChPhalloidin', 2, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('ChDAPI', 4, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('TimeIndex', 1, @(x)isnumeric(x)&&isscalar(x));

% Projections
p.addParameter('UseNFkBProj', 'mean', @(s)ischar(s)||isstring(s)); % 'mean'|'sum'
p.addParameter('DapiProjection', 'topkmean', @(s)ischar(s)||isstring(s)); % 'mip'|'mean'|'topkmean'
p.addParameter('TopK', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);

% DAPI segmentation parameters
p.addParameter('AdaptiveSensitivity', 0.62, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
p.addParameter('PreSmoothSigma', 1.0, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MinObjectPx', 50, @(x)isnumeric(x)&&isscalar(x));

% Nuclei filtering
p.addParameter('NucMinAreaPx', 180, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('NucMaxAreaPx', 30000, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MinSolidity', 0.80, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MaxEccentricity', 0.995, @(x)isnumeric(x)&&isscalar(x));

% LoG rescue options
p.addParameter('UseLoGRescue', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('LoGSigma', 2.0, @(x)isnumeric(x)&&isscalar(x));        % tune 1.5~4 for 10x
p.addParameter('LoGSeedPrct', 85, @(x)isnumeric(x)&&isscalar(x));      % keep top seeds
p.addParameter('SeedDilatePx', 3, @(x)isnumeric(x)&&isscalar(x));      % seed -> small disk

% Cytoplasm ring geometry
p.addParameter('CytoInnerRadiusPx', 6, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('CytoOuterRadiusPx', 28, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('BgOuterRadiusPx', 40, @(x)isnumeric(x)&&isscalar(x));

% CellMask
p.addParameter('UseCellMask', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('CellMaskDilatePx', 10, @(x)isnumeric(x)&&isscalar(x));

% Misc
p.addParameter('ClearBorder', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));

p.parse(varargin{:});
prm = p.Results;

% -----------------------
% Load stacks
% -----------------------
r = md.getReader();
t = prm.TimeIndex;

nfkbStack = im2double(r.loadStack(prm.ChNFKB, t));
phalStack = im2double(r.loadStack(prm.ChPhalloidin, t));
dapiStack = im2double(r.loadStack(prm.ChDAPI, t));

% -----------------------
% Projections
% -----------------------
% DAPI projection
switch lower(string(prm.DapiProjection))
    case "mip"
        D = max(dapiStack, [], 3);
    case "mean"
        D = mean(dapiStack, 3);
    case "topkmean"
        % top-K mean per pixel (robust: reduces haze vs sum, avoids MIP hot pixels)
        Ds = sort(dapiStack, 3, 'descend');
        K = min(prm.TopK, size(Ds,3));
        D = mean(Ds(:,:,1:K), 3);
    otherwise
        error("DapiProjection must be 'mip','mean','topkmean'.");
end

% NFkB projection
switch lower(string(prm.UseNFkBProj))
    case "mean"
        R = mean(nfkbStack, 3);
    case "sum"
        R = sum(nfkbStack, 3);
    otherwise
        error("UseNFkBProj must be 'mean' or 'sum'.");
end

G = mean(phalStack, 3);

% Analysis images (denoise)
Dg = imgaussfilt(D, prm.PreSmoothSigma);
Rg = imgaussfilt(R, 1);
Gg = imgaussfilt(G, 2);

% -----------------------
% 1) Initial nuclei mask via adaptive threshold
% -----------------------
bw = imbinarize(Dg, "adaptive", "Sensitivity", prm.AdaptiveSensitivity);
bw = imfill(bw, "holes");
bw = bwareaopen(bw, prm.MinObjectPx);
if prm.ClearBorder
    bw = imclearborder(bw);
end

% -----------------------
% 2) LoG seed rescue to reduce missed nuclei
% -----------------------
if prm.UseLoGRescue
    sigma = prm.LoGSigma;
    h = fspecial('log', ceil(sigma*6), sigma);
    LoG = imfilter(Dg, h, 'replicate');
    LoG = -LoG; % blobs -> positive peaks

    seedCand = imregionalmax(LoG);

    % keep strong maxima only
    vals = LoG(seedCand);
    if ~isempty(vals)
        thr = prctile(vals, prm.LoGSeedPrct);
        seedCand = seedCand & (LoG >= thr);
    end

    % remove seeds already inside bw
    seedCand(bw) = 0;

    % convert seeds to small disks and OR with bw
    seedDisk = imdilate(seedCand, strel('disk', prm.SeedDilatePx));
    bw = bw | seedDisk;

    % clean again
    bw = imfill(bw, "holes");
    bw = bwareaopen(bw, prm.MinObjectPx);
end

% -----------------------
% 3) Label + filter nuclei candidates (remove fragments but avoid over-pruning)
% -----------------------
nucLabel0 = bwlabel(bw);

props = regionprops(nucLabel0, Dg, ...
    "Area","Solidity","Eccentricity","MeanIntensity");

areas = [props.Area];
sol   = [props.Solidity];
ecc   = [props.Eccentricity];
meanI = [props.MeanIntensity];

keep = areas >= prm.NucMinAreaPx & areas <= prm.NucMaxAreaPx & ...
       sol  >= prm.MinSolidity & ecc <= prm.MaxEccentricity;

% Bright-small fragment suppression (more conservative than v2)
% Only apply to VERY small objects (bottom 15% area)
if any(keep)
    small = areas < prctile(areas, 15);
    if any(small)
        cutBright = prctile(meanI(small), 99.5);
        keep = keep & ~(small & meanI > cutBright);
    end
end

keepIDs = find(keep);
nucKeep = ismember(nucLabel0, keepIDs);
nucLabel = bwlabel(nucKeep);

nCells = max(nucLabel(:));
if nCells == 0
    error("No nuclei detected after filtering. Tune thresholds.");
end

% -----------------------
% 4) CellMask (conservative)
% -----------------------
if prm.UseCellMask
    th = graythresh(Gg);
    cellMask = (Gg > th*0.6);
    cellMask = imdilate(cellMask, strel("disk", prm.CellMaskDilatePx));
    cellMask = imfill(cellMask, "holes");
    if nnz(cellMask) < 0.02*numel(cellMask)
        cellMask = true(size(cellMask));
    end
else
    cellMask = true(size(D));
end

% -----------------------
% 5) Cytoplasm annulus + bg ring per cell (same as v2)
% -----------------------
cytoLabel = zeros(size(nucLabel), "like", nucLabel);
bgRingLabel = zeros(size(nucLabel), "like", nucLabel);

allNuc = nucLabel > 0;

seIn  = strel("disk", prm.CytoInnerRadiusPx);
seOut = strel("disk", prm.CytoOuterRadiusPx);
seBg  = strel("disk", prm.BgOuterRadiusPx);

allDilOut = imdilate(allNuc, seOut);

for i = 1:nCells
    thisNuc = (nucLabel == i);

    inner = imdilate(thisNuc, seIn);
    outer = imdilate(thisNuc, seOut);
    ring  = outer & ~inner;

    ring = ring & cellMask;
    ring(allNuc) = 0;

    otherDil = allDilOut & ~imdilate(thisNuc, seOut);
    ring(otherDil) = 0;

    cytoLabel(ring) = i;

    bgOuter = imdilate(thisNuc, seBg);
    bgRing = bgOuter & ~outer;
    bgRing = bgRing & cellMask;
    bgRing(allNuc) = 0;
    bgRing(otherDil) = 0;

    bgRingLabel(bgRing) = i;
end

% -----------------------
% 6) Intensities + local background subtraction
% -----------------------
Inuc  = nan(nCells,1);
Icyto = nan(nCells,1);
Ibg   = nan(nCells,1);
NC    = nan(nCells,1);
nucArea = nan(nCells,1);
cytoArea = nan(nCells,1);
bgArea = nan(nCells,1);

for i = 1:nCells
    nuc_i  = (nucLabel == i);
    cyto_i = (cytoLabel == i);
    bg_i   = (bgRingLabel == i);

    nucArea(i)  = nnz(nuc_i);
    cytoArea(i) = nnz(cyto_i);
    bgArea(i)   = nnz(bg_i);

    if nucArea(i) < 40 || cytoArea(i) < 40 || bgArea(i) < 40
        continue
    end

    bg = median(Rg(bg_i), "omitnan");
    In = mean(Rg(nuc_i),  "omitnan") - bg;
    Ic = mean(Rg(cyto_i), "omitnan") - bg;

    Inuc(i)  = In;
    Icyto(i) = Ic;
    Ibg(i)   = bg;
    NC(i)    = In / (Ic + eps);
end

valid = isfinite(NC) & Inuc>0 & Icyto>0;

T = table((1:nCells)', nucArea, cytoArea, bgArea, Inuc, Icyto, Ibg, NC, valid, ...
    'VariableNames', {'cellID','nucAreaPx','cytoAreaPx','bgAreaPx','Inuc','Icyto','Ibg','NCratio','valid'});

if prm.Verbose
    fprintf('MD: %s\n', md.movieDataPath_);
    fprintf('DAPI proj: %s (TopK=%d) | AdaptiveSensitivity=%.2f | LoGRescue=%d (sigma=%.2f)\n', ...
        string(prm.DapiProjection), prm.TopK, prm.AdaptiveSensitivity, prm.UseLoGRescue, prm.LoGSigma);
    fprintf('Cells: %d | valid: %d (%.1f%%)\n', nCells, nnz(valid), 100*nnz(valid)/nCells);
end

out = struct();
out.table = T;
out.images = struct('DAPIproj', D, 'NFKBproj', R, 'PHALproj', G, 'LoGmap', []);
out.masks  = struct('nucLabel', nucLabel, 'cytoLabel', cytoLabel, ...
                    'bgRingLabel', bgRingLabel, 'cellMask', cellMask);
out.params = prm;

end