function out = quantifyNFKB_NCratio_MD(MD, varargin)
% quantifyNFKB_NCratio_MD
% Compute NFkB nuclear/cytoplasmic intensity ratio per cell from one MovieData.
%
% INPUT
%   MD : MovieData object (single timepoint IF, z-stack)
%
% NAME-VALUE PARAMS (optional)
%   'ChNFKB'        : channel index for NFkB (default 1)
%   'ChPhalloidin'  : channel index for phalloidin (default 2)
%   'ChDAPI'        : channel index for DAPI (default 4)
%   'TimeIndex'     : time index (default 1)
%   'UseNFkBProj'   : 'mean' (default) or 'sum'
%   'DapiMIP'       : true (default)
%   'NucMinAreaPx'  : default 200
%   'NucMaxAreaPx'  : default 20000
%   'MinSolidity'   : default 0.85
%   'MaxEccentricity': default 0.97
%   'MaxCytoRadiusPx': default 30   % territory cap (10x)
%   'BgRingWidthPx' : default 15
%   'CellMaskDilatePx': default 10
%   'Verbose'       : true/false (default true)
%
% OUTPUT (struct out)
%   out.table        : table (per-cell Inuc, Icyto, NCratio, QC)
%   out.masks        : struct with nucLabel, territoryLabel, cytoLabel, cellMask
%   out.images       : struct with DAPIproj, NFKBproj, PHALproj
%   out.params       : parameter struct used
%
% NOTE
%   This method is optimized for 10x monolayer where true cell boundaries are weak.
%   It uses nuclei-seeded territory segmentation (Voronoi-like) for cytoplasm regions.

% -----------------------
% Parse parameters
% -----------------------
p = inputParser;
p.addParameter('ChNFKB', 1, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('ChPhalloidin', 2, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('ChDAPI', 4, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('TimeIndex', 1, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('UseNFkBProj', 'mean', @(s)ischar(s)||isstring(s));
p.addParameter('DapiMIP', true, @(x)islogical(x)&&isscalar(x));

p.addParameter('NucMinAreaPx', 200, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('NucMaxAreaPx', 20000, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MinSolidity', 0.85, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MaxEccentricity', 0.97, @(x)isnumeric(x)&&isscalar(x));

p.addParameter('MaxCytoRadiusPx', 30, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('BgRingWidthPx', 15, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('CellMaskDilatePx', 10, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
prm = p.Results;

% -----------------------
% Load stacks (z-stack)
% -----------------------
r = MD.getReader();
t = prm.TimeIndex;

nfkbStack = r.loadStack(prm.ChNFKB, t);       % Y X Z
phalStack = r.loadStack(prm.ChPhalloidin, t);
dapiStack = r.loadStack(prm.ChDAPI, t);

% Convert to double
nfkbStack = im2double(nfkbStack);
phalStack = im2double(phalStack);
dapiStack = im2double(dapiStack);

% -----------------------
% Z Projection
% -----------------------
if prm.DapiMIP
    D = max(dapiStack, [], 3);
else
    D = mean(dapiStack, 3);
end

switch lower(string(prm.UseNFkBProj))
    case "mean"
        R = mean(nfkbStack, 3);
    case "sum"
        R = sum(nfkbStack, 3);
    otherwise
        error("UseNFkBProj must be 'mean' or 'sum'.");
end

G = mean(phalStack, 3);

% Mild denoise
Dg = imgaussfilt(D, 1);
Rg = imgaussfilt(R, 1);
Gg = imgaussfilt(G, 2);

% -----------------------
% Nuclei segmentation (DAPI)
% -----------------------
bw = imbinarize(Dg, "adaptive", "Sensitivity", 0.45);
bw = imfill(bw, "holes");
bw = bwareaopen(bw, 50);
bw = imclearborder(bw);

% % Separate touching nuclei via watershed
% Ddist = -bwdist(~bw);
% Ddist(~bw) = -Inf;
% Lw = watershed(Ddist);
% 
% nuc = bw;
% nuc(Lw==0) = 0;
nuc = bw;                     % no watershed

nucLabel0 = bwlabel(nuc);

% Measure nucleus candidates
props = regionprops(nucLabel0, Dg, ...
    "Area","Solidity","Eccentricity","MeanIntensity","PixelIdxList");

areas = [props.Area];
sol   = [props.Solidity];
ecc   = [props.Eccentricity];
meanI = [props.MeanIntensity];

% Fragment removal heuristics:
% 1) area range
keep = areas >= prm.NucMinAreaPx & areas <= prm.NucMaxAreaPx & ...
       sol  >= prm.MinSolidity & ecc <= prm.MaxEccentricity;

% 2) remove "bright tiny outliers" among small objects
if any(keep)
    small = areas < prctile(areas, 30);  % bottom 30% area group
    if any(small)
        cutBright = prctile(meanI(small), 99); % top 1% brightness among small
        keep = keep & ~(small & meanI > cutBright);
    end
end

keepIDs = find(keep);
nucKeep = ismember(nucLabel0, keepIDs);
nucLabel = bwlabel(nucKeep);
nCells = max(nucLabel(:));

if nCells == 0
    error("No nuclei detected after filtering. Check thresholds.");
end

% -----------------------
% Build rough cellMask (optional but helps prevent territory spilling)
% Since phalloidin boundary is weak, we use low threshold + dilation
% -----------------------
th = graythresh(Gg);
cellMask = (Gg > th*0.6);  % conservative; tune if needed
cellMask = imdilate(cellMask, strel("disk", prm.CellMaskDilatePx));
cellMask = imfill(cellMask, "holes");

% If cellMask is empty or too small, fall back to "everything" minus borders
if nnz(cellMask) < 0.02*numel(cellMask)
    cellMask = true(size(cellMask));
end

% -----------------------
% Territory (Voronoi-like) segmentation
% -----------------------
seeds = nucLabel > 0;
[distToSeed, idx] = bwdist(seeds);

territoryMask = cellMask & (distToSeed <= prm.MaxCytoRadiusPx);

territoryLabel = zeros(size(nucLabel), "like", nucLabel);
territoryLabel(territoryMask) = nucLabel(idx(territoryMask));

% Cytoplasm per cell = territory - nucleus
cytoLabel = territoryLabel;
cytoLabel(nucLabel>0) = 0;

% -----------------------
% Local background subtraction + N/C ratio
% -----------------------
Inuc  = nan(nCells,1);
Icyto = nan(nCells,1);
Ibg   = nan(nCells,1);
NC    = nan(nCells,1);
nucArea = nan(nCells,1);
cytoArea = nan(nCells,1);

allNuc = nucLabel > 0;

for i = 1:nCells
    nuc_i  = (nucLabel == i);
    cyto_i = (cytoLabel == i);
    terr_i = (territoryLabel == i);

    nucArea(i)  = nnz(nuc_i);
    cytoArea(i) = nnz(cyto_i);

    if nucArea(i) < 50 || cytoArea(i) < 50
        continue
    end

    % Background ring just outside territory
    ring = imdilate(terr_i, strel("disk", prm.BgRingWidthPx)) & ~terr_i;
    ring(allNuc) = 0;  % exclude nuclei

    bgVals = Rg(ring);
    if isempty(bgVals)
        bg = median(Rg(~cellMask), "omitnan"); % fallback global bg
    else
        bg = median(bgVals, "omitnan");
    end

    In = mean(Rg(nuc_i),  "omitnan") - bg;
    Ic = mean(Rg(cyto_i), "omitnan") - bg;

    Inuc(i)  = In;
    Icyto(i) = Ic;
    Ibg(i)   = bg;

    if (Ic + eps) > 0
        NC(i) = In / (Ic + eps);
    end
end

valid = isfinite(NC) & Inuc>0 & Icyto>0;

% -----------------------
% Output table
% -----------------------
T = table((1:nCells)', nucArea, cytoArea, Inuc, Icyto, Ibg, NC, valid, ...
    'VariableNames', {'cellID','nucAreaPx','cytoAreaPx','Inuc','Icyto','Ibg','NCratio','valid'});

if prm.Verbose
    fprintf('MD: %s\n', MD.movieDataPath_);
    fprintf('Cells detected: %d | valid: %d (%.1f%%)\n', nCells, nnz(valid), 100*nnz(valid)/nCells);
    if nnz(valid) > 0
        fprintf('NCratio median (valid): %.3f\n', median(NC(valid), "omitnan"));
    end
end

% -----------------------
% Package outputs
% -----------------------
out = struct();
out.table = T;
out.images = struct('DAPIproj', D, 'NFKBproj', R, 'PHALproj', G);
out.masks  = struct('nucLabel', nucLabel, 'territoryLabel', territoryLabel, ...
                    'cytoLabel', cytoLabel, 'cellMask', cellMask);
out.params = prm;

end
