function out = quantifyNFKB_NCratio_MD_v2(md, varargin)
% quantifyNFKB_NCratio_MD_v2
% One MovieData -> per-cell NFkB N/C ratio (IF single timepoint z-stack)
%
% Key design (for 10x, weak cell boundary):
% - DAPI MIP nuclei segmentation (NO watershed -> avoid split nuclei)
% - Cytoplasm defined as nucleus-centered annulus ("ring") with overlap exclusion
% - Local background subtraction from an outer ring around the cytoplasm annulus
%
% Channels (default):
%   Ch1 NFkB, Ch2 phalloidin, Ch4 DAPI
%
% Outputs:
%   out.table: cellID, Inuc, Icyto, NCratio, QC fields
%   out.masks: nucLabel, cytoLabel (ring), bgRingLabel, cellMask
%   out.images: projections

% -----------------------
% Parameters
% -----------------------
p = inputParser;
p.addParameter('ChNFKB', 1, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('ChPhalloidin', 2, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('ChDAPI', 4, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('TimeIndex', 1, @(x)isnumeric(x)&&isscalar(x));

p.addParameter('UseNFkBProj', 'mean', @(s)ischar(s)||isstring(s)); % 'mean' or 'sum'
p.addParameter('DapiMIP', true, @(x)islogical(x)&&isscalar(x));

% Nuclei filtering (tune once)
p.addParameter('NucMinAreaPx', 120, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('NucMaxAreaPx', 20000, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MinSolidity', 0.85, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MaxEccentricity', 0.98, @(x)isnumeric(x)&&isscalar(x));

% Cytoplasm ring geometry (10x)
p.addParameter('CytoInnerRadiusPx', 6, @(x)isnumeric(x)&&isscalar(x));   % gap away from nucleus
p.addParameter('CytoOuterRadiusPx', 28, @(x)isnumeric(x)&&isscalar(x));  % ring thickness/extent
% Background ring outside cytoplasm ring
p.addParameter('BgOuterRadiusPx', 40, @(x)isnumeric(x)&&isscalar(x));    % bg ring outer radius

% CellMask (optional; helps avoid ring spilling into empty field)
p.addParameter('UseCellMask', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('CellMaskDilatePx', 10, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('DapiProjMode', 'trimmedMean', @(s)ischar(s)||isstring(s)); 
p.addParameter('DapiTrimTopFrac', 0.02, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<1);

% Display/QC options
p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
prm = p.Results;

% -----------------------
% Load stacks (z-stack)
% -----------------------
r = md.getReader();
t = prm.TimeIndex;

nfkbStack = im2double(r.loadStack(prm.ChNFKB, t));
phalStack = im2double(r.loadStack(prm.ChPhalloidin, t));
dapiStack = im2double(r.loadStack(prm.ChDAPI, t));

% -----------------------
% Projections
% -----------------------
% --- DAPI projection (replace MIP/mean with top-K mean for robustness)
% -----------------------
% Projections (DAPI)
% -----------------------
switch lower(string(prm.DapiProjMode))
    case "topk"
        K = 3;
        Ds = sort(dapiStack, 3, 'descend');
        Rs = sort(nfkbStack, 3, 'descend');
        Gs = sort(phalStack, 3, 'descend');
        K = min(K, size(Ds,3));
        D  = mean(Ds(:,:,1:K), 3);
        R  = mean(Rs(:,:,1:K), 3);
        G  = mean(Gs(:,:,1:K), 3);
        
    case "mip"
        D = max(dapiStack, [], 3);
        R = max(nfkbStack, [], 3);
        G = max(phalStack, [], 3);

    case "mean"
        D = mean(dapiStack, 3);
        R = mean(nfkbStack, 3);
        G = mean(phalStack, 3);

    case "trimmedmean"
        ZsD = sort(dapiStack,3,'ascend');
        ZsR = sort(nfkbStack,3,'ascend');
        ZsG = sort(phalStack,3,'ascend');
        zN = size(ZsD,3);
        hi = max(1, round((1-prm.DapiTrimTopFrac)*zN)); % keep bottom (1-frac)
        D  = mean(ZsD(:,:,1:hi), 3);
        R  = mean(ZsR(:,:,1:hi), 3);
        G  = mean(ZsG(:,:,1:hi), 3);

        % [rangeFocus,D] = findBestFocusFromStack(dapiStack);
        % R = mean(nfkbStack(:,:,rangeFocus),3);
        % G = mean(phalStack(:,:,rangeFocus),3);
        
    otherwise
        error("Unknown DapiProjMode. Use 'topK','mip','mean','trimmedMean'.");
end

% G = mean(phalStack, 3);

% Mild denoise (analysis uses these)
Dg = imgaussfilt(D, 1);
Rg = imgaussfilt(R, 1);
Gg = imgaussfilt(G, 2);

% -----------------------
% 1) Nuclei segmentation (NO watershed)
% -----------------------
bw = imbinarize(Dg, "adaptive", "Sensitivity", 0.60);
bw = imfill(bw, "holes");
bw = bwareaopen(bw, 50);
bw = imclearborder(bw);

% --- LoG blob detection for missing nuclei seeds
sigma = 2.0; % 10x?? ? ??? ?? ??(2~4 px ??)
h = fspecial('log', ceil(sigma*6), sigma);

LoG = imfilter(Dg, h, 'replicate');
LoG = -LoG;  % nuclei centers become positive peaks

% find regional maxima as candidate seeds
seedCand = imregionalmax(LoG);

% keep only reasonably strong seeds
thr = prctile(LoG(seedCand), 80);
seedCand = seedCand & (LoG > thr);

% remove seeds that already fall inside existing nuclei mask
seedCand(bw) = 0;

% turn seeds into small disks and add to nuclei mask, then refit
seedDisk = imdilate(seedCand, strel('disk', 3));
bw2 = bw | seedDisk;

% clean
bw2 = imfill(bw2,'holes');
bw2 = bwareaopen(bw2, 50);

% then label and proceed
nucLabel0 = bwlabel(bw2);

% IMPORTANT: no watershed split
% nucLabel0 = bwlabel(bw);

props = regionprops(nucLabel0, Dg, ...
    "Area","Solidity","Eccentricity","MeanIntensity","PixelIdxList");

areas = [props.Area];
sol   = [props.Solidity];
ecc   = [props.Eccentricity];
meanI = [props.MeanIntensity];

keep = areas >= prm.NucMinAreaPx & areas <= prm.NucMaxAreaPx & ...
       sol  >= prm.MinSolidity & ecc <= prm.MaxEccentricity;

% remove "bright tiny outliers" (BBS fragments)
if any(keep)
    small = areas < prctile(areas, 30);
    if any(small)
        cutBright = prctile(meanI(small), 99);
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
% 2) CellMask (optional, conservative)
% -----------------------
if prm.UseCellMask
    th = graythresh(Gg);
    cellMask = (Gg > th*0.6);
    cellMask = imdilate(cellMask, strel("disk", prm.CellMaskDilatePx));
    cellMask = imfill(cellMask, "holes");
    if nnz(cellMask) < 0.02*numel(cellMask)
        cellMask = true(size(cellMask)); % fallback
    end
else
    cellMask = true(size(D));
end

% -----------------------
% 3) Cytoplasm ring per cell (annulus) + overlap exclusion
% -----------------------
cytoLabel = zeros(size(nucLabel), "like", nucLabel);
bgRingLabel = zeros(size(nucLabel), "like", nucLabel);

allNuc = nucLabel > 0;

seIn  = strel("disk", prm.CytoInnerRadiusPx);
seOut = strel("disk", prm.CytoOuterRadiusPx);
seBg  = strel("disk", prm.BgOuterRadiusPx);

% To exclude overlap: build a "dilated nuclei map" for all nuclei
allDilOut = imdilate(allNuc, seOut);

for i = 1:nCells
    thisNuc = (nucLabel == i);

    % Cytoplasm annulus: outer dilated - inner dilated
    inner = imdilate(thisNuc, seIn);
    outer = imdilate(thisNuc, seOut);
    ring  = outer & ~inner;

    % Limit to cellMask and exclude all nuclei
    ring = ring & cellMask;
    ring(allNuc) = 0;

    % Exclude areas closer to other nuclei (overlap exclusion):
    % If ring pixels also belong to other nuclei dilations, remove them.
    otherDil = allDilOut & ~imdilate(thisNuc, seOut);
    ring(otherDil) = 0;

    cytoLabel(ring) = i;

    % Background ring: outer-bg region just outside cytoplasm annulus
    bgOuter = imdilate(thisNuc, seBg);
    bgRing = bgOuter & ~outer;        % outside cyto outer
    bgRing = bgRing & cellMask;
    bgRing(allNuc) = 0;
    % Exclude overlap zones as well
    bgRing(otherDil) = 0;
    bgRingLabel(bgRing) = i;
end

% -----------------------
% 4) Intensities + local background subtraction
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

    if nucArea(i) < 50 || cytoArea(i) < 50 || bgArea(i) < 50
        continue
    end

    bg = median(Rg(bg_i), "omitnan");
    In = mean(Rg(nuc_i), "omitnan")  - bg;
    Ic = mean(Rg(cyto_i), "omitnan") - bg;

    Inuc(i)  = In;
    Icyto(i) = Ic;
    Ibg(i)   = bg;

    NC(i) = In / (Ic + eps);
end

valid = isfinite(NC) & Inuc>0 & Icyto>0;

T = table((1:nCells)', nucArea, cytoArea, bgArea, Inuc, Icyto, Ibg, NC, valid, ...
    'VariableNames', {'cellID','nucAreaPx','cytoAreaPx','bgAreaPx','Inuc','Icyto','Ibg','NCratio','valid'});

if prm.Verbose
    fprintf('MD: %s\n', md.movieDataPath_);
    fprintf('Cells: %d | valid: %d (%.1f%%)\n', nCells, nnz(valid), 100*nnz(valid)/nCells);
    if nnz(valid)>0
        fprintf('NCratio median(valid)=%.3f\n', median(T.NCratio(valid), "omitnan"));
    end
end

out = struct();
out.table = T;
out.images = struct('DAPIproj', D, 'NFKBproj', R, 'PHALproj', G);
out.masks  = struct('nucLabel', nucLabel, 'cytoLabel', cytoLabel, ...
                    'bgRingLabel', bgRingLabel, 'cellMask', cellMask);
out.params = prm;
end