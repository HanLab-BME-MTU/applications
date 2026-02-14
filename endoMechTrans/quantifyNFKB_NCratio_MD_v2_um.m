function out = quantifyNFKB_NCratio_MD_v2_um(md, ringUm, varargin)
% Calls quantifyNFKB_NCratio_MD_v2 with radii specified in microns (um).
% md.pixelSize_ is assumed nm/pixel.

pxSize_nm = md.pixelSize_;
umPerPx = pxSize_nm / 1000;  % nm->um

toPx = @(um) max(1, round(um / umPerPx));

cyIn  = toPx(ringUm.CytoInner);
cyOut = toPx(ringUm.CytoOuter);
bgOut = toPx(ringUm.BgOuter);

out = quantifyNFKB_NCratio_MD_v2(md, ...
    'CytoInnerRadiusPx', cyIn, ...
    'CytoOuterRadiusPx', cyOut, ...
    'BgOuterRadiusPx',   bgOut, ...
    varargin{:});

out.params_um = ringUm;
out.params_px = struct('CytoInnerRadiusPx',cyIn,'CytoOuterRadiusPx',cyOut,'BgOuterRadiusPx',bgOut,'umPerPx',umPerPx);
end