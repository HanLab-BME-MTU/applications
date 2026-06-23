function [bodyMask, level] = segmentFilopodiaBody(img, p, manualLevel)
%SEGMENTFILOPODIABODY  Cell-body mask for the FilopodiaForcePackage P1 step.
%Shared by segmentMovieFilopodia (batch) and the segmentation settings GUI
%preview, so the preview always matches what the process will produce.
%
%   [bodyMask, level] = segmentFilopodiaBody(img, p)
%   [bodyMask, level] = segmentFilopodiaBody(img, p, manualLevel)
%
% img         : grayscale frame (double)
% p           : funParams struct (GaussianBlurSigma, BodyThreshold,
%               ThreshScale, BodyOpenRadius, BodyClosingRadius, BodyMinArea)
% manualLevel : optional fixed intensity threshold; when given (non-empty),
%               it overrides the method/scale and is applied directly.
% level       : the intensity threshold actually used (after ThreshScale).
% Sangyoon J. Han / 2026

if nargin < 3, manualLevel = []; end

imgB = imgaussfilt(img, gf(p,'GaussianBlurSigma',2));

if ~isempty(manualLevel)
    level = manualLevel;                 % manual override (GUI slider)
else
    switch lower(num2str(gf(p,'BodyThreshold','rosin')))
        case 'otsu',  level = thresholdOtsu(imgB);
        case 'rosin', level = thresholdRosin(imgB);
        otherwise,    level = gf(p,'BodyThreshold',[]);
    end
    ts = gf(p,'ThreshScale',1);
    if ~isempty(ts) && isnumeric(ts) && ts > 0
        level = level * ts;              % <1 recovers dim talin band
    end
end

bodyMask = imgB > level;
bodyMask = imfill(bodyMask, 'holes');
if any(bodyMask(:))
    bodyMask = bwareafilt(bodyMask, 1);              % keep largest object
end
if gf(p,'BodyOpenRadius',0) > 0
    bodyMask = imopen(bodyMask, strel('disk', round(p.BodyOpenRadius)));
end
if gf(p,'BodyClosingRadius',0) > 0
    bodyMask = imclose(bodyMask, strel('disk', round(p.BodyClosingRadius)));
end
bodyMask = imfill(bodyMask, 'holes');
bodyMask = bwareaopen(bodyMask, gf(p,'BodyMinArea',500));
end

% =====================================================================
function v = gf(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
