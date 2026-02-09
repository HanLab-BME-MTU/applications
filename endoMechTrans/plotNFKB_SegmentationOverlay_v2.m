function hFig = plotNFKB_SegmentationOverlay_v2(out, varargin)
% plotNFKB_SegmentationOverlay_v2
% QC overlay with robust contrast stretching for display.
%
% INPUT:
%   out: struct from quantifyNFKB_NCratio_MD_v2
%
% Params:
%   'BaseImage'   : 'nfkb'(default) | 'dapi' | 'phal'
%   'OnlyValid'   : true/false (default false)
%   'ShowCellMask': true/false (default true)
%   'ShowLabels'  : true/false (default false)
%   'ClipPrct'    : [low high] percentiles for display (default [1 99])
%   'Gamma'       : gamma for imadjust (default 0.7)
%   'SavePath'    : '' or file path
%   'Title'       : custom title

p = inputParser;
p.addParameter('BaseImage','nfkb',@(s)ischar(s)||isstring(s));
p.addParameter('OnlyValid',false,@(x)islogical(x)&&isscalar(x));
p.addParameter('ShowCellMask',true,@(x)islogical(x)&&isscalar(x));
p.addParameter('ShowLabels',false,@(x)islogical(x)&&isscalar(x));
p.addParameter('ClipPrct',[1 99],@(x)isnumeric(x)&&numel(x)==2);
p.addParameter('Gamma',0.7,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('SavePath','',@(s)ischar(s)||isstring(s));
p.addParameter('Title','',@(s)ischar(s)||isstring(s));
p.parse(varargin{:});
prm = p.Results;

% Base image selection
baseChoice = lower(string(prm.BaseImage));
switch baseChoice
    case "nfkb"
        I = out.images.NFKBproj;
        baseName = 'NFkB (proj)';
    case "dapi"
        I = out.images.DAPIproj;
        baseName = 'DAPI (proj)';
    case "phal"
        I = out.images.PHALproj;
        baseName = 'Phalloidin (proj)';
    otherwise
        error("BaseImage must be 'nfkb', 'dapi', or 'phal'.");
end
I = im2double(I);

% Robust display scaling (IMPORTANT: display only)
lo = prctile(I(:), prm.ClipPrct(1));
hi = prctile(I(:), prm.ClipPrct(2));
if hi <= lo
    Ishow = mat2gray(I);
else
    Ishow = imadjust(I, [lo hi], [0 1], prm.Gamma);
end

% Masks
T = out.table;
nucLabel = out.masks.nucLabel;
cytoLabel = out.masks.cytoLabel;
bgLabel = out.masks.bgRingLabel;
cellMask = out.masks.cellMask;

if prm.OnlyValid
    validIDs = T.cellID(T.valid);
    nucBW  = ismember(nucLabel, validIDs);
    cytoBW = ismember(cytoLabel, validIDs);
    bgBW   = ismember(bgLabel, validIDs);
else
    nucBW  = nucLabel > 0;
    cytoBW = cytoLabel > 0;
    bgBW   = bgLabel > 0;
end

bNuc  = bwperim(nucBW);
bCyto = bwperim(cytoBW);
bBg   = bwperim(bgBW);
bCell = bwperim(cellMask);

% Compose RGB overlay
RGB = repmat(Ishow,1,1,3);

% Nucleus boundary: white
RGB(:,:,1) = max(RGB(:,:,1), bNuc);
RGB(:,:,2) = max(RGB(:,:,2), bNuc);
RGB(:,:,3) = max(RGB(:,:,3), bNuc);

% Cytoplasm ring boundary: yellow-ish (R+G)
RGB(:,:,1) = max(RGB(:,:,1), bCyto);
RGB(:,:,2) = max(RGB(:,:,2), bCyto);

% Background ring boundary: cyan (G+B) - optional but helpful QC
RGB(:,:,2) = max(RGB(:,:,2), bBg);
RGB(:,:,3) = max(RGB(:,:,3), bBg);

% CellMask boundary: blue
if prm.ShowCellMask
    RGB(:,:,3) = max(RGB(:,:,3), bCell);
end

hFig = figure('Color','w');
imshow(RGB, []);
hold on;

% Optional labels
if prm.ShowLabels
    stats = regionprops(nucLabel, 'Centroid');
    for i = 1:numel(stats)
        if ~isempty(stats(i).Centroid)
            text(stats(i).Centroid(1), stats(i).Centroid(2), sprintf('%d', i), ...
                'Color','y','FontSize',8,'FontWeight','bold');
        end
    end
end

if strlength(prm.Title) > 0
    title(prm.Title, 'Interpreter','none');
else
    title(sprintf('%s + segmentation overlay', baseName), 'Interpreter','none');
end

% Small legend
txt = {
    'White: nuclei boundary'
    'Yellow: cytoplasm ring boundary'
    'Cyan: background ring boundary'
};
if prm.ShowCellMask
    txt{end+1} = 'Blue: cellMask boundary';
end
text(10, 15, txt, 'Color','w','FontSize',9, 'Interpreter','none', ...
    'BackgroundColor',[0 0 0 0.35], 'Margin',4);

hold off;

% Save if requested
if strlength(prm.SavePath) > 0
    [saveDir,~,~] = fileparts(prm.SavePath);
    if ~isempty(saveDir) && ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
    exportgraphics(hFig, prm.SavePath, 'Resolution', 200);
end
end