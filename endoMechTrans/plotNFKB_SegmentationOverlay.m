function hFig = plotNFKB_SegmentationOverlay(out, varargin)
% plotNFKB_SegmentationOverlay
% QC overlay for nuclei + territory/cytoplasm segmentation.
%
% INPUT
%   out: struct returned by quantifyNFKB_NCratio_MD
%
% NAME-VALUE PARAMS
%   'BaseImage'   : 'nfkb' (default) or 'dapi' or 'phal'
%   'ShowCellMask': true/false (default true)
%   'ShowLabels'  : true/false (default false) % show cellID at nucleus centroid
%   'OnlyValid'   : true/false (default false) % show only valid cells
%   'SavePath'    : '' (default) or full file path to save png
%   'Title'       : '' custom title
%
% OUTPUT
%   hFig: figure handle

p = inputParser;
p.addParameter('BaseImage','nfkb',@(s)ischar(s)||isstring(s));
p.addParameter('ShowCellMask',true,@(x)islogical(x)&&isscalar(x));
p.addParameter('ShowLabels',false,@(x)islogical(x)&&isscalar(x));
p.addParameter('OnlyValid',false,@(x)islogical(x)&&isscalar(x));
p.addParameter('SavePath','',@(s)ischar(s)||isstring(s));
p.addParameter('Title','',@(s)ischar(s)||isstring(s));
p.parse(varargin{:});
prm = p.Results;

% pick base image
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
% Ishow = mat2gray(I);  % normalize for display
% Robust contrast for display
pLow  = prctile(I(:), 1);
pHigh = prctile(I(:), 99);
Ishow = imadjust(I, [pLow pHigh], [0 1]);

% ?? ? ???
% Ishow = imadjust(I, [pLow pHigh], [0 1], 0.5); % gamma

nucLabel = out.masks.nucLabel;
terrLabel = out.masks.territoryLabel;
cytoLabel = out.masks.cytoLabel;
cellMask = out.masks.cellMask;

T = out.table;

% Optionally restrict to valid cells
if prm.OnlyValid
    validIDs = T.cellID(T.valid);
    keepTerr = ismember(terrLabel, validIDs);
    keepNuc  = ismember(nucLabel, validIDs);
    keepCyto = ismember(cytoLabel, validIDs);

    terrLabel2 = terrLabel; terrLabel2(~keepTerr) = 0;
    nucLabel2  = nucLabel;  nucLabel2(~keepNuc) = 0;
    cytoLabel2 = cytoLabel; cytoLabel2(~keepCyto) = 0;
else
    terrLabel2 = terrLabel;
    nucLabel2  = nucLabel;
    cytoLabel2 = cytoLabel;
end

% boundaries
nucBW  = nucLabel2 > 0;
terrBW = terrLabel2 > 0;
cytoBW = cytoLabel2 > 0;

bNuc  = bwperim(nucBW);
bTerr = bwperim(terrBW);
bCyto = bwperim(cytoBW);
bCell = bwperim(cellMask);

% Create RGB overlay
RGB = repmat(Ishow,1,1,3);

% Overlay colors (do not specify ?pretty?; keep simple)
% Nuclei boundary -> white
RGB(:,:,1) = max(RGB(:,:,1), bNuc);
RGB(:,:,2) = max(RGB(:,:,2), bNuc);
RGB(:,:,3) = max(RGB(:,:,3), bNuc);

% Territory boundary -> green
RGB(:,:,2) = max(RGB(:,:,2), bTerr);

% Cytoplasm boundary -> red (optional; can be busy)
% Here we use a thinner cue: only cyto boundary where not nucleus boundary
bCytoOnly = bCyto & ~bNuc;
RGB(:,:,1) = max(RGB(:,:,1), bCytoOnly);

% CellMask boundary -> blue (optional)
if prm.ShowCellMask
    RGB(:,:,3) = max(RGB(:,:,3), bCell);
end

% Plot
hFig = figure('Color','w');
imshow(RGB, []);
hold on;

% Optionally show labels at nucleus centroid
if prm.ShowLabels
    stats = regionprops(nucLabel2, 'Centroid');
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

% Legend text (simple)
txt = {
    'White: nucleus boundary'
    'Green: territory boundary'
    'Red: cytoplasm boundary'
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