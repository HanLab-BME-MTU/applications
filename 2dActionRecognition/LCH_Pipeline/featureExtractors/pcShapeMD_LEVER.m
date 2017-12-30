%% pcShapeMD_LEVER - calculates shape features based on LEVER segmentation masks
% Assaf Zaritsky, December 2017

% MD - not used (image just used for size), for future use...
function [] = pcShapeMD_LEVER(MD,params,dirs)

load([dirs.tracking 'cellIdTYX.mat']);% cellTYX

nCells = length(cellTYX); %#ok<USENS>

shapeData.shapeFeats = cell(1,nCells);
shapeData.t0 = nan(1,nCells); % NEW: initial time for the trajectory

shapeFname = [dirs.tracking filesep 'shapeDataLever.mat'];

if exist(shapeFname,'file') && ~params.always
    fprintf(sprintf('Shape features using LEVER ROI for %s exists, finishing\n',shapeFname));
    return;
end

if ~isfield(params,'sTime')
    params.sTime = 1;
end

%% Extract shape features for each cell
for icell = 1 : nCells
    fprintf(sprintf('Extract shape from LEVERs ROI cell %d/%d\n',icell,nCells));
    curCell = cellTYX{icell};
    ntime = length(curCell.ts);
    t0 = curCell.ts(1);
    
    shapeData.t0(icell) = t0; % initial time for the trajectory    
    shapeData.shapeFeats{icell}.feats = nan(ntime,7);
    
    load([dirs.roiLever sprintf('%d',icell) '_roi.mat']); % curCellRoi.roi
    
    if isempty(curCellRoi)
        shapeData.shapeFeats{icell}.feats = [];
        continue;
    end
    
    
    for t = curCell.ts
        curT = t - t0 + 1;
        
        I = MD.getChannel(1).loadImage(curT);
        
        %% ROI        
        ROIBB = curCellRoi{curT}.roi;
                     
        %% Shape 
        shapeData.shapeFeats{icell}.feats(curT,:) = getShapeFeats(ROIBB);
    end
end
save(shapeFname,'lbpData');
end

function shapeFeats = getShapeFeats(ROI)
stats = regionprops(ROI, 'Area', 'Eccentricity', 'Perimeter', 'Solidity', 'MajorAxisLength','MinorAxisLength');
shapeFeats = [stats.Area, stats.Eccentricity, stats.Perimeter, stats.Solidity, stats.MajorAxisLength, stats.MinorAxisLength, stats.MinorAxisLength/stats.MajorAxisLength];
end