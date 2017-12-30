%% pcCellRoiLeverMD - takes ROI masks from LEVER segmentation
% Nov. 2017
% mask taken from dirs.lever and saved to the local file

function [] = pcCellRoiLeverMD(MD,params,dirs)

load([dirs.tracking 'cellIdTYX.mat']);% cellTYX

nCells = length(cellTYX); %#ok<USENS>

structElement = strel('square',round(3.0/params.pixelSize));

%% ROI for cells
for icell = 1 : nCells
    roiLeverFname = [dirs.roiLever sprintf('%d',icell) '_roi.mat'];
    
    if exist(roiLeverFname,'file') && ~params.always
        if params.deepDebug
            fprintf(sprintf('cell %d roi LEVER: continue\n',icell));
        end
        
        continue;
    end
    
    fprintf(sprintf('ROI LEVER cell %d/%d\n',icell,nCells));
    
    curCell = cellTYX{icell};
    ntime = length(curCell.ts);
    
    curCellRoi = cell(1,ntime);
    
    t0 = curCell.ts(1);
    for t = curCell.ts
        curT = t - t0 + 1;
        curX = round(curCell.xs(curT));
        curY = round(curCell.ys(curT));
        
        I = MD.getChannel(1).loadImage(curT);
        
        inLeverFname = [dirs.lever filesep reduceZeros4Lever(dirs.expname) '.LEVER' filesep num2str(t) '.mat'];
        
        if ~exist(inLeverFname,'file')
            continue;
        end
        
        load(inLeverFname); % ROI_LEVER
        
        bby0 = max(1,curY - params.FOVRadius);
        bby1 = min(size(ROI_LEVER,1),curY + params.FOVRadius); % FOVRadius = 35um (RoiRadius = 120 um)
        bbx0 = max(1,curX - params.FOVRadius);
        bbx1 = min(size(ROI_LEVER,2),curX + params.FOVRadius);
        
        curI = I(bby0:bby1,bbx0:bbx1);
        MASK = ROI_LEVER(bby0:bby1,bbx0:bbx1);
        
        yBB = curY - bby0;
        xBB = curX - bbx0;
        
        [labels,nLabels] = bwlabel(MASK,8);
        
        ROI = [];
        for iLabels = 1 : nLabels
            curLabel = (MASK == iLabels);
            curLabel1 = imdilate(curLabel,structElement);
            curLabel2 = imfill(curLabel1,'holes');
            if curLabel2(yBB,xBB)
                ROI = curLabel2;
                break;
            end
        end
        
        if (isempty(ROI) && (curT == 1))                   
            curCellRoi = [];            
            break;
        end
        
        if isempty(ROI)
            ROI = prevROI;
        end
        
        prevROI = ROI;
        
        curCellRoi{curT}.bby0 = bby0;
        curCellRoi{curT}.bby1 = bby1;
        curCellRoi{curT}.bbx0 = bbx0;
        curCellRoi{curT}.bbx1 = bbx1;
        
        curCellRoi{curT}.roi = ROI;
        curCellRoi{curT}.img = curI;
        
        curCellRoi{curT}.x = curX;
        curCellRoi{curT}.y = curY;
        curCellRoi{curT}.xBB = xBB;
        curCellRoi{curT}.yBB = yBB;
        curCellRoi{curT}.t = t;                
    end
    save(roiLeverFname,'curCellRoi');    
end
end

function reducedStr = reduceZeros4Lever(origStr)
if strcmp(origStr(end-1),'0')
    reducedStr = [origStr(1:end-2) origStr(end)];
else
    reducedStr = origStr;
end
end


