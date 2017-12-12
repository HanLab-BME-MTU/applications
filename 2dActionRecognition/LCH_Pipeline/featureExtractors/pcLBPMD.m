%% pcLBPMD - calculates LBP for single cell trajectories
% calculate multi-scale LBP for entire frames, then crops specific cell and extract
% corresponding LBP values for every cell's trajectory

function [] = pcLBPMD(MD,params,dirs)

lbpMapping = getmapping(8,'riu2');

load([dirs.tracking 'cellIdTYX.mat']);% cellTYX

nCells = length(cellTYX); %#ok<USENS>

accumulatedFovPyramidLBP = cell(1,params.nScales); % take LBP for all field of view
accumulatedBckPyramidLBP = cell(1,params.nScales); % take LBP for all bck in FOV
accumulatedFwdPyramidLBP = cell(1,params.nScales); % take LBP for all fwd in FOV
for iScale = 1 : params.nScales
    accumulatedFovPyramidLBP{iScale} = [];
end

lbpData.fov = cell(1,nCells);
lbpData.bck = cell(1,nCells);
lbpData.fwd = cell(1,nCells);

lbpData.t0 = nan(1,nCells); % NEW: initial time for the trajectory

lbpFname = [dirs.tracking filesep 'lbpData.mat'];

if exist(lbpFname,'file') && ~params.always
    fprintf(sprintf('LBP for %s exists, finishing\n',lbpFname));
    return;
end

if ~isfield(params,'sTime')
    params.sTime = 1;
end

%% LBP for frames - do not execute even if "always"!
for t = params.sTime : params.nTime - params.frameJump
    lbpTFname = [dirs.lbp sprintf('%03d',t) '_lbp.mat'];
    
    % Just to correct a bug
    assert(exist([dirs.lbp sprintf('%03d',t) '_mf.mat'],'file'));
    if exist([dirs.lbp sprintf('%03d',t) '_mf.mat'],'file')
        unix(sprintf('mv %s %s',[dirs.lbp sprintf('%03d',t) '_mf.mat'],lbpTFname));
    end
    
    if exist(lbpTFname,'file')
        if params.deepDebug
            fprintf(sprintf('LBP frame %d: continue\n',t));
        end
        continue;
    end
    
    fprintf(sprintf('LBP frame %d\n',t));
    
    I = MD.getChannel(1).loadImage(t);
    
    pyramidLBP = cell(1,params.nScales);
    for iScale = 1 : params.nScales
        curScale = params.scales(iScale);
        curI = imresize(I,curScale);
        
        IlbpTmp = lbp(curI,1,8,lbpMapping,'');
        Ilbp = nan(size(curI));
        Ilbp(2:end-1,2:end-1) = IlbpTmp;
        pyramidLBP{iScale} = Ilbp;        
    end
    
    save(lbpTFname,'pyramidLBP');
end

I = MD.getChannel(1).loadImage(t);

%% LBP for cells
for icell = 1 : nCells    
    fprintf(sprintf('LBP cell %d/%d\n',icell,nCells));
    curCell = cellTYX{icell};
    ntime = length(curCell.ts);
    t0 = curCell.ts(1);
    
    lbpData.t0(icell) = t0; % NEW: initial time for the trajectory
    
    lbpData.fov{icell}.pyramidLBP = cell(1,params.nScales);
    lbpData.bck{icell}.pyramidLBP = cell(1,params.nScales);
    
%     % ROI
%     pyramidROI = cell(1,params.nScales);
%     for iScale = 1 : params.nScales
%         curScale = params.scales(iScale);
%         curI = imresize(I,curScale);
%         
%         IlbpTmp = lbp(curI,1,8,lbpMapping,'');
%         Ilbp = nan(size(curI));
%         Ilbp(2:end-1,2:end-1) = IlbpTmp;
%         pyramidLBP{iScale} = Ilbp;
%     end
    load([dirs.roiData sprintf('%d',icell) '_roi.mat']); % curCellRoi for the bck calculations  

    for iScale = 1 : params.nScales
        curScale = params.scales(iScale);
        lbpData.fov{icell}.pyramidLBP{iScale}.lbp = nan(ntime,10);  
        lbpData.bck{icell}.pyramidLBP{iScale}.lbp = nan(ntime,10);  
        
        curFOVRadius = round(params.FOVRadius*curScale);
        
        for t = curCell.ts
            curT = t - t0 + 1;
            curX = round(curCell.xs(curT)*curScale);            
            curY = round(curCell.ys(curT)*curScale);
            
            bby0 = max(1,curY - curFOVRadius);
            bby1 = min(size(I,1),curY + curFOVRadius);
            bbx0 = max(1,curX - curFOVRadius);
            bbx1 = min(size(I,2),curX + curFOVRadius);
                
            lbpTFname = [dirs.lbp sprintf('%03d',t) '_lbp.mat'];
            load(lbpTFname); % pyramidLBP
            Ilbp = pyramidLBP{iScale};
            IlbpBB = Ilbp(bby0:bby1,bbx0:bbx1);
            
            %         % if segmentation
            %         lbpBin = lbpBB(cells.seg{icells});
        
            lbpDesc = hist(IlbpBB(:),0:9);
            lbpDesc = lbpDesc ./ sum(lbpDesc);
            accumulatedFovPyramidLBP{iScale} = [accumulatedFovPyramidLBP{iScale},lbpDesc'];       
            lbpData.fov{icell}.pyramidLBP{iScale}.lbp(curT,:) = lbpDesc;   
            
            %% TODO ROI
            roibby0 = curCellRoi{curT}.bby0; %#ok<USENS>
            roibby1 = curCellRoi{curT}.bby1;
            roibbx0 = curCellRoi{curT}.bbx0;
            roibbx1 = curCellRoi{curT}.bbx1;
            
            ROIBB = curCellRoi{curT}.roi;
            ROI = false(size(I));
            ROI(roibby0:roibby1,roibbx0:roibbx1) = ROIBB;
            ROIScale = imresize(ROI,curScale,'nearest');
            
            Ibb = Ilbp(bby0:bby1,bbx0:bbx1);
            ROIBBbb = ROIScale(bby0:bby1,bbx0:bbx1);
            IlbpBck =  Ibb(~ROIBBbb);
            IlbpFwd =  Ibb(ROIBBbb);
            
            lbpBckDesc = hist(IlbpBck(:),0:9);
            lbpBckDesc = lbpBckDesc ./ sum(lbpBckDesc);
            
            lbpFwdDesc = hist(IlbpFwd(:),0:9);
            lbpFwdDesc = lbpFwdDesc ./ sum(lbpFwdDesc);
            
            accumulatedBckPyramidLBP{iScale} = [accumulatedBckPyramidLBP{iScale},lbpBckDesc'];               
            lbpData.bck{icell}.pyramidLBP{iScale}.lbp(curT,:) = lbpBckDesc;               
            
            accumulatedFwdPyramidLBP{iScale} = [accumulatedFwdPyramidLBP{iScale},lbpFwdDesc'];               
            lbpData.fwd{icell}.pyramidLBP{iScale}.lbp(curT,:) = lbpFwdDesc;               
        end
    end    
end
save(lbpFname,'lbpData','accumulatedFovPyramidLBP','accumulatedBckPyramidLBP','accumulatedFwdPyramidLBP');
end