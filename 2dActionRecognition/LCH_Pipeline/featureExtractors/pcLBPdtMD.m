%% pcLBPdtMD - calculates LBPdt = LBP(I1-I0) for single cell trajectories
% calculate multi-scale LBPdt for entire frames, then crops specific cell and extract
% corresponding LBPdt values for every cell's trajectory

function [] = pcLBPdtMD(MD,params,dirs)

lbpMapping = getmapping(8,'riu2');

load([dirs.tracking 'cellIdTYX.mat']);% cellTYX

nCells = length(cellTYX); %#ok<USENS>

accumulatedFovPyramidLBPdt = cell(1,params.nScales); % take LBP for all field of view
for iScale = 1 : params.nScales
    accumulatedFovPyramidLBPdt{iScale} = [];
end

lbpDtData.fov = cell(1,nCells);
lbpDtData.t0 = nan(1,nCells);

lbpDtFname = [dirs.tracking filesep 'lbpDtData.mat'];

if exist(lbpDtFname,'file') && ~params.always
    fprintf(sprintf('LBPdt for %s exists, finishing\n',lbpDtFname));
    return;
end

if ~isfield(params,'sTime')
    params.sTime = 1;
end

for t = params.sTime : params.nTime - params.frameJump
    lbpDtTFname = [dirs.lbpDt sprintf('%03d',t) '_lbpDt.mat'];
    
    if exist(lbpDtTFname,'file')
        continue;
    end
    
    fprintf(sprintf('LBP Dt frame %d\n',t));
    
    I0 = MD.getChannel(1).loadImage(t);
    I1 = MD.getChannel(1).loadImage(t+1);
    
    I0 = (I0 - mean(I0(:)))./std(double(I0(:)));
    I1 = (I1 - mean(I1(:)))./std(double(I1(:)));
    
    Idt = I0 - I1;
    
    pyramidLbpDt = cell(1,params.nScales);
    for iScale = 1 : params.nScales
        curScale = params.scales(iScale);
        curIdt = imresize(Idt,curScale);
        
        ILbpDtTmp = lbp(curIdt,1,8,lbpMapping,'');
        ILbpDt = nan(size(curIdt));
        ILbpDt(2:end-1,2:end-1) = ILbpDtTmp;
        pyramidLbpDt{iScale} = ILbpDt;
    end
    
    save(lbpDtTFname,'pyramidLbpDt');
end

% I = MD.getChannel(1).loadImage(t);

%% dtLBP for cells
for icell = 1 : nCells
    fprintf(sprintf('dtLBP cell %d/%d\n',icell,nCells));
    curCell = cellTYX{icell};
    ntime = length(curCell.ts);
    t0 = curCell.ts(1);
    
    lbpDtData.t0(icell) = t0; % NEW: initial time for the trajectory
    
    lbpDtData.fov{icell}.pyramidLbpDt = cell(1,params.nScales);
    %     lbpDtData.bck{icell}.pyramidLbpDt = cell(1,params.nScales);
    
    %     % ROI
    %     pyramidROI = cell(1,params.nScales);
    %     for iScale = 1 : params.nScales
    %         curScale = params.scales(iScale);
    %         curI = imresize(I,curScale);
    %
    %         IlbpTmp = lbp(curI,1,8,lbpMapping,'');
    %         Ilbp = nan(size(curI));
    %         Ilbp(2:end-1,2:end-1) = IlbpTmp;
    %         pyramidLbpDt{iScale} = Ilbp;
    %     end
    %     load([dirs.roiData sprintf('%d',icell) '_roi.mat']); % curCellRoi for the bck calculations
    
    for iScale = 1 : params.nScales
        curScale = params.scales(iScale);
        lbpDtData.fov{icell}.pyramidLbpDt{iScale}.lbpDt = nan(ntime,10);
        %         lbpDtData.bck{icell}.pyramidLbpDt{iScale}.lbpDt = nan(ntime,10);
        
        curFOVRadius = round(params.FOVRadius*curScale);
        
        for t = curCell.ts
            curT = t - t0 + 1;
            curX = round(curCell.xs(curT)*curScale);
            curY = round(curCell.ys(curT)*curScale);            
            
            lbpDtTFname = [dirs.lbpDt sprintf('%03d',t) '_lbpDt.mat'];
            load(lbpDtTFname); % pyramidLbpDt
            IlbpDt = pyramidLbpDt{iScale};            
            
            bby0 = max(1,curY - curFOVRadius);
            bby1 = min(size(IlbpDt,1),curY + curFOVRadius);
            bbx0 = max(1,curX - curFOVRadius);
            bbx1 = min(size(IlbpDt,2),curX + curFOVRadius);
            
            IlbpDtBB = IlbpDt(bby0:bby1,bbx0:bbx1);
            
            %         % if segmentation
            %         lbpBin = lbpBB(cells.seg{icells});
            
            lbpDtDesc = hist(IlbpDtBB(:),0:9);
            lbpDtDesc = lbpDtDesc ./ sum(lbpDtDesc);
            accumulatedFovPyramidLBPdt{iScale} = [accumulatedFovPyramidLBPdt{iScale},lbpDtDesc'];
            lbpDtData.fov{icell}.pyramidLbpDt{iScale}.lbpDt(curT,:) = lbpDtDesc;
            
            %             %% TODO ROI
            %             roibby0 = curCellRoi{curT}.bby0; %#ok<USENS>
            %             roibby1 = curCellRoi{curT}.bby1;
            %             roibbx0 = curCellRoi{curT}.bbx0;
            %             roibbx1 = curCellRoi{curT}.bbx1;
            %
            %             ROIBB = curCellRoi{curT}.roi;
            %             ROI = false(size(I));
            %             ROI(roibby0:roibby1,roibbx0:roibbx1) = ROIBB;
            %             ROIScale = imresize(ROI,curScale,'nearest');
            %
            %             Ibb = IlbpDt(bby0:bby1,bbx0:bbx1);
            %             ROIBBbb = ROIScale(bby0:bby1,bbx0:bbx1);
            %             IlbpBck =  Ibb(~ROIBBbb);
            %             IlbpFwd =  Ibb(ROIBBbb);
            %
            %             lbpBckDesc = hist(IlbpBck(:),0:9);
            %             lbpBckDesc = lbpBckDesc ./ sum(lbpBckDesc);
            %
            %             lbpFwdDesc = hist(IlbpFwd(:),0:9);
            %             lbpFwdDesc = lbpFwdDesc ./ sum(lbpFwdDesc);
            %
            %             accumulatedBckPyramidLBP{iScale} = [accumulatedBckPyramidLBP{iScale},lbpBckDesc'];
            %             lbpDtData.bck{icell}.pyramidLbpDt{iScale}.lbpDt(curT,:) = lbpBckDesc;
            %
            %             accumulatedFwdPyramidLBP{iScale} = [accumulatedFwdPyramidLBP{iScale},lbpFwdDesc'];
            %             lbpDtData.fwd{icell}.pyramidLbpDt{iScale}.lbpDt(curT,:) = lbpFwdDesc;
        end
    end
end
% save(lbpDtFname,'lbpDtData','accumulatedFovPyramidLBPdt','accumulatedBckPyramidLBP','accumulatedFwdPyramidLBP');
save(lbpDtFname,'lbpDtData','accumulatedFovPyramidLBPdt');
end