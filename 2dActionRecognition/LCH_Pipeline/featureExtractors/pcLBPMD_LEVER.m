%% pcLBPMD_LEVER - calculates LBP for single cell trajectories, within LEVER segmentation masks
% calculate multi-scale LBP for entire frames, then crops specific cell and extract
% corresponding LBP values for every cell's trajectory
% Assaf Zaritsky, December 2017

% MD - not used (image just used for size), for future use...
function [] = pcLBPMD_LEVER(MD,params,dirs)

load([dirs.tracking 'cellIdTYX.mat']);% cellTYX

curScale = 1; % not doing multi-scale analysis
% curFOVRadius = round(params.FOVRadius*curScale);

nCells = length(cellTYX); %#ok<USENS>

% accumulatedFovLBP = []; % take LBP for all field of view
% accumulatedBckLBP = []; % take LBP for all bck in FOV
% accumulatedFwdLBP = []; % take LBP for all fwd in FOV

lbpData.fov = cell(1,nCells);
lbpData.bck = cell(1,nCells);
lbpData.fwd = cell(1,nCells);

lbpData.corr = cell(1,nCells);

lbpData.t0 = nan(1,nCells); % NEW: initial time for the trajectory

lbpFname = [dirs.tracking filesep 'lbpDataLever.mat'];

if exist(lbpFname,'file') && ~params.always
    fprintf(sprintf('LBP features in LEVER ROI for %s exists, finishing\n',lbpFname));
    return;
end

if ~isfield(params,'sTime')
    params.sTime = 1;
end

%% LBP for cells
for icell = 1 : nCells
    fprintf(sprintf('LBP in LEVERs ROI cell %d/%d\n',icell,nCells));
    curCell = cellTYX{icell};
    ntime = length(curCell.ts);
    t0 = curCell.ts(1);
    
    lbpData.t0(icell) = t0; % initial time for the trajectory
    
    lbpData.fov{icell}.lbp = nan(ntime,10);
    lbpData.bck{icell}.lbp = nan(ntime,10);
    lbpData.fwd{icell}.lbp = nan(ntime,10);
    
    load([dirs.roiLever sprintf('%d',icell) '_roi.mat']); % curCellRoi.roi for the bck calculations
    
    if isempty(curCellRoi)
        lbpData.fov{icell} = [];
        lbpData.bck{icell} = [];
        lbpData.fwd{icell} = [];
        lbpData.corr{icell} = [];
        continue;
    end
    
    
    for t = curCell.ts
        curT = t - t0 + 1;
        %         curX = round(curCell.xs(curT)*curScale);
        %         curY = round(curCell.ys(curT)*curScale);
        
        I = MD.getChannel(1).loadImage(curT);
        
        %% ROI
        roibby0 = curCellRoi{curT}.bby0; %#ok<USENS>
        roibby1 = curCellRoi{curT}.bby1;
        roibbx0 = curCellRoi{curT}.bbx0;
        roibbx1 = curCellRoi{curT}.bbx1;
        
        ROIBB = curCellRoi{curT}.roi;
        
        lbpTFname = [dirs.lbp sprintf('%03d',t) '_lbp.mat'];
        load(lbpTFname); % pyramidLBP
        Ilbp = pyramidLBP{curScale}; % first scale                
        
        %% FOV
        ROI_FOV = false(size(I));
        ROI_FOV(roibby0:roibby1,roibbx0:roibbx1) = true;
        ROI_FOV_Scale = imresize(ROI_FOV,curScale,'nearest');
        
        IlbpFov = Ilbp(ROI_FOV_Scale);
        lbpDescFov = hist(IlbpFov(:),0:9);
        lbpDescFov = lbpDescFov ./ sum(lbpDescFov);
        
        %         accumulatedFovLBP = [accumulatedFovLBP,lbpDescFov'];
        lbpData.fov{icell}.lbp(curT,:) = lbpDescFov;
        
        %% FWD
        ROI_FWD = false(size(I));
        ROI_FWD(roibby0:roibby1,roibbx0:roibbx1) = ROIBB;
        ROI_FWD_Scale = imresize(ROI_FWD,curScale,'nearest');
        
        IlbpFwd = Ilbp(ROI_FWD_Scale);
        lbpDescFwd = hist(IlbpFwd,0:9);
        lbpDescFwd = lbpDescFwd ./ sum(lbpDescFwd);
        
        %         accumulatedFwdLBP = [accumulatedFwdLBP,lbpDescFwd'];
        lbpData.fwd{icell}.lbp(curT,:) = lbpDescFwd;
        
        %% BCK
        ROI_BCK = false(size(I));
        ROI_BCK(roibby0:roibby1,roibbx0:roibbx1) = ~ROIBB;
        ROI_BCK_Scale = imresize(ROI_BCK,curScale,'nearest');
        
        IlbpBck = Ilbp(ROI_BCK_Scale);
        lbpDescBck = hist(IlbpBck,0:9);
        lbpDescBck = lbpDescBck ./ sum(lbpDescBck);
        
        %         accumulatedBckLBP = [accumulatedBckLBP,lbpDescBck'];
        lbpData.bck{icell}.lbp(curT,:) = lbpDescBck;         
    end      
    
    if ~isempty(lbpData)
        %% FOV - dLBP/dT
        dLbp = abs(lbpData.fov{icell}.lbp(1:end-1,:) - lbpData.fov{icell}.lbp(2:end,:));
        lbpData.fov{icell}.dlbp10d = dLbp;
        lbpData.fov{icell}.dlbp1d = sum(lbpData.fov{icell}.dlbp10d,2);
        %% FWD - dLBP/dT
        dLbp = abs(lbpData.fwd{icell}.lbp(1:end-1,:) - lbpData.fwd{icell}.lbp(2:end,:));
        lbpData.fwd{icell}.dlbp10d = dLbp;
        lbpData.fwd{icell}.dlbp1d = sum(lbpData.fwd{icell}.dlbp10d,2);
        %% BCK - dLBP/dT
        dLbp = abs(lbpData.bck{icell}.lbp(1:end-1,:) - lbpData.bck{icell}.lbp(2:end,:));
        lbpData.bck{icell}.dlbp10d = dLbp;
        lbpData.bck{icell}.dlbp1d = sum(lbpData.bck{icell}.dlbp10d,2);
        %% Correlations
        lbpData.corr{icell}.fov_fwd = corr(lbpData.fov{icell}.dlbp1d,lbpData.fwd{icell}.dlbp1d); % ~0
        lbpData.corr{icell}.fov_bck = corr(lbpData.fov{icell}.dlbp1d,lbpData.bck{icell}.dlbp1d); % > 0
        lbpData.corr{icell}.fwd_bck = corr(lbpData.fwd{icell}.dlbp1d,lbpData.bck{icell}.dlbp1d); % ~0
    end
end

%% Verifying that this is a non-empty task
nullanize = true;
for icell = 1 : nCells
    if ~isempty(lbpData.fov{icell})
        nullanize = false;
    end
end
if nullanize
    lbpData = [];
end
%%

% save(lbpFname,'lbpData','accumulatedFovLBP','accumulatedBckLBP','accumulatedFwdLBP');
save(lbpFname,'lbpData');
end