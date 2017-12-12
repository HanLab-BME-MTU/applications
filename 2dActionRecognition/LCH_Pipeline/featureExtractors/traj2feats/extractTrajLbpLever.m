%% extractTrajLbpLever - calculates LBP with LEVER's segmentation mask

function curTrajFeats = extractTrajLbpLever(curCell,MD,params,dirs)

% curFOVRadius = round(params.FOVRadius);

ntime = length(curCell.ts);
t0 = curCell.ts(1);

curTrajFeats.xs = curCell.xs; 
curTrajFeats.ys = curCell.ys;
curTrajFeats.ts = curCell.ts;
curTrajFeats.icell = curCell.icell;

curTrajFeats.feats = [];

leverFname = [dirs.roiLever filesep num2str(curCell.icell) '_roi.mat']; 
if ~exist(leverFname,'file')
    return; % empty feature - lever segmentation does not exist for this cell
end
load(leverFname); % curCellRoi{curT}.roi

curTrajFeats.rawFeats = nan(10,length(curTrajFeats.ts));

for t = curCell.ts
    curT = t - t0 + 1;
    
    LEVER_ROI = curCellRoi{curT}.roi;
    
    if sum(LEVER_ROI(:) == 0)
        curTrajFeats.feats = [];
        return;
    end        
    
    curX = round(curCell.xs(curT));
    curY = round(curCell.ys(curT));
    
    lbpTFname = [dirs.lbp sprintf('%03d',t) '_lbp.mat'];
    load(lbpFname); % pyramidLbp
    Ilbp = pyramidLbp{1};    
        
    %     bby0 = max(1,curY - curFOVRadius);
    %     bby1 = min(size(Ilbp,1),curY + curFOVRadius);
    %     bbx0 = max(1,curX - curFOVRadius);
    %     bbx1 = min(size(Ilbp,2),curX + curFOVRadius);
    %
    %     IlbpBB = Ilbp(bby0:bby1,bbx0:bbx1);
    
    %         % if segmentation
    %         lbpBin = lbpBB(cells.seg{icells});
    IlbpBB_Lever = Ilbp(LEVER_ROI);    
    
    lbpLeverDesc = hist(IlbpBB_Lever(:),0:9);
    lbpLeverDesc = lbpLeverDesc ./ sum(lbpLeverDesc);
    
    curTrajFeats.rawFeats(:,curT) = lbpLeverDesc;      
end
curTrajFeats.feats = ; % check how it is done...
end