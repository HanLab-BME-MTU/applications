%% pcLBPMD - calculates LBP for single cell trajectories
% calculate LBP for image, then crops specific cell and extract
% corresponding LBP values

function [] = pcLBPMD(MD,params,dirs)

lbpMapping = getmapping(8,'riu2');

load([dirs.tracking 'cellIdTYX.mat']);% cellTYX

nCells = length(cellTYX);

accumulatedFovLbp = []; % take LBP for all field of view
lbpData.fov = cell(1,nCells);

lbpFname = [dirs.tracking filesep 'lbpData.mat'];

if exist(lbpFname,'file') && ~params.always
    fprintf(sprintf('LBP for %s exists, finishing\n',lbpFname));
    return;
end

%% LBP for frames
for t = 1 : params.nTime - params.frameJump
    lbpTFname = [dirs.lbp sprintf('%03d',t) '_lbp.mat'];       
    
    % Just to correct a bug
    if exist([dirs.lbp sprintf('%03d',t) '_mf.mat'],'file')
        unix(sprintf('mv %s %s',[dirs.lbp sprintf('%03d',t) '_mf.mat'],lbpTFname));
    end
    
    if exist(lbpTFname,'file') && ~params.always
        if params.deepDebug
            fprintf(sprintf('LBP frame %d: continue\n',t)); 
        end;
        continue;
    end
    
    fprintf(sprintf('LBP frame %d\n',t));
    
    I = MD.getChannel(1).loadImage(t);
    IlbpTmp = lbp(I,1,8,lbpMapping,'');
    Ilbp = nan(size(I));
    Ilbp(2:end-1,2:end-1) = IlbpTmp;
    
    save(lbpTFname,'Ilbp');    
end

%% LBP for cells
for icell = 1 : nCells
    fprintf(sprintf('LBP cell %d/%d\n',icell,nCells));
    curCell = cellTYX{icell};
    ntime = length(curCell.ts);
    t0 = curCell.ts(1);
    lbpData.fov{icell}.lbp = nan(ntime,10);
    for t = curCell.ts
        curT = t - t0 + 1;                
        
        bby0 = max(1,curCell.ys(curT) - params.FOVRadius);
        bby1 = min(size(I,1),curCell.ys(curT) + params.FOVRadius);
        bbx0 = max(1,curCell.xs(curT) - params.FOVRadius);
        bbx1 = min(size(I,2),curCell.xs(curT) + params.FOVRadius);
                
        lbpTFname = [dirs.lbp sprintf('%03d',t) '_lbp.mat'];       
        load(lbpTFname); % Ilbp
        IlbpBB = Ilbp(bby0:bby1,bbx0:bbx1);
                
        %         % if segmentation
        %         lbpBin = lbpBB(cells.seg{icells});
        
        lbpDesc = hist(IlbpBB(:),0:9);
        lbpDesc = lbpDesc ./ sum(lbpDesc);
        accumulatedFovLbp = [accumulatedFovLbp,lbpDesc'];
        lbpData.fov{icell}.lbp(curT,:) = lbpDesc;
    end    
end
save(lbpFname,'lbpData','accumulatedFovLbp');
end