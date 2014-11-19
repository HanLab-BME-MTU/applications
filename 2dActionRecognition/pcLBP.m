function [] = pcLBP(params,dirs)

lbpMapping = getmapping(8,'riu2');

accumulatedLBP = []; % TODO: efficient implementation
accumulatedCellsID = [];

for t = 2 : params.nTime - params.frameJump - 1
    
    fprintf(sprintf('LBP %d\n',t));
    
    
    lbpFname = [dirs.lbp pad(t,3) '_lbp.mat']; 
    
    if exist(lbpFname,'file') && ~params.always
        load(lbpFname);
    else
        
        I = imread([dirs.images pad(t,3) '.tif']);
        
        Ilbp0 = lbp(I,1,8,lbpMapping,'');
        Ilbp = nan(size(I));
        Ilbp(2:end-1,2:end-1) = Ilbp0;
        
        cellsFname = [dirs.detectData pad(t,3) '_cells.mat'];
        load(cellsFname);
        ncells = length(cells.xs);
        
        accumulatedLbpFrame = [];
        accumulatedCellIDFrame = [];
        for icells = 1 : ncells
            lbpBB = Ilbp(cells.bb{icells}.bbys,cells.bb{icells}.bbxs);
            lbpBin = lbpBB(cells.seg{icells});
            lbpDesc = hist(lbpBin,0:9);
            lbpDesc = lbpDesc ./ sum(lbpDesc);
            cells.lbp{icells} = lbpDesc;
            accumulatedLbpFrame = [accumulatedLbpFrame,lbpDesc'];
            accumulatedCellIDFrame = [accumulatedCellIDFrame,[t,icells]']; % (time, #cell)
        end
        
        save(lbpFname,'Ilbp','cells','accumulatedLbpFrame','accumulatedCellIDFrame');
    end
    accumulatedLBP = [accumulatedLBP, accumulatedLbpFrame];
    accumulatedCellsID = [accumulatedCellsID, accumulatedCellIDFrame];
end

save([dirs.results dirs.expname '_lbpPerFrame.mat'],'accumulatedLBP','accumulatedCellsID');

end