function [] = pcMorphLocalSpeed(params,dirs)

accumulatedSizes = [];
accumulatedEccentricity = [];
accumulatedLocalSpeed = [];
accumulatedLocalMatchingScore = [];
accumulatedCellsID = [];

for t = 2 : params.nTime - params.frameJump - 1
    
    fprintf(sprintf('MorphLocalSpeed %d\n',t));
    
    
    localMorphDynamicsFname = [dirs.localMorphDynam pad(t,3) '_localMorphDynam.mat']; % TODO: create dirs.localMorphDynam
    
    if exist(localMorphDynamicsFname,'file') && ~params.always
        load(localMorphDynamicsFname);
    else
        load([dirs.mfData pad(t,3) '_mf.mat']); % dxs,dys
        speed = sqrt(dxs.^2+dys.^2);
        load([dirs.fineResScores pad(t,3) '_fineScores.mat']); % fineScores, BIN_SCORES
        matchScores = nan(size(speed));
        matchScores(BIN_SCORES) = fineScores(BIN_SCORES); % nanmean
        
        
        cellsFname = [dirs.detectData pad(t,3) '_cells.mat'];
        load(cellsFname); % cells
        ncells = length(cells.xs);
        
        accumulatedSizesFrame = [];
        accumulatedEccentricityFrame = [];
        accumulatedLocalSpeedFrame = [];
        accumulatedLocalMatchingScoreFrame = [];
        accumulatedCellIDFrame = [];
        for icells = 1 : ncells
            roi = cells.seg{icells};
            stats = regionprops(roi, 'Area', 'Eccentricity');
            
            ROI = false(size(I));
            ROI(cells.bb{icells}.bbys,cells.bb{icells}.bbxs) = roi;
            accumulatedSizesFrame = [accumulatedSizesFrame stats.Area];
            accumulatedEccentricityFrame = [accumulatedEccentricityFrame stats.Eccentricity];
            accumulatedLocalSpeedFrame = [accumulatedLocalSpeedFrame nanmean(speed(ROI))];
            accumulatedLocalMatchingScoreFrame = [accumulatedLocalMatchingScoreFrame nanmean(matchScores(ROI))];
            accumulatedCellIDFrame = [accumulatedCellIDFrame,[t,icells]']; % (time, #cell)
        end
        
        save(localMorphDynamicsFname,'accumulatedSizesFrame','accumulatedEccentricityFrame',...
            'accumulatedLocalSpeedFrame','accumulatedLocalMatchingScoreFrame','accumulatedCellIDFrame');
    end
    accumulatedSizes = [accumulatedSizes, accumulatedSizesFrame];
    accumulatedEccentricity = [accumulatedEccentricity, accumulatedEccentricityFrame];
    accumulatedLocalSpeed = [accumulatedLocalSpeed, accumulatedLocalSpeedFrame];
    accumulatedLocalMatchingScore = [accumulatedLocalMatchingScore, accumulatedLocalMatchingScoreFrame];
    accumulatedCellsID = [accumulatedCellsID, accumulatedCellIDFrame];
    
end

save([dirs.results dirs.expname '_localMorphDynam.mat'],'accumulatedSizes','accumulatedEccentricity',...
    'accumulatedLocalSpeed','accumulatedLocalMatchingScore','accumulatedCellsID');

end