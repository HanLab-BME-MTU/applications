function [] = pcFineMatchingScoresMD(MD,params,dirs)

% Rosin thresholding of the cross-correlation-based matching scores for
% finding regions for matching at a finer resolution to localize the hits.

% Assaf Zaritsky, June 2015


for t = 1 : params.nTime - params.frameJump
    
    fprintf(sprintf('fine resolution frame %d\n',t));
    
    fineScoresFname = [dirs.fineResScores pad(t,3) '_fineScores.mat'];
    
    if exist(fineScoresFname,'file') && ~params.always
        continue;
    end
    
    mfFname = [dirs.mfData pad(t,3) '_mf.mat'];
    load(mfFname); % scores
    
    %     fineScores = scores;
    
    % Threhold the matching score
    scores1 = scores(~isnan(scores));
    scores2 = 1 - scores1;
    percentile2 = prctile(scores2(:),2);
    percentile98 = prctile(scores2(:),98);    
    centers = percentile2:(percentile98-percentile2)/50:percentile98;
    %     centers = 0 : 0.00001 : 0.001;
    [nelements, centers] = hist(scores2,centers);
    [~,thMatching] =  cutFirstHistMode(nelements,centers); % rosin threshold
    scoresNoNans = inpaint_nans(scores,4);
    BIN_SCORES = scoresNoNans < 1 - thMatching;
    BIN_SCORES = imdilate(BIN_SCORES,strel('square',params.patchSize));% added by Assaf 2015.06.23
    
    clear scores scoresNoNans;
    
    %     imgFname0 = [dirs.images pad(t,3) '.tif'];
    %     imgFname1 = [dirs.images pad(t+params.frameJump,3) '.tif'];
    %     I0 = imread(imgFname0);
    %     I1 = imread(imgFname1);
    I0 = MD.getChannel(1).loadImage(t);
    I1 = MD.getChannel(1).loadImage(t + params.frameJump);
    
    %     [dydx, dys, dxs, scores] = blockMatching(I0, I1, params.fineResolution,params.fineSearchRadius,BIN_SCORES);
    [dydx, dys, dxs, fineScores] = blockMatching(I0, I1, params.fineResolution,params.fineSearchRadius,BIN_SCORES);    
    
    %     fineScores(~isnan(scores)) = scores(~isnan(scores));
    %     inpaint_nans(fineScores); % added by Assaf 2015.06.23
    %     fineScores((BIN_SCORES & isnan(scores)) | isnan(fineScores)) = 1; % this cause problems, better in_paint nans.
    
    save(fineScoresFname,'fineScores','BIN_SCORES','thMatching');
    
    h = figure('visible','off');
    imagescnan(fineScores); title(sprintf('frame %d match score',t));
    caxis([percentile2,percentile98]); colorbar;
    outputFile = [dirs.fineResScoresVis pad(t,3) '_fineScores.jpg'];
    saveas(h,outputFile);
    %     eval(sprintf('print -djpeg %s', outputFile));
    close all;
end
end