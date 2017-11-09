function [] = pcDetectMotionEvents(params,dirs)

windowSize = round(10 / (params.pixelSize * params.fineResolution)); % 10 um
G = fspecial('gaussian',[windowSize windowSize],2);
maskLocalMax = round(40 / (params.pixelSize * params.fineResolution)); % 40 um
for t = 2 : params.nTime - params.frameJump - 1
    
    fprintf(sprintf('detect motion events frame %d\n',t));
    
    detectFname = [dirs.detectData pad(t,3) '_detections.mat'];
    
    if exist(detectFname,'file') && ~params.always
        continue;
    end
    
    fineScoresFname0 = [dirs.fineResScores pad(t-1,3) '_fineScores.mat'];
    fineScoresFname1 = [dirs.fineResScores pad(t,3) '_fineScores.mat'];
%     fineScoresFname2 = [dirs.fineResScores pad(t+1,3) '_fineScores.mat'];
    
    load(fineScoresFname0); % fineScores,BIN_SCORES,thMatching
    fineScores0 = fineScores;
    load(fineScoresFname1);
    fineScores1 = fineScores;
    mask = BIN_SCORES;
    %     load(fineScoresFname2);
    %     fineScores2 = fineScores;
    %     combinedScores = min(fineScores0 + fineScores1 + fineScores2,fineScores1*3); % why take the min?
    combinedScores = fineScores0 + fineScores1;
    % fix 2 to inpaint
    combinedScores(combinedScores == 2) = nan;
    combinedScores = inpaint_nans(combinedScores,4);
    
    % Local maxima detection and filter by BIN_SCORES
    maskLowRes = imresize(mask,1/params.fineResolution);
    scoresLowRes = imresize(combinedScores,1/params.fineResolution);
    scoresLowResFilter = imfilter(scoresLowRes,G,'same');
    %     scoresLowResFilter = inpaint_nans(scoresLowResFilter,4);
    %     scoresLowResFilter(isnan(scoresLowResFilter)) = 2;%3;
    keepFlat = 1;
    socresLocalMaxima = locmax2d(2-scoresLowResFilter, maskLocalMax, keepFlat);%3
    socresLocalMaxima(dilate(isnan(scoresLowResFilter),3)) = nan;
    socresLocalMaxima(~maskLowRes) = nan;
    [X,Y] = meshgrid(1:size(socresLocalMaxima,2),1:size(socresLocalMaxima,1));
    detectionsLowRes.xs = X(socresLocalMaxima > 0);
    detectionsLowRes.ys = Y(socresLocalMaxima > 0);
    
    detections.xs = detectionsLowRes.xs * params.fineResolution-floor(params.fineResolution/2);
    detections.ys = detectionsLowRes.ys * params.fineResolution-floor(params.fineResolution/2);
    
    I = imread([dirs.images pad(t,3) '.tif']);
    I(dilate(bwperim(mask),5)) = max(255,max(I(:)));
    h = figure('visible','off');
    imagesc(I); colormap(gray);    
    hold on;
    plot(detections.xs,detections.ys,'sw','MarkerFaceColor','w','MarkerSize',3);
    haxes = get(h,'CurrentAxes');
    set(haxes,'XTick',[]);
    set(haxes,'XTickLabel',[]);
    set(haxes,'YTick',[]);
    set(haxes,'YTickLabel',[]);
    hold off;
    saveas(h,[dirs.detectVis pad(t,3) '_detection1.jpg']);    
    close(h);
       
    
    h = figure('visible','off');    
    combinedScoresDisplay = combinedScores;    
    combinedScoresDisplay(dilate(bwperim(mask),5)) = nan;    
    %     imagescnan(combinedScoresDisplay);
    imagesc(combinedScoresDisplay);
    caxis([1.999,max(prctile(combinedScores(~isnan(combinedScores)),1.9991),90)]);    
    hold on;
    plot(detections.xs,detections.ys,'sw','MarkerFaceColor','w','MarkerSize',6);    
    haxes = get(h,'CurrentAxes');
    set(haxes,'XTick',[]);
    set(haxes,'XTickLabel',[]);
    set(haxes,'YTick',[]);
    set(haxes,'YTickLabel',[]);    
    hold off;    
    saveas(h,[dirs.detectVis pad(t,3) '_detection2.jpg']);    
    close(h);
    
    save(detectFname,'I','detections','combinedScores','BIN_SCORES','socresLocalMaxima');
end
end