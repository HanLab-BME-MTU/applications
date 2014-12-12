function [] = whLocalMotionEstimation(params,dirs)

for t = 1 : params.nTime - params.frameJump
    mfFname = [dirs.mfData pad(t,3) '_mf.mat'];
    
    if exist(mfFname,'file') && ~params.always
        continue;
    end
        
    fprintf(sprintf('motion estimation frame %d\n',t));
    imgFname0 = [dirs.images pad(t,3) '.tif'];
    imgFname1 = [dirs.images pad(t+params.frameJump,3) '.tif'];
    I0 = imread(imgFname0);
    I1 = imread(imgFname1);
    
    [dydx, dys, dxs, scores] = blockMatching(I0, I1, params.patchSize,params.searchRadiusInPixels,true(size(I0))); % block width, search radius,
    
    if params.fixGlobalMotion
        meanDxs = mean(dxs(~isnan(dxs)));
        meanDys = mean(dxs(~isnan(dxs)));
        if abs(meanDxs) > 0.5
            dxs = dxs - meanDxs;
        end
        if abs(meanDys) > 0.5
            dys = dys - meanDys;
        end
    end
    
    save(mfFname,'dxs', 'dys','scores');
    
    figure;
    imagesc(scores); title(sprintf('frame %d match score',t));
    caxis([0.995,1]); colorbar;
    outputFile = [dirs.mfScores pad(t,3) '_score.jpg'];
    eval(sprintf('print -djpeg %s', outputFile));
    close all;
end
end