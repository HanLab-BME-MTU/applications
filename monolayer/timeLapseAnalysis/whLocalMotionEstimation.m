function [] = whLocalMotionEstimation(params,dirs)

if exist([dirs.mfDataOrig filesep '001_mf.mat'],'file') && params.always   
    % unix(sprintf('rm %s',[dirs.mfDataOrig '*.mat']));
    delete([dirs.mfDataOrig '*.mat']);    
end

if exist([dirs.mfData filesep '001_mf.mat'],'file') && params.always    
    % unix(sprintf('rm %s',[dirs.mfData '*.mat']));    
    delete([dirs.mfData '*.mat']);
end

for t = 1 : params.nTime - params.frameJump
    mfFname = [dirs.mfData sprintf('%03d',t) '_mf.mat'];
    
    if exist(mfFname,'file') && ~params.always
        fprintf(sprintf('fetching motion estimation frame %d\n',t));
        continue;
    end
        
    fprintf(sprintf('motion estimation frame %d\n',t));
    imgFname0 = [dirs.images sprintf('%03d',t) '.tif'];
    imgFname1 = [dirs.images sprintf('%03d',t+params.frameJump) '.tif'];
    I0 = imread(imgFname0);
    
    if ~exist(imgFname1,'file') % create from previous 
        % unix(sprintf('cp %s %s',imgFname0,imgFname1));
        copyfile(imgFname0, imgFname1);
    end
    
    I1 = imread(imgFname1);
    
    % Assumeing 3 same channels
    if size(I0,3) > 1
        tmp = I0(:,:,1) - I0(:,:,2);
        assert(sum(tmp(:)) == 0);
        I0 = I0(:,:,1);
        I1 = I1(:,:,1);
    end
    
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
    
    %     eval(sprintf('print -djpeg %s', outputFile));
    % if isunix
    %     outputFile = [dirs.mfScores sprintf('%03d',t) '_score.eps'];
    %     export_fig_biohpc(outputFile);
    % else
        disp('Exporting figure via MATLAB print -dpdf (instead of export_fig)');
        outputFile = [dirs.mfScores sprintf('%03d',t) '_score'];
        print(outputFile, '-dpdf')
    % end
    close all;
end
end