
function [] = whLocalMotionEstimationMD(MD,params,dirs)

% PIV-based local motion estimation. 
% Outputs for each frame the motion vectors dxs, dys and the maximal cross
% correlation score to be later used for segmentation / detection.

% Assaf Zaritsky, June 2015


if params.always
    unix(sprintf('rm %s*_mf.mat',[dirs.mfData filesep]));
    unix(sprintf('rm %s*.jpg',[dirs.mfScores filesep]));
    unix(sprintf('rm %s*.eps',[dirs.mfScores filesep]));
    fprintf('motion estimation (always): clean output directories\n'); 
end

for t = 1 : params.nTime - params.frameJump
    mfFname = [dirs.mfData sprintf('%03d',t) '_mf.mat'];
    
    fprintf(sprintf('motion estimation frame %d\n',t));
    
    if exist(mfFname,'file') && ~params.always
        if params.deepDebug
            fprintf('motion estimation: continue\n'); 
        end;
        continue;
    end
            
    %     imgFname0 = [dirs.images pad(t,3) '.tif'];
    %     imgFname1 = [dirs.images pad(t+params.frameJump,3) '.tif'];
    %     I0 = imread(imgFname0);
    %     I1 = imread(imgFname1);
    I0 = MD.getChannel(1).loadImage(t);
    I1 = MD.getChannel(1).loadImage(t + params.frameJump);
    
    if params.deepDebug
            fprintf('motion estimation: before block matching\n'); 
    end;
    
    [dydx, dys, dxs, scores] = blockMatching(I0, I1, params.patchSize,params.searchRadiusInPixels,true(size(I0))); % block width, search radius,
    
    if params.deepDebug
            fprintf('motion estimation: after block matching\n'); 
    end;
    
    %     if params.fixGlobalMotion
    %         meanDxs = mean(dxs(~isnan(dxs)));
    %         meanDys = mean(dxs(~isnan(dxs)));
    %         if abs(meanDxs) > 0.5
    %             dxs = dxs - meanDxs;
    %         end
    %         if abs(meanDys) > 0.5
    %             dys = dys - meanDys;
    %         end
    %     end
        
    
    if params.deepDebug
            fprintf('motion estimation: before figures\n'); 
    end;
    
    figure;
    imagesc(scores); title(sprintf('frame %d match score',t));
    %     caxis([prctile(scores(:),2),prctile(scores(:),98)]); colorbar;
    outputFile = [dirs.mfScores sprintf('%03d',t) '_score.jpg'];
    eval(sprintf('print -djpeg %s', outputFile));    
    
    % Illumination correction??
    maxY = 120;
    x_ = mean(I0,1); x_1 = min(max(x_ - prctile(x_,10),0),maxY);
    y_ = mean(I0,2); y_1 = min(max(y_ - prctile(y_,10),0),maxY);
    figure; plot(1:size(I0,2),x_1); ylim([0,maxY]);  
    title('X','FontSize',32);
    export_fig_biohpc([dirs.mfScores sprintf('%03d',t) '_illumX.eps']);
    %     eval(sprintf('print -djpeg %s', [dirs.mfScores pad(t,3) '_illumX.jpg']));
    
    figure; plot(1:size(I0,1),y_1); ylim([0,maxY]);
    title('Y','FontSize',32);
    export_fig_biohpc([dirs.mfScores sprintf('%03d',t) '_illumY.eps']);
    %     eval(sprintf('print -djpeg %s', [dirs.mfScores sprintf('%03d',t) '_illumY.eps']));
    close all;        
    
    save(mfFname,'dxs', 'dys','scores','x_','y_');
    
    if params.deepDebug
            fprintf('motion estimation: after figures\n'); 
    end;
end
end