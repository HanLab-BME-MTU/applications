
function [] = pcDetectionMD(MD,params,dirs)

% Cell detection
% Input: raw image, cross-corrlation matching scores
% Output: detected cells (+ filtered)

% Algorithm:
%   1. Down sample images by patchSize
%   2. Gaussian filter images
%   3. Normalize images
%   4. Create combined image from both channels
%   5. pointSourceDetection

% Note: tends to over-detect (false-positives), which should be solved in
% the tracking

% Assaf Zaritsky, June 2015

if params.always
    unix(sprintf('rm %s*.mat',[dirs.detectData filesep]));
    unix(sprintf('rm %s*.eps',[dirs.detectVis filesep]));
    fprintf('cell detection (always): clean output directories\n'); 
end

pointSourceDetectionScale = 1.85; % optimized for our dataset
gaussianSigma = 0.6; % optimized for our dataset

gaussianFilter = fspecial('gaussian',[3 3], gaussianSigma);

for t = 1 : params.nTime - params.frameJump
    detectFname = [dirs.detectData sprintf('%03d',t) '_detections.mat'];
    combinedImageFname = [dirs.detectVis sprintf('%03d',t) '_combinedImage.eps'];
    detectVisFname = [dirs.detectVis sprintf('%03d',t) '_detections.eps'];
    mfFname = [dirs.mfData sprintf('%03d',t) '_mf.mat'];
       
    fprintf(sprintf('cell detection frame %d\n',t));
    
    if exist(detectFname,'file') && ~params.always        
        if params.deepDebug
            fprintf('cell detection: continue\n'); 
        end;
        
        continue;
    end            
    
    I0 = MD.getChannel(1).loadImage(t);
    load(mfFname); % scores
    
    if params.deepDebug
            fprintf('cell detection: before combined image\n'); 
    end;
    
    scoresInv = 1 - scores;
    resizedScore = imresize(scoresInv,1.0/params.patchSize);
    resizedScoreFiltered = imfilter(inpaint_nans(resizedScore), gaussianFilter, 'replicate');
    p02 = double(prctile(resizedScoreFiltered(:),0.2));
    p998 = prctile(resizedScoreFiltered(:),99.8);
    resizedScoreFilteredNorm = min(1,max(0,((double(resizedScoreFiltered - p02)) ./ double(p998-p02))));
    
    I0Resized = imresize(I0,1.0/params.patchSize);
    I0ResizedFiltered = imfilter(I0Resized, gaussianFilter, 'replicate');
    p02 = double(prctile(I0ResizedFiltered(:),0.2));
    p998 = prctile(I0ResizedFiltered(:),99.8);
    I0ResizedFilteredNorm = min(1,max(0,(double(I0ResizedFiltered - p02) ./ double(p998-p02))));
    
    combinedImage = (resizedScoreFilteredNorm + I0ResizedFilteredNorm)/2;
    [pstruct,mask,imgLM,imgLoG]=pointSourceDetection(combinedImage,pointSourceDetectionScale,'FitMixtures', false,'Alpha',0.05);
    
    if params.deepDebug
            fprintf('cell detection: after combined image\n'); 
    end;
    
    if isempty(pstruct)
        
        if params.deepDebug
            fprintf('cell detection: no detections (warning)\n');
        end;
        
        warning(['pstruct empty in pcDetectionMD frame ' num2str(t)]);
        
        detections.y = [];
        detections.x = [];
        
        detections.yFiltered = [];
        detections.xFiltered = [];
        
        detections.yLowRes = [];
        detections.xLowRes = [];
        
        detections.yFilteredLowRes = [];
        detections.xFilteredLowRes = [];
        
        save(detectFname,'detections','combinedImage','pstruct');        
        continue;
    end
    
    detections.y = (pstruct.y(~pstruct.isPSF)-0.5) .* params.patchSize;
    detections.x = (pstruct.x(~pstruct.isPSF)-0.5) .* params.patchSize;
        
    detections.yFiltered = (pstruct.y(pstruct.isPSF)-0.5) .* params.patchSize;
    detections.xFiltered = (pstruct.x(pstruct.isPSF)-0.5) .* params.patchSize;
    
    detections.yLowRes = pstruct.y(~pstruct.isPSF);
    detections.xLowRes = pstruct.x(~pstruct.isPSF);
        
    detections.yFilteredLowRes = pstruct.y(pstruct.isPSF);
    detections.xFilteredLowRes = pstruct.x(pstruct.isPSF);
    
    if params.deepDebug
            fprintf('cell detection: before outputs\n'); 
    end;
    
    %% outputs
    
    % image given as input to the detection algorithm
    h = figure;
    imagesc(combinedImage);
    colormap(gray);
    haxes = get(h,'CurrentAxes');
    set(haxes,'XTick',[]);
    set(haxes,'YTick',[]);
    set(haxes,'XTickLabel',[]);
    set(haxes,'YTickLabel',[]);
    set(h,'Color','w');
    axis equal;
    axis off;
    axis tight;
    export_fig_biohpc(combinedImageFname);
    
    % overlay detections on original image
    h = figure();
    imagesc(I0); colormap(gray);    
    hold on
    plot(detections.x,detections.y,'g+','MarkerSize',5,'LineWidth',2);
    plot(detections.xFiltered,detections.yFiltered,'m+','MarkerSize',5,'LineWidth',2);
    haxes = get(h,'CurrentAxes');
    set(haxes,'XTick',[]);
    set(haxes,'YTick',[]);
    set(haxes,'XTickLabel',[]);
    set(haxes,'YTickLabel',[]);
    set(h,'Color','w');
    axis equal;
    axis off;
    axis tight;
    hold off;
    export_fig_biohpc(detectVisFname);
    
    close all;
    
    save(detectFname,'detections','combinedImage','pstruct');
    
    if params.deepDebug
            fprintf('cell detection: after figures\n'); 
    end;
end
end