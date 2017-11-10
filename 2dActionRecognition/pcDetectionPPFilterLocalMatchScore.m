
function [] = pcDetectionPPFilterLocalMatchScore(MD,params,dirs)

% pcDetectionPPFilterLocalMatchScore
% Input: detected cells, processed cross-corrlation matching scores
% (include phase contrast) for forground and background
% Output: filter FP based on 
%   (1) matching score forward- backward <= 0.1 (markerd in 'r' for debug)
%   (2) huge size of connected components (never really happens)

% Assaf Zaritsky, June 2015
if params.always
    unix(sprintf('rm %s*_filtered.mat',[dirs.detectPPData filesep]));
    unix(sprintf('rm %s*_filtered.eps',[dirs.detectPPVis filesep]));
    fprintf('detection PP filter (always): clean output directories\n'); 
end
   
for t = 1 : params.nTime - params.frameJump
    PPFilteredFname = [dirs.detectPPData sprintf('%03d',t) '_filtered.mat'];
    PPFilteredVisFname = [dirs.detectPPVis sprintf('%03d',t) '_filtered.eps'];
    
    fprintf(sprintf('detection filtering frame %d\n',t));            
    
    if exist(PPFilteredFname,'file') && ~params.always
        if params.deepDebug
            fprintf('detection PP filter: file exist, continue\n');
        end;
        
        continue;
    end
        
    detectPPFname = [dirs.detectPPData sprintf('%03d',t) '_detectionStats.mat'];
    
    load(detectPPFname);% 'detections','combinedImage'                          
        
    % filter based on detection-background differences
        
    detections.yFilteredBackgroundLowRes = detections.yLowRes(detections.LowRes.diffScore <= params.diffDetectionBackgroundScoreTH);
    detections.xFilteredBackgroundLowRes = detections.xLowRes(detections.LowRes.diffScore <= params.diffDetectionBackgroundScoreTH);
            
    detections.yFilteredBackground = (detections.yFilteredBackgroundLowRes-0.5) .* params.patchSize;
    detections.xFilteredBackground = (detections.xFilteredBackgroundLowRes-0.5) .* params.patchSize;
    
    detections.yLowRes = detections.yLowRes(detections.LowRes.diffScore > params.diffDetectionBackgroundScoreTH);
    detections.xLowRes = detections.xLowRes(detections.LowRes.diffScore > params.diffDetectionBackgroundScoreTH);
    
    %     % filter by size
    %     [mask] = otsuMultTH(combinedImage,2);
    %     [L,nL] = bwlabel(mask==2,4);
    %
    %     filterMask = false(size(combinedImage));
    %     for l = 1 : nL
    %         curMask = L == l;
    %         if sum(curMask) * (params.patchSize*params.pixelSize).^2 > 80 * 30; % 80 x 30 um
    %             filterMask = filterMask | curMask;
    %         end
    %     end
    %
    %     n = length(detections.xLowRes);
    %     filtered = false(1,n);
    %     if sum(filterMask(:)) > 0
    %         for d = 1 : n
    %             if filterMask(round(detections.yLowRes(d)),round(detections.xLowRes(d)))
    %                 filtered(d) = true;
    %             end
    %         end
    %     end
    %
    %     detections.yFilteredSizeLowRes = detections.yLowRes(filtered);
    %     detections.xFilteredSizeLowRes = detections.xLowRes(filtered);
    %
    %     detections.yFilteredSize = (detections.yFilteredSizeLowRes-0.5) .* params.patchSize;
    %     detections.xFilteredSize = (detections.xFilteredSizeLowRes-0.5) .* params.patchSize;
    %
    %     detections.yLowRes = detections.yLowRes(~filtered);
    %     detections.xLowRes = detections.xLowRes(~filtered);
       
    detections.y = (detections.yLowRes-0.5) .* params.patchSize;
    detections.x = (detections.xLowRes-0.5) .* params.patchSize;

    I = MD.getChannel(1).loadImage(t);
    
    if params.deepDebug
            fprintf('detection PP filter: done loading\n');
    end;
    
     % overlay detections on original image
    h = figure();
    imagesc(I); colormap(gray);    
    hold on
    plot(detections.x,detections.y,'g+','MarkerSize',5,'LineWidth',1.5);
    plot(detections.xFiltered,detections.yFiltered,'m+','MarkerSize',5,'LineWidth',1.5);
    if ~isempty(detections.xFilteredBackground)
        plot(detections.xFilteredBackground,detections.yFilteredBackground,'r+','MarkerSize',5,'LineWidth',1.5);
    end
    %     if ~isempty(detections.xFilteredSize)
    %         plot(detections.xFilteredSize,detections.yFilteredSize,'y+','MarkerSize',5,'LineWidth',1.5);
    %     end
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
    export_fig_biohpc(PPFilteredVisFname);
    
    close all;
    
    
    save(PPFilteredFname,'detections','combinedImage','pstruct');
    
    if params.deepDebug
            fprintf('detection PP filter: done\n');
    end;
end
end