
function [] = pcDetectionMovie(MD,params,dirs)

% Save detection movies to /project/cellbiology/gdanuser/melanomaModel/Analysis/Movies/detectionMovies
% Note: over 2 minutes per frame from getFrame!!

% Assaf Zaritsky, November 2015

outdir = '/project/cellbiology/gdanuser/december/assaf/LCH_Analysis/Movies/detectionMovies/';
% '/project/cellbiology/gdanuser/melanomaModel/Analysis/Movies/detectionMovies/';

movieFname = [outdir dirs.expname '_detection.avi'];

if exist(movieFname,'file') && ~params.always
    fprintf(sprintf('Detection movie %s exists, finishing\n',movieFname));
    return;
end

vwriter = VideoWriter(movieFname,'Uncompressed AVI');
vwriter.FrameRate = 10;
open(vwriter);

W = nan; H = nan;
for t = 1 : params.nTime - params.frameJump
    detectPPFname = [dirs.detectPPData sprintf('%03d',t) '_backtrack.mat'];
    
    if ~exist(detectPPFname,'file')
        continue;
    end
    
    I = MD.getChannel(1).loadImage(t);
    
    load(detectPPFname);%'detections','combinedImage'
    
    h = figure('visible','off');
    subplot(1,2,1)
    imagesc(combinedImage); colormap(gray);
    hold on;
    plot(detections.xLowRes,detections.yLowRes,'go','MarkerSize',3,'LineWidth',1);
    plot(detections.xBacktrackLowRes,detections.yBacktrackLowRes,'co','MarkerSize',3,'LineWidth',1);
    haxes = get(h,'CurrentAxes');
    set(haxes,'XTick',[]);
    set(haxes,'XTickLabel',[]);
    set(haxes,'YTick',[]);
    set(haxes,'YTickLabel',[]);
    axis equal; axis off;
    hold off;
    
    subplot(1,2,2)
    imagesc(I); colormap(gray);
    hold on;
    plot(detections.x,detections.y,'go','MarkerSize',3,'LineWidth',1);
    plot(detections.xBacktrack,detections.yBacktrack,'co','MarkerSize',3,'LineWidth',1);
    text(size(I,1)-1000,size(I,2)-550,sprintf('%d minutes',round(t*params.timePerFrame)),'color','w','FontSize',15);
    haxes = get(h,'CurrentAxes');
    set(haxes,'XTick',[]);
    set(haxes,'XTickLabel',[]);
    set(haxes,'YTick',[]);
    set(haxes,'YTickLabel',[]);
    axis equal; axis off;
    
    hold off;
    
    drawnow;
    pause(0.05);
    
    movieFrame = getframe(h);
    
    if isnan(W)
        [H,W,~] = size(movieFrame.cdata);
        minH = H;
        maxH = H;
        minW = W;
        maxW = W;
    end
    
    if H ~= size(movieFrame.cdata,1) || W ~= size(movieFrame.cdata,2)
        minH = min(H,size(movieFrame.cdata,1));
        maxH = max(H,size(movieFrame.cdata,2));
        minW = min(W,size(movieFrame.cdata,1));
        maxW = max(W,size(movieFrame.cdata,2));
    end
    
    movieFrameResized = uint8(zeros(H,W,3));
    movieFrameResized(:,:,1) = imresize(movieFrame.cdata(:,:,1),[H,W]);
    movieFrameResized(:,:,2) = imresize(movieFrame.cdata(:,:,2),[H,W]);
    movieFrameResized(:,:,3) = imresize(movieFrame.cdata(:,:,3),[H,W]);
    movieFrame.cdata = movieFrameResized;
    
    writeVideo(vwriter,movieFrame);
    close all;
    fprintf(sprintf('frame %d\n',t));
end

close(vwriter);

fprintf(sprintf('Done creating detection movie (H: %d-%d, W:%d-%d)\n,',minH,maxH,minW,maxW));
end