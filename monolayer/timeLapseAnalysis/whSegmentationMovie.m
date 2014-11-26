function [] = whSegmentationMovie(params,dirs)
close all;
fprintf('start segmentation movie\n')

segmentationFname = [dirs.segmentation dirs.expname '_segmentation.avi'];

aviFname = [dirs.roiVis dirs.expname '_segmentation.avi'];
aviobj = avifile(segmentationFname,'fps',3,'compression','None');

for t = 1 : params.nTime - params.frameJump
    
    load([dirs.roiData pad(t,3) '_roi.mat']); % ROI
    I = imread([dirs.images pad(t,3) '.tif']);
    perimWidth = round(max(size(I)) / 200);
    I(dilate(bwperim(ROI,8),perimWidth)) = max(255,max(I(:)));
    
    h = figure('visible','off'); imagesc(I); colormap(gray);
    text(size(I,1)-500,size(I,2)-100,sprintf('%d minutes',round(t*params.timePerFrame)),'color','w','FontSize',15);
    haxes = get(h,'CurrentAxes');
    set(haxes,'XTick',[]);
    set(haxes,'XTickLabel',[]);
    set(haxes,'YTick',[]);
    set(haxes,'YTickLabel',[]);
    %     drawnow;
    movieFrame = getframe(h);
    aviobj = addframe(aviobj, movieFrame);
    close all;
end
aviobj = close(aviobj);
fprintf('start segmentation movie\n')
end