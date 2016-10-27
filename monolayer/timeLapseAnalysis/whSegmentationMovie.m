function [] = whSegmentationMovie(params,dirs)
close all;
fprintf('start segmentation movie\n')

segmentationFname = [dirs.segmentation dirs.expname '_segmentation.avi'];

if exist(segmentationFname,'file') && ~params.always
    return;
end

vwriter = VideoWriter(segmentationFname);
vwriter.FrameRate = 3;
open(vwriter);

% % aviFname = [dirs.roiVis dirs.expname '_segmentation.avi'];
% aviobj = avifile(segmentationFname,'fps',3,'compression','None');
W = nan; H = nan;
for t = 1 : params.nTime - params.frameJump
    
    load([dirs.roiData pad(t,3) '_roi.mat']); % ROI
    I = imread([dirs.images pad(t,3) '.tif']);
    
    % Assumeing 3 same channels
    if size(I,3) > 1
        tmp = I(:,:,1) - I(:,:,2);
        assert(sum(tmp(:)) == 0);
        I = I(:,:,1);
    end
    
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
    drawnow; pause(0.01);
    movieFrame = getframe(h);
    %     aviobj = addframe(aviobj, movieFrame);
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
end
% aviobj = close(aviobj);
close(vwriter);
fprintf('finish segmentation movie\n')
end

%% UTILS
function [Aer] = erode(A,maskSize)

if nargin < 2
    error('whTemporalBasedSegmentation: erode missing mask size');
end

mask = ones(maskSize);
se1 = strel(mask);

Aer = imerode(A,se1);

end

function [Aer] = dilate(A,maskSize)

if nargin < 2
    error('whTemporalBasedSegmentation: dilate missing mask size');
end

mask = ones(maskSize);
se1 = strel(mask);

Aer = imdilate(A,se1);

end
