function [] = whVectorFieldsVisualization(params,dirs)

time = 1 : 91 - params.frameJump - 1; % -1 because segmentation is based on motion

for t = time
    outFname = [dirs.mfVis filesep sprintf('%.3d',t) '_mfVis.eps'];
    
    if exist(outFname,'file')
        fprintf(sprintf('fetching motion visualization frame %d\n',t));
        return;
    end
    
    fprintf(sprintf('fetching motion visualization frame %d\n',t));
    
    load([dirs.mfData sprintf('%.3d',t) '_mf.mat']); % dxs, dys
    load([dirs.roiData sprintf('%.3d',t) '_roi.mat']); % ROI
    I = imread([dirs.images sprintf('%.3d',t) '.tif']); % images data
        
    %% Visualize motion flows and stretching events!    
    [h] = visualizeMotionFields(I,ROI,dxs,dys,params.patchSize,0.5);% last parameter is reduced resolution
    haxes = get(h,'CurrentAxes');
    set(haxes,'XTick',[]);
    set(haxes,'YTick',[]);
    set(haxes,'XTickLabel',[]);
    set(haxes,'YTickLabel',[]);
    set(h,'Color','w');
    axis equal; % no image stretching
    % set(haxes,'box','off');
    axis off;
    axis tight;
    hold off;
    export_fig(outFname);  
    close all;    
end

end