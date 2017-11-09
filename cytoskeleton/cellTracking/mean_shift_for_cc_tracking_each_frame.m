% function mean_shift_for_cc_tracking_each_frame(currentImg,iFrame, max_97_percent_int, bandwidth)
Mask  = currentImg>max_97_percent_int;

[inda,indb] = find(Mask>0);

ptRawData = [(indb) (inda)];

% Run mean-shift
% bandwidth =70;
[clusterInfo,pointToClusterMap] = MeanShiftClustering(ptRawData, bandwidth, ... 
                                                      'flagDebug', false, ...
                                                      'kernel', 'gaussian', ...
                                                      'flagUseKDTree', true);
% plot result

h2 = figure(2);hold off;
imagesc(currentImg,[min_int max_int]); axis off; axis image;
title(['Detection on Frame: ', num2str(iFrame)]);
    
colormap(gray);
hold on;

sigma_x=[];
sigma_y=[];
center_x=[];
center_y=[];
clustering_nCell = 0;

for iCell = 1:numel(clusterInfo)
    ptCurClusterCenter = clusterInfo(iCell).ptClusterCenter;
    plot( ptRawData(pointToClusterMap==iCell, 1), ...
        ptRawData(pointToClusterMap==iCell, 2), ...
        '.', 'Color', color_array(iCell,:));
    
    plot(ptCurClusterCenter(1),ptCurClusterCenter(2),'.','MarkerFaceColor',color_array(iCell,:), 'MarkerSize',10)
    
    X = ptRawData(pointToClusterMap==iCell, 1);
    Y = ptRawData(pointToClusterMap==iCell, 2);
    
%     [prmVect prmStd C res J] = fitGaussian2D([X; Y]);
    
    sigma_x(iCell) = std(X);
    sigma_y(iCell) = std(Y);
    center_x(iCell) = ptCurClusterCenter(1);
    center_y(iCell) = ptCurClusterCenter(2);
    cell_width = 2*round(sigma_x(iCell)*2.5)+1;
    cell_height = 2*round(sigma_y(iCell)*2.5)+1;
    
    position(1) = round(center_x(iCell)) - round(sigma_x(iCell)*2.5 );
    position(2) = round(center_y(iCell)) - round(sigma_y(iCell)*2.5);
    position(3) = cell_width;
    position(4) = cell_height;  
    if(cell_width>min_size && length(X)>min_size*min_size/3)
        clustering_nCell = clustering_nCell+1;
        position_array_from_clustering{iFrame}(clustering_nCell,1:4) = position;
        cluster_xy_s_from_clustering{iFrame,clustering_nCell} = [X Y];        
        cluster_sigma_s_from_clustering{iFrame,clustering_nCell} = ...
            [sigma_x(iCell);    sigma_y(iCell);    center_x(iCell);    center_y(iCell);    cell_width;    cell_height;];
        
        mean_intensity_from_clustering(iFrame,clustering_nCell) = ...
             sum(sum(double(Mask(...
             max(1,round(position(2))):...
             min(round(position(2)+position(4)),size(Mask,1)), ...
             max(1,round(position(1))):...
             min(round(position(1)+position(3)),size(Mask,2))...
             ))));
         current_position = position;

     
        X = [current_position(1);...
            current_position(1)+ current_position(3);...
            current_position(1)+current_position(3);...
            current_position(1);...
            current_position(1);];
        
        Y = [current_position(2);...
            current_position(2);...
            current_position(2)+ current_position(4);...
            current_position(2)+current_position(4);...
            current_position(2);];
        
        h2= figure(2);        
        plot(X,Y,'color',color_array(clustering_nCell,:));
        

    end
end


mean_intensity_matrix = zeros(nFrame, clustering_nCell);
VIF_intensity = zeros(nFrame, clustering_nCell);

% title(sprintf('Mean-shift clustering result - %d clusters were found', numel(clusterInfo)));
saveas(h2,[MD.outputDirectory_,filesep,'CC_tracking',filesep,'mean_shift_clustering_frame_',num2str(iFrame),'.tif']);


