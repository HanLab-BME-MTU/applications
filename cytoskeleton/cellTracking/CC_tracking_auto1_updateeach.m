function [position_array, present_cells, mean_intensity_matrix,VIF_intensity] = ...
    CC_tracking_auto1_updateeach(MD, Tracking_iChannel, VIF_iChannel, bandwidth, min_size)
% Cross Correlation based tracking with automatic initialization
% Input: MD the movieData object
%        bandwidth: for mean shift clustering, set as the mean distance for the cells in pixels
% Output: None, save images to the CC_tracking folder in the same directory as the MD object

% Update after each frame for new cell and cell shape changes
% Liya
% Nov 29, 2012

% % To get the index for different processes
% package_process_ind_script;

nFrame = MD.nFrames_;

color_array= [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1;0 1 1 ;rand(150,3)];
iFrame = 1;
if(~exist([MD.outputDirectory_,'/CC_tracking/'],'dir'))
    mkdir([MD.outputDirectory_,'/CC_tracking/']);
end

intensity_pool = [];
min_int = 2^16;
max_int = 0;

for iFrame = 1 : nFrame;
    currentImg = double(MD.channels_(Tracking_iChannel).loadImage(iFrame));
    currentImg = currentImg(1:10:end,1:10:end);
    intensity_pool = [intensity_pool;currentImg(:)];
end
[hist_int, bin_int] =  hist(intensity_pool,100);


for bin_ind = 100:-1:1
    if(sum(hist_int(bin_ind:end))/sum(hist_int(1:end))>1/100)
        max_int = max(0, bin_int(bin_ind)+(bin_int(2)-bin_int(1))*2);
        break;
    end
end

for bin_ind = 100:-1:1
    if(sum(hist_int(bin_ind:end))/sum(hist_int(1:end))>1/100)
        max_99_percent_int = max(0, bin_int(bin_ind)+(bin_int(2)-bin_int(1))*2);
        break;
    end
end


for bin_ind = 100:-1:1
    if(sum(hist_int(bin_ind:end))/sum(hist_int(1:end))>2/100)
        max_98_percent_int = max(0, bin_int(bin_ind)+(bin_int(2)-bin_int(1))*2);
        break;
    end
end


for bin_ind = 100:-1:1
    if(sum(hist_int(bin_ind:end))/sum(hist_int(1:end))>3/100)
        max_97_percent_int = max(0, bin_int(bin_ind)+(bin_int(2)-bin_int(1))*2);
        break;
    end
end


for bin_ind = 100:-1:1
    if(sum(hist_int(bin_ind:end))/sum(hist_int(1:end))>20/100)
        max_80_percent_int = max(0, bin_int(bin_ind)+(bin_int(2)-bin_int(1))*2);
        break;
    end
end

for bin_ind = 100:-1:1
    if(sum(hist_int(bin_ind:end))/sum(hist_int(1:end))>60/100)
        max_60_percent_int = max(0, bin_int(bin_ind)+(bin_int(2)-bin_int(1))*2);
        break;
    end
end

for bin_ind = 1:1:100
    if(sum(hist_int(1:bin_ind))/sum(hist_int(1:end))>1/100)
        min_int = bin_int(bin_ind) - (bin_int(2)-bin_int(1))*2;
        break;
    end
end



iFrame=1;
currentImg = double(MD.channels_(Tracking_iChannel).loadImage(iFrame));

for_first_threshold_img = currentImg;
for_first_threshold_img(find(for_first_threshold_img>max_int*1.5)) = max_int*1.5;

std_muler=0.3;
level1 = thresholdOtsu(for_first_threshold_img(for_first_threshold_img>mean2(for_first_threshold_img)+std_muler*std2(for_first_threshold_img)));

% threshold out the bright nucleas part
Mask  = currentImg>max_99_percent_int;

[inda,indb] = find(Mask>0);

ptRawData = [(indb) (inda)];

% Run mean-shift
% bandwidth =70;
[clusterInfo,pointToClusterMap] = MeanShiftClustering(ptRawData, bandwidth, ...
    'flagDebug', false, ...
    'kernel', 'gaussian', ...
    'flagUseKDTree', true);
% plot result
currentImg = double(MD.channels_(Tracking_iChannel).loadImage(iFrame));

h1 = figure(1);hold off;
imagesc(currentImg,[min_int max_int]); axis off; axis image;
colormap(gray);
hold on;

sigma_x=[];
sigma_y=[];
center_x=[];
center_y=[];
nCell = 0;



for iCell = 1:numel(clusterInfo)
    ptCurClusterCenter = clusterInfo(iCell).ptClusterCenter;
X_clustered = ptRawData(pointToClusterMap==iCell, 1);
if(isempty(X_clustered))
stophere=1;
end
if(size(ptRawData,2)==1)
stophere=1;
end
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
        nCell = nCell+1;
        position_array{iFrame}(nCell,1:4) = position;
        mean_intensity_start(nCell) = ...
            sum(sum(double(Mask(...
            max(1,round(position(2))):...
            min(round(position(2)+position(4)),size(Mask,1)), ...
            max(1,round(position(1))):...
            min(round(position(1)+position(3)),size(Mask,2))...
            ))));
        
    end
end

mean_intensity_matrix = zeros(nFrame, nCell*2);
VIF_intensity = zeros(nFrame, nCell*2);

title( sprintf('Mean-shift clustering result - %d clusters were found', numel(clusterInfo)));
saveas(h1,[MD.outputDirectory_,'/CC_tracking/mean_shift_clustering.tif']);

% mask_cells = zeros(size(currentImg,1),size(currentImg,2),nCell);

dmax = [10 10];
subpix = 'none';
d0=[0 0];
img_width = size(currentImg,2);
img_height = size(currentImg,1);
pad_xy = [img_width/2 img_height/2];

present_cells = cell(1,nFrame);
present_cells{1} = ones(1,nCell);

position_array_from_clustering = cell(1,nFrame);
cluster_xy_s_from_clustering= cell(nFrame,nCell);
cluster_sigma_s_from_clustering= cell(nFrame,nCell);
mean_intensity_from_clustering= zeros(nFrame,nCell);

for iFrame = 1 : nFrame;
    
    previoucurrentImg = currentImg;
    currentImg = double(MD.channels_(Tracking_iChannel).loadImage(iFrame));
    
    
    currentVIFImg = double(MD.channels_(VIF_iChannel).loadImage(iFrame));
    
    current_level1 = thresholdOtsu(currentImg(currentImg>mean2(currentImg)+std_muler*std2(currentImg)));
    
    current_Mask  = currentImg>level1;
    
    
    h3 = figure(3);hold off;
    imagesc(currentImg,[min_int max_int]); axis off; axis image;
    colormap(gray);
    hold on;
    
    h13 = figure(13);hold off;
    imagesc(currentImg,[min_int max_int]); axis off; axis image;
    colormap(gray);
    hold on;
    
    if iFrame >1
        present_cells{iFrame} =present_cells{iFrame-1};
    end
    
    
    ind_cell = find(present_cells{iFrame}>0);
    
    for iCell = ind_cell
        if iFrame >1
            
            previous_position = position_array{iFrame-1}(iCell,1:4);
            
            tPos = previous_position(1:2)+(previous_position(3:4)+1)/2;
            tDim = previous_position(3:4);
            
            [displacement, CC_map] = ccbased_track(pad_boundary(previoucurrentImg),...
                [tPos(1) tPos(2)]+pad_xy,[tDim(1) tDim(2)],pad_boundary(currentImg),dmax,subpix,d0);
            
            current_position(3:4) = previous_position(3:4);
            current_position(1:2) = previous_position(1:2)+displacement;
            
            if current_position(1)<=0
                current_position(3) = 2*(round((current_position(3) + current_position(1) -2)/2)) + 1;
                current_position(1) = 1;
            end
            
            if current_position(2)<=0
                current_position(4) = 2*(round((current_position(4) + current_position(2) -2)/2)) + 1;
                current_position(2) = 1;
            end
            
            if current_position(1)+ current_position(3) > size(currentImg,2)-1
                current_position(3) = 2*(round((size(currentImg,2)-2 - current_position(1))/2)) + 1;
            end
            
            if current_position(2)+ current_position(4) > size(currentImg,1)-1
                current_position(4) = 2*(round((size(currentImg,1)-2 - current_position(2))/2)) + 1;
            end
            
            tracked_pos = current_position(1:2)+(current_position(3:4)+1)/2;
            
            position_array{iFrame}(iCell,1:4) = current_position;
            
        else
            current_position = position_array{iFrame}(iCell,1:4);
            tracked_pos = current_position(1:2)+(current_position(3:4)+1)/2;
        end
        %         h3 = figure(3);
        %         plot(tracked_pos(1),tracked_pos(2),'.','color',color_array(iCell,:));
        %         h13 = figure(13);
        %         plot(tracked_pos(1),tracked_pos(2),'.','color',color_array(iCell,:));
        %
        %         X = [current_position(1);...
        %             current_position(1)+ current_position(3);...
        %             current_position(1)+current_position(3);...
        %             current_position(1);...
        %             current_position(1);];
        %
        %         Y = [current_position(2);...
        %             current_position(2);...
        %             current_position(2)+ current_position(4);...
        %             current_position(2)+current_position(4);...
        %             current_position(2);];
        %
        %         h13 = figure(13);
        %         plot(X,Y,'color',color_array(iCell,:));
        %         h3 = figure(3);
        %         plot(X,Y,'color',color_array(iCell,:));
        
        if(current_position(3)<5||current_position(4)<5)
            present_cells{iFrame}(iCell)=0;
        end
        
        mean_intensity_matrix(iFrame, iCell)...
            = sum(sum(double(current_Mask(...
            max(1,round(current_position(2))):...
            min(round(current_position(2)+current_position(4)),size(currentImg,1)), ...
            max(1,round(current_position(1))):...
            min(round(current_position(1)+current_position(3)),size(currentImg,2))...
            ))));
        
        
        VIF_intensity(iFrame, iCell)...
            = mean2(double(currentVIFImg(...
            max(1,round(current_position(2))):...
            min(round(current_position(2)+current_position(4)),size(currentImg,1)), ...
            max(1,round(current_position(1))):...
            min(round(current_position(1)+current_position(3)),size(currentImg,2))...
            )));
        
        if(mean_intensity_matrix(iFrame, iCell) < mean_intensity_start(iCell)/20)
            present_cells{iFrame}(iCell)=0;
        end
        
        %     These are for level set segmentation
        %     Not ready yet
        
        %         mask_cells( round(current_position(2)):round(current_position(2))+round(current_position(4)),...
        %             round(current_position(1)):round(current_position(1))+round(current_position(3)),iCell)=1;
        
    end
    %     title(['Tracking Frame ',num2str(iFrame)]);
    %     print(h3,'-dtiff',[MD.outputDirectory_,'/CC_tracking/tracking_',num2str(iFrame),'.tif']);
    %    print(h13,'-dtiff',[MD.outputDirectory_,'/CC_tracking/VIF_',num2str(iFrame),'.tif']);
    
    
    
    
    %     These are for level set segmentation
    %     Not ready yet
    
    %     I = double(currentImg);
    %     I = (I-min(min(I)))/(max(max(I))-min(min(I)));
    
    %     chenvese(currentImg, mask_cells,1000,0.05,'multiphase');
    
    mean_shift_for_cc_tracking_each_frame;
    
    tracked_cell_taken = zeros(1, nCell);
    clustered_cell_taken = zeros(1,clustering_nCell);
    
    for iCell_tracking = 1 : nCell
        if present_cells{iFrame}(iCell_tracking)>0
            
            for iCell_clustering = 1 : clustering_nCell
                
                tracking_position = position_array{iFrame}(iCell_tracking,1:4);
                clustering_position = position_array_from_clustering{iFrame}(iCell_clustering,1:4);
                
                center_clustering = clustering_position(1:2)+(clustering_position(3:4)+1)/2;
                center_tracking = tracking_position(1:2)+(tracking_position(3:4)+1)/2;
                
                
                if  clustered_cell_taken(iCell_clustering)==0 && center_clustering(1)> tracking_position(1) && center_clustering(1)< tracking_position(1) +tracking_position(3) ...
                        && center_clustering(2)> tracking_position(2) && center_clustering(2)< tracking_position(2) +tracking_position(4)
                    
                    for_display_XY{iFrame,iCell_tracking} = cluster_xy_s_from_clustering{iFrame,iCell_clustering};
                    center_tracking = tracking_position(1:2)+(tracking_position(3:4)+1)/2;
                    size_new = cluster_sigma_s_from_clustering{iFrame,iCell_clustering}(1:2)
                    
                    cell_width_new = 2*round(size_new(1)*2.5)+1;
                    cell_height_new = 2*round(size_new(2)*2.5)+1;
                    
                    position_updated(1) = round(center_tracking(1) - round(size_new(1)*2.5) );
                    position_updated(2) = round(center_tracking(2)  - round(size_new(2)*2.5));
                    position_updated(3) = cell_width_new;
                    position_updated(4) = cell_height_new;
                    
                    position_array{iFrame}(iCell_tracking,1:4) = position_updated;
                    
                    tracked_cell_taken(iCell_tracking)=1;
                    clustered_cell_taken(iCell_clustering)=1;
                    break;
                end
                
            end
        end
    end
    
    
    for iCell_clustering = 1 : clustering_nCell
        
        clustering_position = position_array_from_clustering{iFrame}(iCell_clustering,1:4);
        center_clustering = clustering_position(1:2)+(clustering_position(3:4)+1)/2;
        size_new = cluster_sigma_s_from_clustering{iFrame,iCell_clustering}(1:2);
        
        if  clustered_cell_taken(iCell_clustering)==0 && mean(size_new) > bandwidth/2 && mean_intensity_from_clustering(iFrame,iCell_clustering)>  max_80_percent_int
            
            cell_width_new = 2*round(size_new(1)*2.5)+1;
            cell_height_new = 2*round(size_new(2)*2.5)+1;
            
            position_updated(1) = round(center_tracking(1) - round(size_new(1)*2.5) );
            position_updated(2) = round(center_tracking(2)  - round(size_new(2)*2.5));
            position_updated(3) = cell_width_new;
            position_updated(4) = cell_height_new;
            
            nCell = nCell + 1;
            position_array{iFrame}(nCell,1:4) = position_updated;
            for_display_XY{iFrame,nCell} = cluster_xy_s_from_clustering{iFrame,iCell_clustering};
            tracked_cell_taken(nCell)=1;
            present_cells{iFrame}(nCell)=1;
            mean_intensity_start(nCell) = mean_intensity_from_clustering(iFrame,iCell_clustering);
        end
        
    end
    
    for iCell_tracking = 1 : nCell
        if present_cells{iFrame}(iCell_tracking)>0
            
            clustering_position = position_array_from_clustering{iFrame}(iCell_clustering,1:4);
            center_clustering = clustering_position(1:2)+(clustering_position(3:4)+1)/2;
            size_new = cluster_sigma_s_from_clustering{iFrame,iCell_clustering}(1:2);
            
            if  tracked_cell_taken(iCell_tracking)==0
                display_img_cluster = currentImg*0;
                display_img_cluster(max(1,round(tracking_position(2))):...
             min(round(tracking_position(2)+tracking_position(4)),size(currentImg,1)), ...
             max(1,round(tracking_position(1))):...
             min(round(tracking_position(1)+tracking_position(3)),size(currentImg,2)))=1;
      
                display_img_cluster = display_img_cluster.*currentImg;
                
                [ind_y, ind_x] = find(display_img_cluster>max_98_percent_int);
                for_display_XY{iFrame,iCell_tracking} = [ind_x ind_y];
            end
        end
    end
    
    
       ind_cell = find(present_cells{iFrame}>0);
    currentImg = double(MD.channels_(Tracking_iChannel).loadImage(iFrame));
    
    for iCell = ind_cell
        current_position = position_array{iFrame}(iCell,1:4);
        tracked_pos = current_position(1:2)+(current_position(3:4)+1)/2;
        
        dots_XY = for_display_XY{iFrame,iCell};
        
        h3 = figure(3);
        plot(tracked_pos(1),tracked_pos(2),'.','color',color_array(iCell,:));
        h13 = figure(13);
        plot(tracked_pos(1),tracked_pos(2),'.','color',color_array(iCell,:));        
        
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
        
        h13 = figure(13);       title(['Tracking Frame ',num2str(iFrame)]);
        plot(X,Y,'color',color_array(iCell,:));
        h3 = figure(3);         title(['Tracking Frame ',num2str(iFrame)]);
        
        hold on; plot(dots_XY(:,1),dots_XY(:,2),'.','color',color_array(iCell,:));
        
     end
           print(h13,'-dtiff',[MD.outputDirectory_,'/CC_tracking/tracking_',num2str(iFrame),'.tif']);
        print(h3,'-dtiff',[MD.outputDirectory_,'/CC_tracking/clustered_tracking_',num2str(iFrame),'.tif']);

end



save([MD.outputDirectory_,'/CC_tracking/cc_tracking_results.mat'], ...
    'position_array','present_cells','mean_intensity_matrix','VIF_intensity');

currentImg = double(MD.channels_(1).loadImage(1));
currentVIFImg = double(MD.channels_(VIF_iChannel).loadImage(iFrame));

% h3 = figure(3);hold off; imagescc(currentImg);hold on;
% h13 = figure(13);hold off; imagescc(currentVIFImg);hold on;
%
for iFrame = 1 : 1: nFrame
    ind_cell = find(present_cells{iFrame}>0);
    currentImg = double(MD.channels_(Tracking_iChannel).loadImage(iFrame));
    
    for iCell = ind_cell
        current_position = position_array{iFrame}(iCell,1:4);
        tracked_pos = current_position(1:2)+(current_position(3:4)+1)/2;
        
        dots_XY = for_display_XY(iFrame,nCell);
        
        h3 = figure(3);
        plot(tracked_pos(1),tracked_pos(2),'.','color',color_array(iCell,:));
        h13 = figure(13);
        plot(tracked_pos(1),tracked_pos(2),'.','color',color_array(iCell,:));
        
        
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
        
        h13 = figure(13);       title(['Tracking Frame ',num2str(iFrame)]);
        plot(X,Y,'color',color_array(iCell,:));
        h3 = figure(3);         title(['Tracking Frame ',num2str(iFrame)]);
        plot(dots_XY(:,1),dots_XY(:,2),'.','color',color_array(iCell,:));
        
        print(h13,'-dtiff',[MD.outputDirectory_,'/CC_tracking/tracking_',num2str(iFrame),'.tif']);
        print(h3,'-dtiff',[MD.outputDirectory_,'/CC_tracking/clustered_tracking_',num2str(iFrame),'.tif']);
        
        
        %         if iFrame >1
        %             previous_position = position_array{iFrame-1}(iCell,1:4);
        %             tracked_old_pos = previous_position(1:2)+(previous_position(3:4)+1)/2;
        %             quiver(tracked_old_pos(1), tracked_old_pos(2), tracked_pos(1)-tracked_old_pos(1), tracked_pos(2)-tracked_old_pos(2),'color',color_array(iCell,:));
        %         end
    end
end
