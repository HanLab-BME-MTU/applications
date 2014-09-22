function [position_array, present_cells, mean_intensity_matrix] = CC_tracking_autoinitialization(MD, iChannel, bandwidth, min_size)
% Cross Correlation based tracking with automatic initialization
% Input: MD the movieData object
%        bandwidth: for mean shift clustering, set as the mean distance for the cells in pixels        
% Output: None, save images to the CC_tracking folder in the same directory as the MD object
% Liya 
% Oct 23, 2012

% % To get the index for different processes
% package_process_ind_script;

nFrame = MD.nFrames_;

color_array= [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1;0 1 1 ;rand(50,3)];
iFrame = 1;
if(~exist([MD.outputDirectory_,'/CC_tracking/'],'dir'))
    mkdir([MD.outputDirectory_,'/CC_tracking/']);
end

intensity_pool = [];
min_int = 2^16;
max_int = 0;

for iFrame = 1 : nFrame;       
    currentImg = double(MD.channels_(iChannel).loadImage(iFrame));
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

for bin_ind = 1:1:100
    if(sum(hist_int(1:bin_ind))/sum(hist_int(1:end))>1/100)
        min_int = bin_int(bin_ind) - (bin_int(2)-bin_int(1))*2;
        break;
    end
end

iFrame=1;
currentImg = double(MD.channels_(iChannel).loadImage(iFrame));

std_muler=1;
level1 = thresholdOtsu(currentImg(currentImg>mean2(currentImg)+std_muler*std2(currentImg)));

% threshold out the bright nucleas part
Mask  = currentImg>level1;

[inda,indb] = find(Mask>0);

ptRawData = [(indb) (inda)];

% Run mean-shift
% bandwidth =70;
[clusterInfo,pointToClusterMap] = MeanShiftClustering(ptRawData, bandwidth, ... 
                                                      'flagDebug', false, ...
                                                      'kernel', 'gaussian', ...
                                                      'flagUseKDTree', true);
% plot result
currentImg = double(MD.channels_(iChannel).loadImage(iFrame));

h1 = figure(1);hold off;
imagesc(currentImg,[min_int max_int]); axis off; axis image;
colormap(gray);
hold on;

sigma_x=[];
sigma_y=[];
center_x=[];
center_y=[];
nCell = 0;

mean_intensity_start = [];

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

mean_intensity_matrix = zeros(nFrame, nCell);

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
present_cells{1}=ones(1,nCell);


for iFrame = 1 : nFrame;   
    
    previoucurrentImg = currentImg;
    currentImg = double(MD.channels_(iChannel).loadImage(iFrame));
    
    current_level1 = thresholdOtsu(currentImg(currentImg>mean2(currentImg)+std_muler*std2(currentImg)));

    current_Mask  = currentImg>level1;


    h3 = figure(3);hold off;    
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
        plot(X,Y,'color',color_array(iCell,:));

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

         if(mean_intensity_matrix(iFrame, iCell) < mean_intensity_start(iCell)/20)
               present_cells{iFrame}(iCell)=0;
         end
        
%     These are for level set segmentation
%     Not ready yet
        
%         mask_cells( round(current_position(2)):round(current_position(2))+round(current_position(4)),...
%             round(current_position(1)):round(current_position(1))+round(current_position(3)),iCell)=1;
                
    end
    title(['Tracking Frame ',num2str(iFrame)]);
    print(h3,'-dtiff',[MD.outputDirectory_,'/CC_tracking/tracking_',num2str(iFrame),'.tif']);

%     These are for level set segmentation
%     Not ready yet

%     I = double(currentImg);
%     I = (I-min(min(I)))/(max(max(I))-min(min(I)));
    
%     chenvese(currentImg, mask_cells,1000,0.05,'multiphase');
    
end

save([MD.outputDirectory_,'/CC_tracking/cc_tracking_results.mat'],'position_array','present_cells','mean_intensity_matrix');

currentImg = double(MD.channels_(1).loadImage(1));

h3 = figure(3);hold off; imagescc(currentImg);hold on;
   


for iFrame = 1 : nFrame    
    ind_cell = find(present_cells{iFrame}>0);
    
    for iCell = ind_cell     
        current_position = position_array{iFrame}(iCell,1:4);
        tracked_pos = current_position(1:2)+(current_position(3:4)+1)/2;
        
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
        plot(X,Y,'color',color_array(iCell,:));
        if iFrame >1 
            previous_position = position_array{iFrame-1}(iCell,1:4);
            tracked_old_pos = previous_position(1:2)+(previous_position(3:4)+1)/2;
            quiver(tracked_old_pos(1), tracked_old_pos(2), tracked_pos(1)-tracked_old_pos(1), tracked_pos(2)-tracked_old_pos(2),'color',color_array(iCell,:));
        
        end
    end
    

end

saveas(h3,[MD.outputDirectory_,'/CC_tracking/whole_tracking.tif']);
