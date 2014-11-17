function VIF_outgrowth_measurement(movieData)

% Created June 2012 by Liya Ding, Matlab R2011b

% Find the package of Filament Analysis
nPackage = length(movieData.packages_);

indexFilamentPackage = 0;
for i = 1 : nPackage
    if(strcmp(movieData.packages_{i}.getName,'FilamentAnalysis')==1)
        indexFilamentPackage = i;
        break;
    end
end

if(indexFilamentPackage==0)
    msgbox('Need to be in Filament Package for now.')
    return;
end


%%

nProcesses = length(movieData.processes_);

indexFilamentSegmentationProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Filament Segmentation')==1)
        indexFilamentSegmentationProcess = i;
        break;
    end
end

if indexFilamentSegmentationProcess==0
    msgbox('Please set parameters for Filament Segmentation and run.')
    return;
end

funParams=movieData.processes_{indexFilamentSegmentationProcess}.funParams_;

selected_channels = funParams.ChannelIndex;
Combine_Way = funParams.Combine_Way;
Cell_Mask_ind = funParams.Cell_Mask_ind;
VIF_Outgrowth_Flag = funParams.VIF_Outgrowth_Flag;

indexSteerabeleProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Steerable filtering')==1)
        indexSteerabeleProcess = i;
        break;
    end
end

if indexSteerabeleProcess==0 && Combine_Way~=2
    msgbox('Please run steerable filtering first.')
    return;
else
    funParams_st = movieData.processes_{indexSteerabeleProcess}.funParams_;
    ImageFlattenFlag = funParams_st.ImageFlattenFlag;
end

indexFlattenProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Image Flatten')==1)
        indexFlattenProcess = i;
        break;
    end
end

if indexFlattenProcess==0  && ImageFlattenFlag == 2
    display('The setting shows you want to use flattened image for steerable filtering. Please set parameters for Image Flatten and run.')
    return;
end

indexCellSegProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Mask Refinement')==1)
        indexCellSegProcess = i;
        break;
    end
end

if indexCellSegProcess == 0 && Cell_Mask_ind == 1
    msgbox('Please run segmentation and refinement first.')
    return;
end


nFrame = movieData.nFrames_;

% If the user set an cell ROI read in
if(exist([movieData.outputDirectory_,filesep,'MD_ROI.tif'],'file'))
    user_input_mask = imread([movieData.outputDirectory_,filesep,'MD_ROI.tif']);
end


%% Output Directories

FilamentSegmentationProcessOutputDir  = [movieData.packages_{indexFilamentPackage}.outputDirectory_, filesep 'FilamentSegmentation'];
if (~exist(FilamentSegmentationProcessOutputDir,'dir'))
    mkdir(FilamentSegmentationProcessOutputDir);
end

for iChannel = selected_channels
    FilamentSegmentationChannelOutputDir = [FilamentSegmentationProcessOutputDir,'/Channel',num2str(iChannel)];
    if (~exist(FilamentSegmentationChannelOutputDir,'dir'))
        mkdir(FilamentSegmentationChannelOutputDir);
    end
    
    movieData.processes_{indexFilamentSegmentationProcess}.setOutImagePath(iChannel,FilamentSegmentationChannelOutputDir);
end


%%


% Gets the screen size
scrsz = get(0,'ScreenSize');
save_fig_flag = 0;

for iChannel = selected_channels
    
    % Get frame number from the title of the image, this not neccesarily
    % the same as iFrame due to some shorting problem of the channel
    
       Channel_FilesNames = movieData.channels_(iChannel).getImageFileNames(1:movieData.nFrames_);
    
    filename_short_strs = uncommon_str_takeout(Channel_FilesNames);

    
    % Make output directory for the steerable filtered images
    FilamentSegmentationChannelOutputDir = [FilamentSegmentationProcessOutputDir,'/Channel',num2str(iChannel)];
    if (~exist(FilamentSegmentationChannelOutputDir,'dir'))
        mkdir(FilamentSegmentationChannelOutputDir);
    end
    
    SteerableChannelOutputDir = movieData.processes_{indexSteerabeleProcess}.outFilePaths_{iChannel};
    
    HeatOutputDir = [FilamentSegmentationChannelOutputDir,'/HeatOutput'];
    
    if (~exist(HeatOutputDir,'dir'))
        mkdir(HeatOutputDir);
    end
    
    HeatEnhOutputDir = [HeatOutputDir,'/Enh'];
    
    if (~exist(HeatEnhOutputDir,'dir'))
        mkdir(HeatEnhOutputDir);
    end
    
    HeatEnhBoundOutputDir = [HeatOutputDir,'/Enh_bound'];
    
    if (~exist(HeatEnhBoundOutputDir,'dir'))
        mkdir(HeatEnhBoundOutputDir);
    end
    
    DataOutputDir = [FilamentSegmentationChannelOutputDir,'/DataOutput'];
    
    if (~exist(DataOutputDir,'dir'))
        mkdir(DataOutputDir);
    end
    

    H_close = fspecial('disk',39);
    H_close = H_close>0;
    
    
    
    % this line in commandation for shortest version of filename
    filename_shortshort_strs = all_uncommon_str_takeout(Channel_FilesNames);
    
    
    for iFrame = 1 : nFrame
        disp(['Frame: ',num2str(iFrame)]);
        
        currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        currentImgFlatten = currentImg;
        % Read in the intensity image, flattened or original
        if indexFlattenProcess > 0
            try
            currentImgFlatten = imread([movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel}, filesep, 'flatten_',filename_short_strs{iFrame},'.tif']);
            catch
               try
                   % this is for the old version
                    currentImgFlatten = imread([movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel}, filesep, 'flatten_',filename_shortshort_strs{iFrame},'.tif']);           
               end
            end
        end
        
        currentImg = currentImgFlatten;
        
        
        try
            load([SteerableChannelOutputDir, filesep, 'steerable_',...
                filename_short_strs{iFrame},'.mat']);            
        catch
            % in the case of only having the short-old version
            load([SteerableChannelOutputDir, filesep, 'steerable_',...
                filename_shortshort_strs{iFrame},'.mat']);            
        end
        
        
        
        if funParams.Cell_Mask_ind == 1
            MaskCell = movieData.processes_{indexCellSegProcess}.loadChannelOutput(iChannel,iFrame);
        else
            if funParams.Cell_Mask_ind == 2
                MaskCell = user_input_mask>0;
            else
                MaskCell = ones(size(currentImg,1),size(currentImg,2));
            end
        end
        
        try
            load([DataOutputDir,'/steerable_vote_',...
                filename_short_strs{iFrame},'.mat'],...
                'orienation_map_filtered', ...
                'MAX_st_res', 'current_seg','SteerabelRes_Segment');
        catch
            % in the case of only having the short-old version
            load([DataOutputDir,'/steerable_vote_',...
                filename_shortshort_strs{iFrame},'.mat'],...
                'orienation_map_filtered', ...
                'MAX_st_res', 'current_seg','SteerabelRes_Segment');
            
        end        
        
        if iFrame==1
            % Read in the initial circle from the 'start_ROI.tif' file
            MaskFirstFrame = imread([movieData.outputDirectory_,'/start_ROI.tif']);
            MaskFirstFrame = (MaskFirstFrame)>0;
            
            RoiYX = bwboundaries(MaskFirstFrame);
            RoiYX = RoiYX{1};
        end
        
        current_seg = SteerabelRes_Segment.*MaskCell;
        %current_seg = SteerabelRes_Segment;
        
        incircle_seg = SteerabelRes_Segment.*MaskFirstFrame;
        
        incircle_seg_dilate = imdilate(incircle_seg, H_close);
        
        all_seg_dilate = imdilate(current_seg, H_close);
        
        ind_incircle = find(incircle_seg_dilate>0);
        labelMask = bwlabel(all_seg_dilate);
        dilate_keep = labelMask.*0;
        
        ind_new_labels = unique(labelMask(ind_incircle));
        
        for li = 1 :  length(ind_new_labels)
            dilate_keep(find(labelMask == ind_new_labels(li))) = ind_new_labels(li);
        end
        
        %         labelMask(ind_incircle) = keep_largest_area(all_seg_dilate);
        
        current_seg = dilate_keep.*current_seg;
        
        labelMask = bwlabel(current_seg);
        
        ob_prop = regionprops(labelMask,'Area','MajorAxisLength','Eccentricity','Centroid');
        
        if length(ob_prop) > 1
            obAreas = [ob_prop.Area];
            obLongaxis = [ob_prop.MajorAxisLength];
            obEccentricity = [ob_prop.Eccentricity];
            obCentroid = [ob_prop.Centroid];
            
            for i_area = 1 : length(obAreas)
                centroid_x = round(obCentroid(2*i_area-1));
                centroid_y = round(obCentroid(2*i_area));
                
                if( MaskFirstFrame(centroid_y, centroid_x) == 0)
                    if obAreas(i_area) < 12 || obLongaxis(i_area) < 10
                        labelMask(labelMask==i_area) = 0;
                    end
                else
                    if obAreas(i_area) < 7 || obLongaxis(i_area) < 5
                        labelMask(labelMask==i_area) = 0;
                    end
                end
            end
        end
        
        current_seg = labelMask > 0;
        
        if iFrame==1
            % Get the first segmented results for the base of comparison
            current_seg_firstframe = current_seg;
            
            current_seg_inside_firstframe = current_seg_firstframe.*(double(MaskFirstFrame));
            seg_sum_inside_firstframe = sum(sum(current_seg_inside_firstframe));
        end
        
        
        current_seg_outside = current_seg.*(1-MaskFirstFrame);
        
        Hue = (-orienation_map_filtered(:)+pi/2)/(pi)-0.2;
        Hue(find(Hue>=1)) = Hue(find(Hue>=1)) -1;
        Hue(find(Hue<0)) = Hue(find(Hue<0)) +1;
        
        Sat = Hue*0+1;
        Value = Hue*0+1;
        RGB_seg_orient_heat_array = hsv2rgb([Hue Sat Value]);
        R_seg_orient_heat_map = col2im(RGB_seg_orient_heat_array(:,1),[1 1],[size(current_seg,1) size(current_seg,2)]);
        G_seg_orient_heat_map = col2im(RGB_seg_orient_heat_array(:,2),[1 1],[size(current_seg,1) size(current_seg,2)]);
        B_seg_orient_heat_map = col2im(RGB_seg_orient_heat_array(:,3),[1 1],[size(current_seg,1) size(current_seg,2)]);
        
        enhanced_im_r = currentImg;
        enhanced_im_g = currentImg;
        enhanced_im_b = currentImg;
        
        enhanced_im_r(find(current_seg>0)) = 255*R_seg_orient_heat_map(find(current_seg>0));
        enhanced_im_g(find(current_seg>0)) = 255*G_seg_orient_heat_map(find(current_seg>0));
        enhanced_im_b(find(current_seg>0)) = 255*B_seg_orient_heat_map(find(current_seg>0));
        
        RGB_seg_orient_heat_map(:,:,1) = enhanced_im_r;
        RGB_seg_orient_heat_map(:,:,2) = enhanced_im_g;
        RGB_seg_orient_heat_map(:,:,3) = enhanced_im_b;
        
        %         RGB_seg_orient_heat_map = uint8(RGB_seg_orient_heat_map*255);
        
        h12 = figure(12);
        hold off;
        
        imagesc(RGB_seg_orient_heat_map);
        axis image; axis off;
        title(['Frame ',num2str(iFrame),', percentage of out growth:', ...
            num2str(100*sum(sum(current_seg_outside))/seg_sum_inside_firstframe), '%'],'FontSize',15);
        
        saveas(h12,[HeatEnhOutputDir,'/Enh_VIF_heat_display_',...
            filename_short_strs{iFrame},'.tif']);
        if(save_fig_flag==1)
            saveas(h12,[HeatEnhOutputDir,'/Enh_VIF_heat_display_',...
                filename_short_strs{iFrame},'.fig']);
        end
        
        hold on; plot(RoiYX(:,2),RoiYX(:,1),'m');
        
        saveas(h12,[HeatEnhOutputDir,'_bound/Enh_Bound_VIF_heat_display_',...
            filename_short_strs{iFrame},'.tif']);
        if(save_fig_flag==1)
            saveas(h12,[HeatEnhOutputDir,'_bound/Enh_Bound_VIF_heat_display_',...
                filename_short_strs{iFrame},'.fig']);
        end
        seg_outside_current(iChannel, iFrame) = sum(sum(current_seg_outside));
        ratio_outside_firstframeinside(iChannel, iFrame) = sum(sum(current_seg_outside))/seg_sum_inside_firstframe;
        
    end
    display('The outgrowth percentage results:');
    display(ratio_outside_firstframeinside'*100);
    
    % Save the outgrowth results
    save([FilamentSegmentationChannelOutputDir,'/seg_outside.mat'],'seg_outside_current','seg_sum_inside_firstframe','ratio_outside_firstframeinside');
end


for iChannel = selected_channels
    
    %% Simply use intensity as outgrowth calculation
    
        org_int_mean_inside_array=zeros(1,nFrame);
        org_int_sum_inside_array=zeros(1,nFrame);
        flatten_int_mean_inside_array=zeros(1,nFrame);
        flatten_int_sum_inside_array=zeros(1,nFrame);
        
        org_int_mean_outside_array=zeros(1,nFrame);
        org_int_sum_outside_array=zeros(1,nFrame);
        flatten_int_mean_outside_array=zeros(1,nFrame);
        flatten_int_sum_outside_array=zeros(1,nFrame);

    
    for iFrame = 1 : nFrame
        disp(['Frame: ',num2str(iFrame)]);
      
        % Read in the intensity image, flattened or original        
        currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        currentImgFlatten = currentImg;
        if indexFlattenProcess > 0
            currentImgFlatten = imread([movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel}, filesep, 'flatten_',filename_short_strs{iFrame},'.tif']);
        end
        
       
        
        currentImg_inside = currentImg(find(MaskFirstFrame>0));
        currentImg_outside = currentImg(find(MaskFirstFrame==0));
        currentImgFlatten_inside = currentImgFlatten(find(MaskFirstFrame>0));
        currentImgFlatten_outside = currentImgFlatten(find(MaskFirstFrame==0));
        
        if iFrame==1
            % Get the first segmented results for the base of comparison
            current_img_firstframe = currentImg;
            int_sum_inside_firstframe = sum(sum(current_img_firstframe));
           org_int_mean_inside_firstframe = sum(currentImg_inside);
           org_int_sum_inside_firstframe = mean(currentImg_inside);
           flatten_int_mean_inside_firstframe = sum(currentImgFlatten_inside);
           flatten_int_sum_inside_firstframe = mean(currentImgFlatten_inside);
        end
        
        org_int_mean_inside_array(1,iFrame) = sum(currentImg_inside);
        org_int_sum_inside_array(1,iFrame)= mean(currentImg_inside);
        flatten_int_mean_inside_array(1,iFrame) = sum(currentImgFlatten_inside);
        flatten_int_sum_inside_array(1,iFrame)= mean(currentImgFlatten_inside);
        
        org_int_mean_outside_array(1,iFrame) = sum(currentImg_outside);
        org_int_sum_outside_array(1,iFrame)= mean(currentImg_outside);
        flatten_int_mean_outside_array(1,iFrame) = sum(currentImgFlatten_outside);
        flatten_int_sum_outside_array(1,iFrame)= mean(currentImgFlatten_outside);
      
        level1 = thresholdOtsu(currentImg);
        currentImg = double(currentImg>level1).*double(currentImg);
        
        current_img_outside = double(currentImg).*(double(1-MaskFirstFrame));
        
        h12 = figure(12);
        hold off;
        
        imagesc(currentImg); colormap(gray);
        axis image; axis off;
        title(['Frame ',num2str(iFrame),', percentage of out growth based on intensity: ', ...
            num2str(100*sum(sum(current_img_outside))/int_sum_inside_firstframe), '%'],'FontSize',15);
        
        hold on; plot(RoiYX(:,2),RoiYX(:,1),'m');
        
        saveas(h12,[HeatEnhOutputDir,'_bound/Int_display_',num2str(iFrame),'.tif']);
        
        int_outside_current_intensity(iChannel, iFrame) = sum(sum(current_img_outside));
        ratio_int_outside_firstframeinside(iChannel, iFrame) = sum(sum(current_img_outside))/int_sum_inside_firstframe;
        
    end
    
    display('The outgrowth percentage results based on intensity:');
    display(ratio_int_outside_firstframeinside'*100);
    
    % Save the outgrowth results
    save([FilamentSegmentationChannelOutputDir,'/seg_outside.mat'], ...
        'seg_outside_current','seg_sum_inside_firstframe','ratio_outside_firstframeinside',...
        'ratio_int_outside_firstframeinside');
    
   
    xlswrite([FilamentSegmentationProcessOutputDir,'/channel_',num2str(iChannel),'_seg_growth.xls'], ...
        {'Frame';'First Frame Filament Amount in ROI';...
        'Each Frame Filament Amount outside ROI';...
        'Outgrowth Ratio(with Filament Segmentation)';...
        'Outgrowth Ratio(with Intensity Segmentation)'; ...
        'Intensity Sum in ROI';'Intensity Sum out ROI';...
        'Intensity Mean in ROI';'Intensity Mean out ROI';...
        'Intensity Sum in ROI(Flattened Image)';'Intensity Sum out ROI(Flattened Image)';...
        'Intensity Mean in ROI(Flattened Image)';'Intensity Mean out ROI(Flattened Image)';}, 1,'A1');
    
    xlswrite([FilamentSegmentationProcessOutputDir,'/channel_',num2str(iChannel),'_seg_growth.xls'], ...
       (1:nFrame), 1,'B1');
    
   xlswrite([FilamentSegmentationProcessOutputDir,'/channel_',num2str(iChannel),'_seg_growth.xls'], ...
       [repmat(seg_sum_inside_firstframe(iChannel),1,nFrame); ...
       seg_outside_current(iChannel,:); ...
       ratio_outside_firstframeinside(iChannel,:);...
       ratio_int_outside_firstframeinside(iChannel,:);...
       org_int_mean_inside_array;...
       org_int_mean_outside_array;...
       org_int_sum_inside_array; ...
       org_int_sum_outside_array;...
       flatten_int_mean_inside_array; ...
       flatten_int_mean_outside_array;...
       flatten_int_sum_inside_array;...
       flatten_int_sum_outside_array], 1,'B2');
   
        
    % Display the curve
    h2 = figure(2);hold off;    
%     [AX,H1,H2] =  plotyy(1:nFrame,[repmat(seg_sum_inside_firstframe(iChannel),1,nFrame); ...
%         seg_outside_current(iChannel,:)]',1:nFrame,[ratio_outside_firstframeinside(iChannel,:);...
%         ratio_int_outside_firstframeinside(iChannel,:)]');
%     legend('Base Filament to compare','Outgrowth Filament','Ratio of Outgrowing Filament','Ratio of Outgrow if using intensity');
     
    plot([ratio_outside_firstframeinside(iChannel,:);...
        ratio_int_outside_firstframeinside(iChannel,:)]'*100);hold on;
    plot([ratio_outside_firstframeinside(iChannel,:);...
        ratio_int_outside_firstframeinside(iChannel,:)]'*100,'.');

    ylabel('Percentage(%)','FontSize',15);
    set(gca,'FontSize',15);
    if(isempty(movieData.timeInterval_))
      xlabel('Frame','FontSize',15);
      set(gca,'xtick',1:1:nFrame);  
    else
      ylabel('Frame','FontSize',15);
      X_label = cell(1,nFrame);
      for iFrame = 1 :nFrame
          iTime = iFrame*movieData.timeInterval_;
        X_label{iFrame}= [num2str(round(iTime/60),'%02d'),':',num2str(mod(iTime,60),'%02d')];
      end      
      set(gca,'xticklabel',X_label,'FontSize',15);  
      xlabel('mm:ss','FontSize',15);
    end
    title(['Channel ',num2str(iChannel),' VIF growth results'],'FontSize',15);
     legend('Based on Filament Segmenation','Based on Intensity');
     
    saveas(h2,[FilamentSegmentationProcessOutputDir,'/channel_',num2str(iChannel),'_results_plot.tif']);
    
end