function steerable_filter_segmentation(MD, varargin)
% steerable_filter_segmentation do multiple scale steerable filtering and local thresholding for VIF images    
%
% SYNOPSIS: steerable_filter_segmentation(MD)
%
% DESCRIPTION: The function gets the VIF image seuqence as in MD VIF channel 
%               and the mask obtained from basic thresholding and post-processing.
%               The fucntion does multiple scale steerable filtering for the
%               image and get the final response of each pixel as the maximum of all
%               responses and set the orientation according to the layer
%               that gives the maximum response, obtained in multiscaleSteerableDetector.
%
%               Local thresholding is applied both to the intensity image
%               and the steerable filtering response image. The final
%               segmentation is the union of the two segmentation.
%
%               The function creates 2 new directories for the segmentation
%               and steerable filtering results and saves the fig and mat
%               files in these directories.
%
% INPUT         MD :   movieData object has the VIF channel that has ran through process of
%                      basic thresholding and post-processing
%               optional:  displayflag: a flag to if show and save the figures. If 1, do.
%                          steerable_base: size of the first steerable
%                          filtering kernal
%                          combine_way: how the segmentation is obtained,
%                          intensity and steerable filtering, or only one
%                          of them: 'int_st_both','st_only','int_only'
%
% OUTPUT:  None, but save the results into mat file for each frame, including:
%                  currentImg:     intensity image
%                  orienation_map: orientation( without segmentation information)
%                  MAX_st_res:     steerable filtering response (maximum)
%                  current_seg:    segmentation of VIF
%                  Intensity_Segment: segmentation of VIF only based on intensity
%                  SteerableRes_Segment: segmentation of VIF only based on steerable filtering response

% Created Jan 2012 by Liya Ding, Matlab R2011b
% Modified according to new movieData structure March 2012, Liya

displayflag = 0;
steerable_base = 1.5;
combine_way = 'int_st_both';

if nargin > 1
    
    if ~(varargin{1} == 1 || varargin{1} == 0)
        error('The first argument should be 0 or 1');
    else
        displayflag = varargin{1};
    end
    
end

if nargin > 2
    
    if ~isfloat(varargin{2})
        error('The second argument should be a positive number.');
    else
        steerable_base = abs(varargin{2});
    end
    
end

if nargin > 3
    
    if ~ischar(varargin{3})
        error('The third argument should be a string.');
    else
        combine_way = varargin{3};
    end
    
end

nFrame = MD.nFrames_;

% Make two directories for the segmentation and steerable filtering results
% if there are not existing.
SegmentOutputDir = [MD.outputDirectory_,'/VIF_segmentation'];
SteerableOutputDir = [MD.outputDirectory_,'/VIF_SteerableDetector'];

if (~exist(SegmentOutputDir,'dir'))
    mkdir(SegmentOutputDir);
end

if (~exist(SteerableOutputDir,'dir'))
    mkdir(SteerableOutputDir);
end

% Find the process of segmentation mask refinement.
nProcesses = length(MD.processes_);

indexVIFSegProcess = 0;
for i = 1 : nProcesses
    if(strcmp(MD.processes_{i}.getName,'Mask Refinement')==1)
        indexVIFSegProcess = i;
        break;
    end
end

if indexVIFSegProcess==0
    msgbox('Please run segmentation and refinement first.')
    return;
end

% The process index for segmentation is the same for both channel

indexMTSegProcess = indexVIFSegProcess;

% % Define the channel index of VIF.
% indexVIFChannel = 0;
% nChannels = length(MD.channels_);
%
% for i = 1 : nChannels
%     if(strcmp(MD.channel_{i}.imageType_, 'VIF')==1)
%         indexVIFChannel = i;
%         break;
%     end
% end
% if indexVIFChannel ==0
%     msgbox('VIF channel needed to be defined.')
% end

% For now, it it not defined, so directly set as 2
indexVIFChannel = 2;

% And MT as the 1st channel
indexMTChannel = 1;

for iFrame = 1 : nFrame
    disp(['Frame: ',num2str(iFrame)]);
    
    if (~exist([SteerableOutputDir,'/steerable_segment_all_',num2str(nFrame),'.mat'],'file'))
        % Read in the intensity image.
        currentImg = MD.channels_(indexVIFChannel).loadImage(iFrame);
          currentImg = double(currentImg);
        
        currentImg = log(currentImg);
        
        currentImg = currentImg - min(min(currentImg)) +1;
        
        currentImg = sqrt(currentImg);
        
        currentImg = imfilter(currentImg, fspecial('gaussian',11),'replicate','same');
        
        % Get the thresholding results after post-processing
        Mask = MD.processes_{indexVIFSegProcess}.loadChannelOutput(indexVIFChannel,iFrame);
        MaskMT = MD.processes_{indexMTSegProcess}.loadChannelOutput(indexMTChannel,iFrame);
       
        
        Mask = Mask.*MaskMT;
        
        H_close = fspecial('disk',5);
        H_close = H_close>0;
        
%         MaskMT = imerode(MaskMT,H_close);
        Mask = imdilate(Mask, ones(5,5),'same');
        TightMask = Mask.*MaskMT;
                        
        % Make the mask bigger in order to include all
        MaskBig = imdilate(Mask, ones(25,25),'same');
        
        % Make another mask to estimate the background intensity
        MaskHuge = imdilate(Mask, ones(65,65),'same');
        MaskDiff = MaskHuge-Mask;
        
        % The intensity of out of the targetted cell region is set to
        % the mean intensity of the pixels in the ring region around the cell
        % boundary.
        background_gray = mean(currentImg(find(MaskDiff>0)));
        currentImg = currentImg.*MaskBig + background_gray*(1- MaskBig);
        
        % Steerable filtering using four scales one doubling the previous one.         
        % function multiscaleSteerableDetector will automatically merge the results
        [MAX_st_res, orienation_map, nms, scaleMap] = multiscaleSteerableDetector(currentImg, 4, steerable_base*[1 2 3]);
        
%         [MAX_st_resB, orienation_mapB, nmsB, scaleMapB] = multiscaleSteerableDetector(currentImg, 1, steerable_base*[1 2]);
       
        % Local thresholding for both intensity and steerabel filtering
        % response.        
        
        patch_size = round(steerable_base*3)*2+1;       
        pace_size = round(steerable_base*1)*2+1;
        
%         load([SteerableOutputDir,'/steerable_',num2str(iFrame),'.mat'],...
%             'Intensity_Segment','SteerabelRes_Segment');
%         current_seg = or(Intensity_Segment,SteerabelRes_Segment);
        
               
%         [level3, SteerabelRes_Segment_B ] = thresholdOtsu_local(MAX_st_resB,patch_size,pace_size,0);
 
        switch combine_way
            case 'int_st_both'
                
                 level0 = thresholdOtsu(MAX_st_res);
               thresh_Segment = MAX_st_res > level0;
                                
                [level1, SteerabelRes_Segment ] = thresholdOtsu_local(MAX_st_res,patch_size,pace_size,0);
                [level2, Intensity_Segment ] = thresholdOtsu_local(currentImg,patch_size,pace_size,0);
                
                %         [level3, NonMaxSup_Segment ] = thresholdOtsu_local(NonMaxSup_MAX_res,patch_size,pace_size,0);
                % The segmentation is set as the union of two segmentation.
                current_seg = or(Intensity_Segment,SteerabelRes_Segment);
                
            case 'st_only'                
%                 [level1, SteerabelRes_Segment ] = thresholdOtsu_local(MAX_st_res,patch_size,pace_size,0);
                
                [level1, SteerabelRes_Segment ] = thresholdOtsu_local_with_mask(MAX_st_res,MaskBig, patch_size,pace_size,0);
     
                
                % The segmentation is set as from steerable filtering
                current_seg = SteerabelRes_Segment; 
                Intensity_Segment = current_seg*0;
                
            case 'int_only'
                [level2, Intensity_Segment ] = thresholdOtsu_local(currentImg,patch_size,pace_size,0);
                % The segmentation is set as from intensity segmentation
                current_seg = Intensity_Segment; 
                SteerabelRes_Segment = current_seg*0;
            otherwise
                warning('Use the default of union');
                [level1, SteerabelRes_Segment ] = thresholdOtsu_local(MAX_st_res,patch_size,pace_size,0);
                [level2, Intensity_Segment ] = thresholdOtsu_local(currentImg,patch_size,pace_size,0);
                %         [level3, NonMaxSup_Segment ] = thresholdOtsu_local(NonMaxSup_MAX_res,patch_size,pace_size,0);
                % The segmentation is set as the union of two segmentation.
                current_seg = or(Intensity_Segment,SteerabelRes_Segment);
        end
        
%         SEG(:,:,1) = SteerabelRes_Segment_B;
%         level0 = thresholdOtsu(MAX_st_res);
%         thresh_Segment = MAX_st_res > level0;
%         [level1, SteerabelRes_Segment ] = thresholdOtsu_local(MAX_st_res,patch_size,pace_size,0);
%         SEG(:,:,2) = SteerabelRes_Segment;
%         SEG(:,:,3) = SteerabelRes_Segment*0;
%         figure;imagesc(SEG)
        

        current_seg = current_seg.*TightMask;        
        
        current_seg_big = imdilate(current_seg, ones(5,5));
        current_seg_big = imfill(current_seg_big,'holes');
        
        %Label all objects in the mask
        labelMask = bwlabel(current_seg_big);
        
        %Get their area
        obAreas = regionprops(labelMask,'Area');       %#ok<MRPBW>
        
        %First, check that there are objects to remove
        if length(obAreas) > 1
            obAreas = [obAreas.Area];
            %Sort by area
            [dummy,iSort] = sort(obAreas,'descend'); %#ok<ASGLU>
            %Keep only the largest requested number
            current_seg_big = labelMask == iSort(1);
        end
        
        current_seg_big = imerode(current_seg_big,ones(5,5));
        
        current_seg = current_seg_big.*current_seg;
        
         
        %%
%        
%         %Label all objects in the segmentation
%         labelMask = bwlabel(current_seg);
%        
%         %Get their area
%         ob_prop = regionprops(labelMask,'Area','MajorAxisLength','Orientation');       %#ok<MRPBW>
%         
%         %First, check that there are objects to remove
%         if length(ob_prop) > 1
%             obAreas = [ob_prop.Area];
%             obALongaxis = [ob_prop.MajorAxisLength];
%             obOrientation = [ob_prop.Orientation];
%             for i_o = 1 : length(ob_prop)
%                 figure(7); imagescc(labelMask == i_o);
%             end
%             
%         end      
%         
        
%         current_seg = current_seg.*SteerabelRes_Segment;
        
        %%
        % A smoothing done only at the steerable filtering results        
        orienation_map_filtered = OrientationSmooth(orienation_map, SteerabelRes_Segment);
        
        %%
        % Voting of the orientation field for the non-steerable filter
        % segmented places.
        tic
        OrientationVoted = OrientationVote(orienation_map,SteerabelRes_Segment,3,25);
        toc
        intensity_addon = current_seg - SteerabelRes_Segment ==1;
        if (~isempty(max(max(intensity_addon))>0))
            orienation_map_filtered(find(intensity_addon>0)) = OrientationVoted(find(intensity_addon>0));
        end
                 
        %%        
                
        if(displayflag==1)
            % If want to display and save images
            % figure 1, for image and the responses
            
            scrsz = get(0,'ScreenSize');
            h1 = figure(1);
            set(h1,'Position',scrsz);
            subplot(131);
            imagesc(currentImg);colormap(gray); axis image; axis off;
            title(['Frame: ',num2str(iFrame)],'Fontsize',18);
            
            subplot(132);imagesc(MAX_st_res);colormap(gray); axis image; axis off;
            title('Maximum response of Steerable Filtering','Fontsize',18);
            
              saveas(h1,[SteerableOutputDir,'/st_',num2str(iFrame),'.tif']);
%              saveas(h1,[SteerableOutputDir,'/st_',num2str(iFrame),'.fig']);
 
             
            % figure 2, for the segmentations
            h2 = figure(2);
            set(h2,'Position',scrsz);
            figure(2); subplot(221);imagesc(currentImg);colormap(gray); axis image; axis off;
            title(['Frame: ',num2str(iFrame)],'Fontsize',18);
            figure(2); subplot(223);imagesc(Intensity_Segment);colormap(gray); axis image; axis off;
            title('Segmentation using intensity','Fontsize',18);
            figure(2); subplot(224);imagesc(SteerabelRes_Segment);colormap(gray); axis image; axis off;
            title('Segmentation using steerable filtering','Fontsize',18);
            figure(2); subplot(222);
            title('Unioned Segmentations','Fontsize',18);
            imagesc(current_seg);colormap(gray); axis image; axis off;
            saveas(h2,[SegmentOutputDir,'/seg_',num2str(iFrame),'.tif']);
%             saveas(h2,[SegmentOutputDir,'/seg_',num2str(iFrame),'.fig']);

%             aa= [800 1070 420 620];
%             subplot(221);axis(aa);
%             subplot(222);axis(aa);
%             subplot(223);axis(aa);
%             subplot(224);axis(aa);
%             
%             aa= [500 800 50 300];
%            subplot(221);axis(aa);
%             subplot(222);axis(aa);
%             subplot(223);axis(aa);
%             subplot(224);axis(aa);
%           
%             aa= [530 800 300 500];
%             subplot(221);axis(aa);
%             subplot(222);axis(aa);
%             subplot(223);axis(aa);
%             subplot(224);axis(aa);
%           
%             aa= [530 800 750 950];
%            subplot(221);axis(aa);
%             subplot(222);axis(aa);
%             subplot(223);axis(aa);
%             subplot(224);axis(aa);           

%            For the direction vector for the orientation.
            A = cos(orienation_map_filtered);
            B = -sin(orienation_map_filtered);

            % figure 3, for direction map, displayed only at segmented
            % target pixel at a grid of 8.
            [X,Y] = meshgrid(1:size(currentImg,2),1:size(currentImg,1));
            A = A(find(current_seg>0));
            B = B(find(current_seg>0));
            X = X(find(current_seg>0));
            Y = Y(find(current_seg>0));
            
            h3 = figure(3); set(h3,'Position',scrsz);
            subplot(121);hold off; imagesc(currentImg);colormap(gray); axis image; axis off;
            title(['Original Image Frame: ',num2str(iFrame)],'Fontsize',18);
            subplot(122);hold off; imagesc(currentImg);colormap(gray); axis image; axis off; hold on;
            
            quiver_grid = 1;
            quiver(X(1:quiver_grid:end,1:quiver_grid:end),Y(1:quiver_grid:end,1:quiver_grid:end), ...
                8*B(1:quiver_grid:end,1:quiver_grid:end,1),8*A(1:quiver_grid:end,1:quiver_grid:end,1));
            
            title('Direction Map','Fontsize',18);
            saveas(h3,[SteerableOutputDir,'/quiver_',num2str(iFrame),'.tif']);
%              saveas(h3,[SteerableOutputDir,'/quiver_',num2str(iFrame),'.fig']);
%             
%             aa= [800 1070 420 620];
%             subplot(121);axis(aa);
%             subplot(122);axis(aa);
%             
%             aa= [500 800 50 300];
%             subplot(121);axis(aa);
%             subplot(122);axis(aa);
%            
%             aa= [530 800 300 500];
%             subplot(121);axis(aa);
%             subplot(122);axis(aa);
%            
%             aa= [530 800 750 950];
%             subplot(121);axis(aa);
%             subplot(122);axis(aa);
%            
%             
%             
        end
              imwrite(current_seg,[SegmentOutputDir,'/segment_',num2str(iFrame),'.tif']);
      
        % Save the results into mat file for each frame, including:
        %       currentImg:     intensity image
        %       orienation_map: orientation( without segmentation information)
        %       MAX_st_res:     steerable filtering response (maximum)
        %       current_seg:    segmentation of VIF
        %       Intensity_Segment: segmentation of VIF only based on intensity
        %       SteerableRes_Segment: segmentation of VIF only based on steerable filtering response
        save([SteerableOutputDir,'/steerable_',num2str(iFrame),'.mat'],...
            'currentImg','orienation_map_filtered','OrientationVoted','orienation_map', ...
            'MAX_st_res', 'current_seg','Intensity_Segment','SteerabelRes_Segment');
        
    end
end
