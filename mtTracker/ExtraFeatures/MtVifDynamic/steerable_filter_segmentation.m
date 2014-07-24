function steerable_filter_segmentation(MD, displayflag)
% steerable_filter_segmentation do multiple scale steerable filtering and
% local thresholding for VIF images
%
% SYNOPSIS: steerable_filter_segmentation(MD)
%
% DESCRIPTION: The function gets the image seuqence as in
%               MD.outputDirectory_ with the mask obtained from basic
%               thresholding and post-processing.
%               The fucntion does multiple (4) scale steerable filtering for the
%               image and get the final response of each pixel as the maximum of all
%               responses and set the orientation according to the layer
%               that gives the maximum response.
%
%               Local thresholding is applied both to the intensity image
%               and the steerable filtering response image. The final
%               segmentation is the union of the two segmentation.
%
%               The function creates 2 new directories for the segmentation
%               and steerable filtering results and saves the fig and mat
%               files in these directories.
%
% INPUT  MD :   movieData object ran through process of
%               basic thresholding and post-processing
%        displayflag: a flag to if show and save the figures. If 1, do.
%
% OUTPUT:  None, but save the results into mat file for each frame, including:
%                  currentImg:     intensity image
%                  orienation_map: orientation( without segmentation information)
%                  MAX_st_res:     steerable filtering response (maximum)
%                  current_seg:    segmentation of VIF
%                  Intensity_Segment: segmentation of VIF only based on intensity
%                  SteerableRes_Segment: segmentation of VIF only based on steerable filtering response

% Created Jan 2012 by Liya Ding, Matlab R2011b

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

tic

if (0)%exist([SteerableOutputDir,'/steerable_segment_all_',num2str(nFrame),'.mat'],'file'))
    % Originally attempted to store all the sequences into one big mat file,
    % but this will affect the total memory and speed, so all the comment about
    % the cells are for this reason.
    load([SteerableOutputDir,'/steerable_segment_all_',num2str(nFrame),'.mat']);
else
    %     DirectionMap = cell(1,nFrame);
    %     Segmentation = cell(1,nFrame);
    %     SteerableMaxresponse = cell(1,nFrame);
    %     SteerableMap = cell(1,nFrame);
    %
    for iFrame = 1: nFrame
        
        if (exist([SteerableOutputDir,'/steerable_segment_all_',num2str(nFrame),'.mat'],'file'))
            % load([SteerableOutputDir,'/steerable_',num2str(iFrame),'.mat'],'currentImg','orienation_map','MAX_st_res','current_seg','Intensity_Segment','SteerabelRes_Segment');
            %              Segmentation{1,iFrame} = current_seg;
            %              SteerableMaxresponse{1,iFrame} = MAX_st_res;
            %              SteerableMap{1,iFrame} = orienation_map;
            %              A = cos(orienation_map);
            %              B = sin(orienation_map);
            %              DirectionMap{1,iFrame}(:,:,1) = A;
            %              DirectionMap{1,iFrame}(:,:,2) = B;
        else
            
            %         tic
            
            % Read in the intensity image.
            currentImg = MD.channels_.loadImage(iFrame);
            
            % Get the thresholding results after post-processing
            Mask = MD.processes_{2}.loadChannelOutput(1,iFrame);
            
            % Make the mask bigger in order to include all
            MaskBig = imdilate(Mask, ones(35,35));
            
            % Make another mask to estimate the background intensity
            MaskHuge = imdilate(Mask, ones(55,55));            
            MaskDiff = MaskHuge-Mask;
            
            % The intensity of out of the targetted cell region is set to
            % the mean intensity of the pixels in the ring region around the cell
            % boundary.        
            background_gray = mean(currentImg(find(MaskDiff>0)));
            currentImg = currentImg.*MaskBig + background_gray*(1- MaskBig);
            
            % Steerable filtering using four scales one doubling the
            % previous one.            
            [st1_1,st1_2,st1_3] = steerableDetector(currentImg,2,1.5);
            [st2_1,st2_2,st2_3] = steerableDetector(currentImg,2,3);
            [st3_1,st3_2,st3_3] = steerableDetector(currentImg,2,6);
            [st4_1,st4_2,st4_3] = steerableDetector(currentImg,2,12);
            
            % Stack the results in big 3D matrice
            ST_res_all = zeros(size(st1_1,1),size(st1_1,2),4);
            ST_res_all(:,:,1) = st1_1;
            ST_res_all(:,:,2) = st2_1;
            ST_res_all(:,:,3) = st3_1;
            ST_res_all(:,:,4) = st4_1;
            
            orienation_maps_all = zeros(size(st1_1,1),size(st1_1,2),4);
            orienation_maps_all(:,:,1) = st1_2;
            orienation_maps_all(:,:,2) = st2_2;
            orienation_maps_all(:,:,3) = st3_2;
            orienation_maps_all(:,:,4) = st4_2;
            
            % Takes the maximum from the all the responses
            [MAX_st_res, IND] = max(ST_res_all,[],3);
            
            % Use the orientation according to the maxmium response (for each pixel)
            for i = 1 : size(currentImg,1)
                for j = 1 : size(currentImg,2)
                    orienation_map(i,j) = orienation_maps_all(i,j,IND(i,j));
                end
            end
            
            % Local thresholding for both intensity and steerabel filtering
            % response.
            %tic
            [level1, Intensity_Segment ] = thresholdOtsu_local(MAX_st_res,15,3,0);
            [level2, SteerabelRes_Segment ] = thresholdOtsu_local(currentImg,15,3,0);
            %toc
            
            % The segmentation is set as the union of two segmentation.
            current_seg = or(Intensity_Segment,SteerabelRes_Segment);
            
            % For the direction vector for the orientation.
            A = cos(orienation_map);
            B = sin(orienation_map);
            %         DirectionMap{1,iFrame}(:,:,1) = A;
            %         DirectionMap{1,iFrame}(:,:,2) = B;
                        
            if(displayflag==1)
            % If want to display and save images
                % figure 1, for image and the responses
                h1 = figure(1);subplot(231);                
                imagesc(currentImg);colormap(gray); axis image; axis off;
                title(['Frame: ',num2str(iFrame)]);
                subplot(232);imagesc(st1_1);colormap(gray); axis image; axis off;
                title('Steerable Filtering with Scale 1.5');
                subplot(233);imagesc(st2_1);colormap(gray); axis image; axis off;
                title('Scale 3');
                subplot(234);imagesc(st3_1);colormap(gray); axis image; axis off;
                title('Scale 6');
                subplot(235);imagesc(st4_1);colormap(gray); axis image; axis off;
                title('Scale 12');
                subplot(236);imagesc(MAX_st_res);colormap(gray); axis image; axis off;
                title('Maximum response');
                saveas(h1,[SteerableOutputDir,'/st_',num2str(iFrame),'.tif']);
                
                % figure 2, for the segmentations
                h2= figure(2);
                figure(2); subplot(221);imagesc(currentImg);colormap(gray); axis image; axis off;
                title(['Frame: ',num2str(iFrame)]);
                figure(2); subplot(223);imagesc(Intensity_Segment);colormap(gray); axis image; axis off;
                title('Segmentation using intensity');
                figure(2); subplot(224);imagesc(SteerabelRes_Segment);colormap(gray); axis image; axis off;
                title('Segmentation using steerable filtering');
                figure(2); subplot(222);
                title('Unioned Segmentations');
                imagesc(current_seg);colormap(gray); axis image; axis off;
                saveas(h2,[SegmentOutputDir,'/seg_',num2str(iFrame),'.tif']);
                
                % figure 3, for direction map, displayed only at segmented
                % target pixel at grid 10.
                [X,Y] = meshgrid(1:size(currentImg,2),1:size(currentImg,1));
                A = A(find(current_seg>0));
                B = B(find(current_seg>0));
                X = X(find(current_seg>0));
                Y = Y(find(current_seg>0));
                
                h3 = figure(3); subplot(121);hold off; imagesc(currentImg);colormap(gray); axis image; axis off;
                title(['Frame: ',num2str(iFrame)]);
                subplot(122);hold off; imagesc(currentImg);colormap(gray); axis image; axis off; hold on;
                quiver(X(1:10:end,1:10:end),Y(1:10:end,1:10:end),8*B(1:10:end,1:10:end,1),8*A(1:10:end,1:10:end,1));
                title('Direction Map');
                saveas(h3,[SteerableOutputDir,'/quiver_',num2str(iFrame),'.tif']);
                saveas(h3,[SteerableOutputDir,'/quiver_',num2str(iFrame),'.fig']);
                imwrite(current_seg,[SegmentOutputDir,'/segment_',num2str(iFrame),'.tif']);
            end
            
            % Save the results into mat file for each frame, including:
            %       currentImg:     intensity image
            %       orienation_map: orientation( without segmentation information)
            %       MAX_st_res:     steerable filtering response (maximum)
            %       current_seg:    segmentation of VIF
            %       Intensity_Segment: segmentation of VIF only based on intensity
            %       SteerableRes_Segment: segmentation of VIF only based on steerable filtering response            
            save([SteerableOutputDir,'/steerable_',num2str(iFrame),'.mat'],...
                'currentImg','orienation_map','MAX_st_res','current_seg','Intensity_Segment','SteerabelRes_Segment');
            
            %         Segmentation{1,iFrame} = current_seg;
            %         SteerableMaxresponse{1,iFrame} = MAX_st_res;
            %         SteerableMap{1,iFrame} = orienation_map;
            %
            %         toc
        end
        % save([SteerableOutputDir,'/steerable_segment_all_',num2str(iFrame),'.mat'],'Segmentation','SteerableMaxresponse','SteerableMap','DirectionMap');        
    end
end
toc