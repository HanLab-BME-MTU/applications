% function movieData = filament_segmentation(movieData, varargin)
% Created 07 2012 by Liya Ding, Matlab R2011b

% input movieData object, with the parameters
%   funParams.ChannelIndex:         the channels to process
%   funParams.Pace_Size:            the parameter to set pace in local segmentation
%   funParams.Patch_Size:           the parameter to set patch size in local segmentation, for the estimation of local threshold
%   funParams.lowerbound_localthresholding: The percentage as the lower bound of local thresholding
%                                    local threshold has to be larger or equal to this percentage of the global threshold
%   funParams.Combine_Way :         The way to combine segmentation results from steerable filtering responce
%                                     and from intensity, default is : only use steerable filtering result
%   funParams.Cell_Mask_ind:        Flag to set if cell mask is used, if 1, use segmentation(refined) results,
%                                     if 2, use the user define ROI as in MD_ROI.tif in movieData folder, if 3, no such limit
%   funParams.VIF_Outgrowth_Flag:   Flag to do VIF_outgrowth or not. This is an option made for Gelfand lab

% a temp flag for saving tif stack image for Gelfand Lab
save_tif_flag=1;


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


nProcesses = length(movieData.processes_);

indexFilamentSegmentationProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Filament Segmentation')==1)
        indexFilamentSegmentationProcess = i;
        break;
    end
end

if indexFilamentSegmentationProcess==0
    msgbox('Please set parameters for steerable filtering.')
    return;
end


funParams=movieData.processes_{indexFilamentSegmentationProcess}.funParams_;

selected_channels = funParams.ChannelIndex;
StPace_Size = funParams.StPace_Size;
StPatch_Size = funParams.StPatch_Size;
Stlowerbound =  funParams.st_lowerbound_localthresholding;
IntPace_Size = funParams.IntPace_Size;
IntPatch_Size = funParams.IntPatch_Size;
Intlowerbound =  funParams.int_lowerbound_localthresholding;

Combine_Way = funParams.Combine_Way;
Cell_Mask_ind = funParams.Cell_Mask_ind;
VIF_Outgrowth_Flag = funParams.VIF_Outgrowth_Flag;
Sub_Sample_Num  = funParams.Sub_Sample_Num;

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
indexSteerabeleProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Steerable filtering')==1)
        indexSteerabeleProcess = i;
        break;
    end
end

if indexSteerabeleProcess==0 &&  strcmp(Combine_Way,'st_only')
    msgbox('Please run steerable filtering first.')
    return;
end
ImageFlattenFlag=0;
if indexSteerabeleProcess>0
    
    funParams_st=movieData.processes_{indexSteerabeleProcess}.funParams_;
    
    BaseSteerableFilterSigma = funParams_st.BaseSteerableFilterSigma;
    Levelsofsteerablefilters = funParams_st.Levelsofsteerablefilters;
    ImageFlattenFlag = funParams_st.ImageFlattenFlag;
end

indexFlattenProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Image Flatten')==1)
        indexFlattenProcess = i;
        break;
    end
end

if indexFlattenProcess == 0 && ImageFlattenFlag==2
    display('Please set parameters for Image Flatten.')
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

%% Prepare the cone masks
cone_size = 15;
cone_angle = 25;
cone_mask = cell(180,1);
cone_zero = zeros(2*cone_size+1,2*cone_size+1);
for ci = 1 : 2*cone_size+1
    for cj = 1 : 2*cone_size+1
        cone_zero(ci,cj) = mod(atan2(ci-cone_size-1,cj-cone_size-1),pi);
    end
end
cone_zero_mask = cone_zero<cone_angle/180*pi | cone_zero>pi - cone_angle/180*pi;

cone_zero_mask(cone_size-3:cone_size+3,:)=1;


for cone_i = 1 :180
    cone_mask{cone_i} = imrotate(cone_zero_mask, cone_i, 'nearest','crop');
end

%%
for iChannel = selected_channels
    
    % Get frame number from the title of the image, this not neccesarily
    % the same as iFrame due to some shorting problem of the channel
   Channel_FilesNames = movieData.channels_(iChannel).getImageFileNames(1:movieData.nFrames_);
    
    filename_short_strs = uncommon_str_takeout(Channel_FilesNames);
    
    % Make output directory for the steerable filtered images
    FilamentSegmentationChannelOutputDir =  movieData.processes_{indexFilamentSegmentationProcess}.outFilePaths_{iChannel};
    
    NMS_FilamentSegmentationChannelOutputDir = [FilamentSegmentationChannelOutputDir,filesep,'NMS'];
    ST_FilamentSegmentationChannelOutputDir = [FilamentSegmentationChannelOutputDir,filesep,'ST'];
    
    if (~exist(FilamentSegmentationChannelOutputDir,'dir'))
        mkdir(FilamentSegmentationChannelOutputDir);
    end
    
    HeatOutputDir = [FilamentSegmentationChannelOutputDir,'/HeatOutput'];
    
    if (~exist(HeatOutputDir,'dir'))
        mkdir(HeatOutputDir);
    end
    
    HeatEnhOutputDir = [HeatOutputDir,'/Enh'];
    
    if (~exist(HeatEnhOutputDir,'dir'))
        mkdir(HeatEnhOutputDir);
    end
    
    DataOutputDir = [FilamentSegmentationChannelOutputDir,'/DataOutput'];
     NMS_DataOutputDir = [NMS_FilamentSegmentationChannelOutputDir,'/DataOutput'];
    ST_DataOutputDir = [ST_FilamentSegmentationChannelOutputDir,'/DataOutput'];
    
    if (~exist(DataOutputDir,'dir'))
        mkdir(DataOutputDir);
    end
    
    
    OrientationOutputDir = [FilamentSegmentationChannelOutputDir,'/OrientImage'];
    
    if (~exist(OrientationOutputDir,'dir'))
        mkdir(OrientationOutputDir);
    end
    
    
    % If steerable filter process is run
    if indexSteerabeleProcess>0
        SteerableChannelOutputDir = movieData.processes_{indexSteerabeleProcess}.outFilePaths_{iChannel};
    end
    
    %     if indexFlattenProcess >0
    %         FileNames = movieData.processes_{indexFlattenProcess}.getOutImageFileNames(iChannel);
    %     end
    %
    display('======================================');
    display(['Current movie: as in ',movieData.outputDirectory_]);
    display(['Start filament segmentation in Channel ',num2str(iChannel)]);
    
    % Segment only the real collected data, but skip the padded ones, which
    % were there just to fill in the time lap to make two channel same
    % number of frames
    Frames_to_Seg = 1:Sub_Sample_Num:nFrame;
    Frames_results_correspondence = im2col(repmat(Frames_to_Seg, [Sub_Sample_Num,1]),[1 1]);
    Frames_results_correspondence = Frames_results_correspondence(1:nFrame);
    
%     indexFlattenProcess=1;
    for iFrame_index = 1 : length(Frames_to_Seg)
        iFrame = Frames_to_Seg(iFrame_index);
        for sub_i = 1 : Sub_Sample_Num
            if iFrame + sub_i-1 <= nFrame
        disp(['Frame: ',num2str(iFrame)]);
        
        ST_results = load([ST_DataOutputDir,'/steerable_vote_', ...
            filename_short_strs{iFrame+ sub_i-1},'.mat'],...
            'currentImg','orienation_map_filtered','current_seg');
        
        NMS_results = load([NMS_DataOutputDir,'/steerable_vote_', ...
            filename_short_strs{iFrame+ sub_i-1},'.mat'],...
            'current_seg','current_model');
        
        orienation_map_filtered = ST_results.orienation_map_filtered;
        current_model = NMS_results.current_model;
        currentImg = ST_results.currentImg;
        
        dilate_NMS_seg = imdilate(NMS_results.current_seg, ones(5,5));
        
        
        current_seg = (ST_results.current_seg).*(dilate_NMS_seg);
        
        orienation_map_filtered(current_seg==0)=NaN;
        
            end
        end
        
        currentImg = uint8(currentImg/1);
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
        
        enhanced_im_r(find(current_seg>0))=255*R_seg_orient_heat_map(find(current_seg>0));
        enhanced_im_g(find(current_seg>0))=255*G_seg_orient_heat_map(find(current_seg>0));
        enhanced_im_b(find(current_seg>0))=255*B_seg_orient_heat_map(find(current_seg>0));
        
        RGB_seg_orient_heat_map(:,:,1 ) = enhanced_im_r;
        RGB_seg_orient_heat_map(:,:,2 ) = enhanced_im_g;
        RGB_seg_orient_heat_map(:,:,3 ) = enhanced_im_b;
        
        
        for sub_i = 1 : Sub_Sample_Num
            if iFrame + sub_i-1 <= nFrame
                imwrite(RGB_seg_orient_heat_map, ...
                    [HeatEnhOutputDir,'/segment_heat_',...
                    filename_short_strs{iFrame+ sub_i-1},'.tif']);
            end
        end
        
        
        %% Save segmentation results
        for sub_i = 1 : Sub_Sample_Num
            if iFrame + sub_i-1 <= nFrame
                save([DataOutputDir,'/steerable_vote_', ...
                    filename_short_strs{iFrame+ sub_i-1},'.mat'],...
                    'currentImg','orienation_map_filtered','current_seg',...
                    'current_model', 'RGB_seg_orient_heat_map');
            end
        end
    end
end




