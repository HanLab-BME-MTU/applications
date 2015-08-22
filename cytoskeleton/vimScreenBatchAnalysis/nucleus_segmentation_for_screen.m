function [movieData, result_flag_matrix, nucleus_number_matrix] = ...
    nucleus_segmentation_for_screen(movieData, plot_save_flag, paramsIn, varargin)
% Main function for the nucleus segmentation for the screen

% Input:     movieData:  movieData object, with the parameters
%            paramsIn:  the parameters to use, if given, overlay
%                    the funParams that comes with movieData. Here are the
%                    fields:
%                      funParams.ChannelIndex:         the channels to process
%                      funParams.Pace_Size:            the parameter to set pace in local segmentation
%                      funParams.Patch_Size:           the parameter to set patch size in local segmentation, for the estimation of local threshold
%                      funParams.lowerbound_localthresholding: The percentage as the lower bound of local thresholding
%                                    local threshold has to be larger or equal to this percentage of the global threshold
%                      funParams.ThresholdWay :         The way to combine segmentation results 
%                      funParams.background_removal:   Flag to do background removal

% Output:    movieData:   updated movieData object. With the segmentation
%                   resulting images saved to the hard disk with corresponing locations.

% Created 2015 by Liya Ding, Matlab R2012b

%% Data Preparation

% before even looking at the nargin, check for the condition and indices
% for the different packages.

% a temp flag for saving tif stack image for Gelfand Lab
save_tif_flag=1;

% Find the package of Filament Analysis
nPackage = length(movieData.packages_);
nProcesses = length(movieData.processes_);

% with no input funparam, use the one the process has on its own
if nargin < 3
    funParams.ChannelIndex = 1;
    if(numel(movieData.channels_)>=3)
         funParams.ChannelIndex = 3;
    end    
    funParams.Pace_Size = 50;
    funParams.Patch_Size = 401;
    funParams.lowerbound=80;
    funParams.ThresholdWay ='OR';
    funParams.background_removal=1;
else
    funParams = paramsIn;
end

result_flag_matrix = nan(numel(movieData.channels_), movieData.nFrames_);
nucleus_number_matrix = nan(numel(movieData.channels_), movieData.nFrames_);


selected_channels = funParams.ChannelIndex;
%% Output Directories

% default steerable filter process output dir
NucleusSegmentationOutputDir = [movieData.outputDirectory_, filesep 'NucleusSegmentation'];

% original:  MD.processes_{2}.outFilePaths_{3}

MD = movieData;
display_msg_flag=1;
package_process_ind_script;

% if there is filamentanalysispackage
if (indexFilamentPackage>0)
    % and a directory is defined for this package
    if (~isempty(movieData.packages_{indexFilamentPackage}.outputDirectory_))
        % and this directory exists
        if (exist(movieData.packages_{indexFilamentPackage}.outputDirectory_,'dir'))
            NucleusSegmentationOutputDir  = [movieData.packages_{indexFilamentPackage}.outputDirectory_, filesep 'NucleusSegmentation'];            
        end
    end
else
    NucleusSegmentationOutputDir = [movieData.outputDirectory_,filesep,'NucleusSegmentation'];    
end


if (~exist(NucleusSegmentationOutputDir,'dir'))
    mkdir(NucleusSegmentationOutputDir);
end

for iChannel = selected_channels
    NucleusSegmentationChannelOutputDir = [NucleusSegmentationOutputDir,filesep,'Channel',num2str(iChannel)];
    if (~exist(NucleusSegmentationChannelOutputDir,'dir'))
        mkdir(NucleusSegmentationChannelOutputDir);
    end
end


%%
nFrame = movieData.nFrames_;
result_flag = nan(numel(movieData.channels_),nFrame);
ImageFlattenFlag=1;
for iChannel =selected_channels
    
    if(indexCellRefineProcess>0)
        MaskrefineChannelOutputDir = movieData.processes_{indexCellRefineProcess}.outFilePaths_{iChannel};    
       else
        MaskrefineChannelOutputDir = NucleusSegmentationChannelOutputDir;
    end
    
    if (~exist(MaskrefineChannelOutputDir,'dir'))
        mkdir(MaskrefineChannelOutputDir);
    end
   
    % Get frame number from the title of the image, this not neccesarily
    % the same as iFrame due to some shorting problem of the channel
    try
        Channel_FilesNames = movieData.channels_(iChannel).getImageFileNames(1:movieData.nFrames_);
    catch
        Channel_FilesNames = movieData.channels_.getImageFileNames(1:movieData.nFrames_);
    end
    filename_short_strs = uncommon_str_takeout(Channel_FilesNames);
        
    %
    display('======================================');
    display(['Current movie: as in ',movieData.outputDirectory_]);
    display(['Start nucleus segmentation in Channel ',num2str(iChannel)]);
    
    % Segment only the real collected data, but skip the padded ones, which
    % were there just to fill in the time lap to make two channel same
    % number of frames
    Sub_Sample_Num = 1 ; % force the sample num to be 1
    Frames_to_Seg = 1:Sub_Sample_Num:nFrame;
    Frames_results_correspondence = im2col(repmat(Frames_to_Seg, [Sub_Sample_Num,1]),[1 1]);
    Frames_results_correspondence = Frames_results_correspondence(1:nFrame);
    
    %     indexFlattenProcess=1;
    for iFrame_index = 1 : length(Frames_to_Seg)
        iFrame = Frames_to_Seg(iFrame_index);
        disp(['Frame: ',num2str(iFrame)]);
        TIC_IC_IF = tic;
        
        % Read in the intensity image.
        if indexFlattenProcess > 0 && ImageFlattenFlag==2
            currentImg = imread([movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel}, filesep, 'flatten_',filename_short_strs{iFrame},'.tif']);
        else
            try
            currentImg = movieData.channels_(iChannel).loadImage(iFrame);
            catch
            currentImg = movieData.channels_.loadImage(iFrame);
            end
        end
        currentImg = double(currentImg);
            
        % this line in commandation for shortest version of filename
        filename_shortshort_strs = all_uncommon_str_takeout(Channel_FilesNames{1});
        
        if  funParams.background_removal > 0 
           currentImg = currentImg - imfilter(currentImg, fspecial('gaussian', 501,200),'replicate','same'); 
        end
        
             
          [level1, Otsu_Segment,level_img_Otsu ] = thresholdLocalSeg(currentImg,'Otsu',funParams.Patch_Size,funParams.Pace_Size,funParams.lowerbound,'showPlots',0);
          [level2, Rosin_Segment,level_img_Rosin ] = thresholdLocalSeg(currentImg,'Rosin',funParams.Patch_Size,funParams.Pace_Size,funParams.lowerbound,'showPlots',0);
     
          level_img_OR = (level_img_Otsu+level_img_Rosin)/2;
          OtsuRosin_Segment = currentImg>level_img_OR;
          
          % if there are lots of saturated area, set the flag as -1;
          % 0 for normal
          if indexFlattenProcess > 0 && ImageFlattenFlag==2        
           saturated_level = double(currentImg).*double(OtsuRosin_Segment) >0.9;
           saturated_img = double(currentImg)>0.9;
          else
           saturated_level = double(currentImg).*double(OtsuRosin_Segment) >3500;   
            saturated_img = double(currentImg)>3500;
         end
          sum(sum(saturated_level))/sum(sum(OtsuRosin_Segment))
          sum(sum(saturated_img))
          
          if(sum(sum(saturated_level))>0.5*sum(sum(OtsuRosin_Segment))) || sum(sum(saturated_img))>100000
              result_flag_matrix(iChannel, iFrame) = -1;
              nucleus_number_matrix(iChannel, iFrame) = -1;         
          else
              result_flag_matrix(iChannel, iFrame) = 0;
              % refine the masl, with separations
              min_size=1000;
              max_num=200;
              close_radius=5;
              OtsuRosin_Segment = ...
                  nucleus_segmentation_refinemask(OtsuRosin_Segment,min_size, max_num, close_radius);
              % count the nubmer of nucleus, with separations
              nucleus_count = nucleus_counting_with_fracs(OtsuRosin_Segment,currentImg,NucleusSegmentationChannelOutputDir, iChannel, iFrame, plot_save_flag);
              nucleus_number_matrix(iChannel, iFrame) = nucleus_count;
          end
          
          
        %% For heat presentation of the segmented filaments
        
        for sub_i = 1 : Sub_Sample_Num
            if iFrame + sub_i-1 <= nFrame
                if(plot_save_flag>0)
                imwrite(Otsu_Segment, ...
                    [NucleusSegmentationChannelOutputDir,filesep,'nucleus_O_segment_',...
                    filename_short_strs{iFrame+ sub_i-1},'.tif']);                
                imwrite(Rosin_Segment, ...
                    [NucleusSegmentationChannelOutputDir,filesep,'nucleus_R_segment_',...
                    filename_short_strs{iFrame+ sub_i-1},'.tif']);   
                imwrite(OtsuRosin_Segment, ...
                    [NucleusSegmentationChannelOutputDir,filesep,'nucleus_OR_segment_',...
                    filename_short_strs{iFrame+ sub_i-1},'.tif']);   
                end
                
                imwrite(OtsuRosin_Segment, ...
                    [MaskrefineChannelOutputDir,filesep,'refined_mask_',...
                    filename_short_strs{iFrame+ sub_i-1},'.tif']);                   
            end
        end
        
        Time_cost = toc(TIC_IC_IF);
        disp(['Frame ', num2str(iFrame), ' nucleus seg costed ',num2str(Time_cost,'%.2f'),'s.']);

    end
end

save([NucleusSegmentationOutputDir,filesep,'result_flag_nucleus_count.mat'],'result_flag_matrix','nucleus_number_matrix');
