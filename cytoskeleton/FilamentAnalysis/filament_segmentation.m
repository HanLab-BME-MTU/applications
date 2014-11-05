function movieData = filament_segmentation(movieData, paramsIn, wholemovie_input_filename, varargin)
% Main function for the segmentation/reconstruction step in filament analysis

% Input:     movieData:  movieData object, with the parameters
%            funParams:  the parameters to use, if given, overlay
%                    the funParams that comes with movieData. Here are the
%                    fields:
%                      funParams.ChannelIndex:         the channels to process
%                      funParams.Pace_Size:            the parameter to set pace in local segmentation
%                      funParams.Patch_Size:           the parameter to set patch size in local segmentation, for the estimation of local threshold
%                      funParams.lowerbound_localthresholding: The percentage as the lower bound of local thresholding
%                                    local threshold has to be larger or equal to this percentage of the global threshold
%                      funParams.Combine_Way :         The way to combine segmentation results from steerable filtering responce
%                                     and from intensity, default is : only use steerable filtering result
%                      funParams.Cell_Mask_ind:        Flag to set if cell mask is used, if 1, use segmentation(refined) results,
%                                     if 2, use the user define ROI as in MD_ROI.tif in movieData folder, if 3, no such limit
%                      funParams.VIF_Outgrowth_Flag:   Flag to do VIF_outgrowth or not. This is an option made for Gelfand lab
%            wholemovie_input_filename: the filename of mat file previously saved for the
%                           whole movie statistics for this or some other
%                           movie. If given, this will overwrite the
%                           wholemovie statistics that comes with the
%                           movieData

% Output:    movieData:   updated movieData object. With the segmentation
%                   resulting images saved to the hard disk with corresponing locations.

% Created 07 2012 by Liya Ding, Matlab R2011b

%% Data Preparation

% before even looking at the nargin, check for the condition and indices
% for the different packages.

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


% with no input funparam, use the one the process has on its own
if nargin < 2
    paramsIn = [];
    funParams = movieData.processes_{indexFilamentSegmentationProcess}.funParams_;
else
    funParams = paramsIn;
end


selected_channels = funParams.ChannelIndex;

StPace_Size_movie = funParams.StPace_Size;
StPatch_Size_movie = funParams.StPatch_Size;
st_lowerbound_localthresholding_movie =  funParams.st_lowerbound_localthresholding;
IntPace_Size_movie = funParams.IntPace_Size;
IntPatch_Size_movie = funParams.IntPatch_Size;
int_lowerbound_localthresholding_movie =  funParams.int_lowerbound_localthresholding;

Combine_Way_movie = funParams.Combine_Way;
Cell_Mask_ind_movie = funParams.Cell_Mask_ind;
VIF_Outgrowth_Flag_movie = funParams.VIF_Outgrowth_Flag;
Sub_Sample_Num_movie  = funParams.Sub_Sample_Num;
Whole_movie_ind_movie  = funParams.Whole_movie_ind;
Rerun_WholeMovie =  funParams.Rerun_WholeMovie;

SaveFigures_movie = funParams.savestepfigures;
ShowDetailMessages_movie = funParams.savestepfigures;

CoefAlpha_movie = funParams.CoefAlpha;
LengthThreshold_movie = funParams.LengthThreshold;
IternationNumber_movie = funParams.IternationNumber;
CurvatureThreshold_movie = funParams.CurvatureThreshold;
Cell_Mask_ind = Cell_Mask_ind_movie;

% if the package was run with old version with only one set of parameter
% for all channels, make a copy of parameter to every channel
if(length(funParams.StPace_Size)==1 && numel(movieData.channels_)>1)
    ones_array = ones(1,numel(movieData.channels_));
    funParams.StPace_Size = funParams.StPace_Size*ones_array;
    funParams.StPatch_Size = funParams.StPatch_Size*ones_array;
    funParams.st_lowerbound_localthresholding = funParams.st_lowerbound_localthresholding*ones_array;
    funParams.IntPace_Size = funParams.IntPace_Size*ones_array;
    funParams.IntPatch_Size = funParams.IntPatch_Size*ones_array;
    funParams.int_lowerbound_localthresholding = funParams.int_lowerbound_localthresholding*ones_array;
    funParams.Cell_Mask_ind = funParams.Cell_Mask_ind*ones_array;
    funParams.Whole_movie_ind = funParams.Whole_movie_ind*ones_array;
    
    Combine_Way = funParams.Combine_Way;
    
    funParams.Combine_Way=cell(1,1);
    
    for iC = 1 : numel(movieData.channels_)
        funParams.Combine_Way{iC}= Combine_Way;
    end
    
    funParams.Classifier_Type_ind = funParams.Classifier_Type_ind*ones_array;
    funParams.LengthThreshold = funParams.LengthThreshold*ones_array;
    funParams.CurvatureThreshold = funParams.CurvatureThreshold*ones_array;
    funParams.IternationNumber = funParams.IternationNumber*ones_array;
    funParams.CoefAlpha = funParams.CoefAlpha*ones_array;
    funParams.training_sample_number = funParams.training_sample_number*ones_array;
end


%% Output Directories

% default steerable filter process output dir
FilamentSegmentationProcessOutputDir = [movieData.outputDirectory_, filesep 'FilamentSegmentation'];

% if there is filamentanalysispackage
if (indexFilamentPackage>0)
    % and a directory is defined for this package
    if (~isempty(movieData.packages_{indexFilamentPackage}.outputDirectory_))
        % and this directory exists
        if (exist(movieData.packages_{indexFilamentPackage}.outputDirectory_,'dir'))
            FilamentSegmentationProcessOutputDir  = [movieData.packages_{indexFilamentPackage}.outputDirectory_, filesep 'FilamentSegmentation'];
        end
    end
end

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

if indexSteerabeleProcess==0 && Combine_Way~=2
    msgbox('Please run steerable filtering first.')
    return;
end

funParams_st=movieData.processes_{indexSteerabeleProcess}.funParams_;

BaseSteerableFilterSigma = funParams_st.BaseSteerableFilterSigma;
Levelsofsteerablefilters = funParams_st.Levelsofsteerablefilters;
ImageFlattenFlag = funParams_st.ImageFlattenFlag;

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

% if user want to use an input whole movie stat result, use it
if nargin >=3
    load(wholemovie_input_filename);
    funParams.Whole_movie_stat_cell = Whole_movie_stat_cell;
else
    % or, calculate it
    %% May 1st 2014, due to change in flattening precedure, this whole movie stat need rerun,
    % if there is already whole movie file and the user didn't ask for
    % rerun the whole movie
    if(exist([FilamentSegmentationProcessOutputDir, filesep, 'whole_movie_stat.mat'],'file')>0 ...
            && Rerun_WholeMovie==0)
        load([FilamentSegmentationProcessOutputDir, filesep, 'whole_movie_stat.mat'],'Whole_movie_stat_cell');
    else
        Whole_movie_stat_cell = whole_movie_stat_function(movieData);
        save([FilamentSegmentationProcessOutputDir, filesep, 'whole_movie_stat.mat'],'Whole_movie_stat_cell');
    end
    
    funParams.Whole_movie_stat_cell = Whole_movie_stat_cell;
    
end

%%
nFrame = movieData.nFrames_;

% If the user set an cell ROI read in
if(exist([movieData.outputDirectory_,filesep,'MD_ROI.tif'],'file'))
    user_input_mask = imread([movieData.outputDirectory_,filesep,'MD_ROI.tif']);
end

%% cones related is not in use
% %% Prepare the cone masks
% cone_size = 15;
% cone_angle = 25;
% cone_mask = cell(180,1);
% cone_zero = zeros(2*cone_size+1,2*cone_size+1);
% for ci = 1 : 2*cone_size+1
%     for cj = 1 : 2*cone_size+1
%         cone_zero(ci,cj) = mod(atan2(ci-cone_size-1,cj-cone_size-1),pi);
%     end
% end
% cone_zero_mask = cone_zero<cone_angle/180*pi | cone_zero>pi - cone_angle/180*pi;
%
% cone_zero_mask(cone_size-3:cone_size+3,:)=1;
%
%
% for cone_i = 1 :180
%     cone_mask{cone_i} = imrotate(cone_zero_mask, cone_i, 'nearest','crop');
% end

%%
for iChannel = selected_channels
    
    if(length(StPace_Size_movie)>=iChannel)
        % in new setting, each channel can have its own setting, different
        % from another channel
        StPace_Size =  StPace_Size_movie(iChannel);
        StPatch_Size = StPatch_Size_movie(iChannel);
        Stlowerbound =  st_lowerbound_localthresholding_movie(iChannel);
        IntPace_Size = IntPace_Size_movie(iChannel);
        IntPatch_Size = IntPatch_Size_movie(iChannel);
        Intlowerbound =  int_lowerbound_localthresholding_movie(iChannel);
        if(~iscell(Combine_Way_movie))
            Combine_Way = Combine_Way_movie;
        else
            Combine_Way = Combine_Way_movie{iChannel};
        end
        Cell_Mask_ind = Cell_Mask_ind_movie(iChannel);
    else
        % in the original situation there is one common setting for all
        % channels, if the MD loaded is in this case, use the same
        StPace_Size =  StPace_Size_movie;
        StPatch_Size = StPatch_Size_movie;
        Stlowerbound =  st_lowerbound_localthresholding_movie;
        IntPace_Size = IntPace_Size_movie;
        IntPatch_Size = IntPatch_Size_movie;
        Intlowerbound =  int_lowerbound_localthresholding_movie;
        Combine_Way = Combine_Way_movie;
        Cell_Mask_ind = Cell_Mask_ind_movie;
    end
    
    VIF_Outgrowth_Flag = VIF_Outgrowth_Flag_movie;
    Sub_Sample_Num  = Sub_Sample_Num_movie;
    
    
    % Get frame number from the title of the image, this not neccesarily
    % the same as iFrame due to some shorting problem of the channel
    Channel_FilesNames = movieData.channels_(iChannel).getImageFileNames(1:movieData.nFrames_);
    
    filename_short_strs = uncommon_str_takeout(Channel_FilesNames);
    
    % Make output directory for the steerable filtered images
    FilamentSegmentationChannelOutputDir =  movieData.processes_{indexFilamentSegmentationProcess}.outFilePaths_{iChannel};
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
        
        disp(['Frame: ',num2str(iFrame)]);
        
        % Read in the intensity image.
        if indexFlattenProcess > 0 && ImageFlattenFlag==2
            currentImg = imread([movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel}, filesep, 'flatten_',filename_short_strs{iFrame},'.tif']);
        else
            currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        end
        currentImg = double(currentImg);
        
        %% %tif stack cost too much memory, comment these
        %
        %         if( save_tif_flag==1 && iFrame==Frames_to_Seg(1)  )
        %             % Gelfand lab needs single file results for tif stack
        %             tif_stack_binary_seg_image_data = uint8(zeros(size(currentImg,1),size(currentImg,2),length(Frames_to_Seg)));
        %             tif_stack_RGB_heat_image_data = uint8(zeros(size(currentImg,1),size(currentImg,2),3,length(Frames_to_Seg)));
        %         end
        %%
        
        % this line in commandation for shortest version of filename
        filename_shortshort_strs = all_uncommon_str_takeout(Channel_FilesNames{1});
        
        try
            load([SteerableChannelOutputDir, filesep, 'steerable_',...
                filename_short_strs{iFrame},'.mat']);
        catch
            % in the case of only having the short-old version
            load([SteerableChannelOutputDir, filesep, 'steerable_',...
                filename_shortshort_strs{iFrame},'.mat']);
        end
        
        
        %%
        
        
        MaskCell = ones(size(currentImg));
        
        if Cell_Mask_ind == 1 % using cell segmentation from same channel
            MaskCell = movieData.processes_{indexCellSegProcess}.loadChannelOutput(iChannel,iFrame);
        else
            if Cell_Mask_ind == 2 % Using input static ROI tiff
                MaskCell = user_input_mask>0;
            else
                if Cell_Mask_ind == 5 % No limit
                    MaskCell = ones(size(currentImg,1),size(currentImg,2));
                else
                    if Cell_Mask_ind == 4 % Combine from both channel directly
                        MaskVIFCell = movieData.processes_{indexCellSegProcess}.loadChannelOutput(2,iFrame);
                        MaskMTCell = movieData.processes_{indexCellSegProcess}.loadChannelOutput(1,iFrame);
                        MaskCell = MaskVIFCell | MaskMTCell;
                        
                    else
                        % Combine from both channel
                        % In this option, the channel need to be 1. MT or Membrame, 2. VIF or Actin
                        MaskVIFCell = movieData.processes_{indexCellSegProcess}.loadChannelOutput(2,iFrame);
                        MaskMTCell = movieData.processes_{indexCellSegProcess}.loadChannelOutput(1,iFrame);
                        
                        H_close_cell = fspecial('disk',5);
                        H_close_cell = H_close_cell>0;
                        
                        MaskMTCell = imerode(MaskMTCell,H_close_cell);
                        TightMask = MaskVIFCell.*MaskMTCell;
                        
                        % Make the mask bigger in order to include all
                        MaskCell = imdilate(TightMask, ones(15,15),'same');
                        
                        clearvars MaskVIFCell MaskMTCell TightMask H_close_cell;
                    end
                end
            end
        end
        
        
        
        %%
        % Correcting the nms ending semicircle due to the aritifact of
        % simulation perfect Gaussian ends
        
        %
        %         median_nms = median(nms(:));
        %         quarter_nms = median(nms(find(nms>median_nms)));
        
        open_nms = imopen(nms, ones(2,2));
        eroded_nms = nms-open_nms;
        
        %use the opened image;
        nms = eroded_nms;
        
        stophere=1;
        
        
        %%
        
        
        NMS_Segment=[];
        Intensity_Segment=[];
        SteerabelRes_Segment=[];
        
        Min_area = 20;
        Min_longaxis = 6;
        
        
        switch Combine_Way
            case 'int_st_both'
                level0 = thresholdOtsu(MAX_st_res);
                thresh_Segment = MAX_st_res > level0;
                % all the local seg was changed back to one without the
                % last parameter; could be a version conflict
                %                 [level1, SteerabelRes_Segment ] = thresholdLocalSeg(MAX_st_res,'Otsu',StPatch_Size,StPace_Size,Stlowerbound,0,Whole_movie_stat_cell{iChannel}.otsu_ST);
                %                 [level2, Intensity_Segment ] = thresholdLocalSeg(currentImg,'Otsu',IntPatch_Size,IntPace_Size,Intlowerbound,0,Whole_movie_stat_cell{iChannel}.otsu_INT);
                [level1, SteerabelRes_Segment ] = thresholdLocalSeg(MAX_st_res,'Otsu',StPatch_Size,StPace_Size,Stlowerbound,'showPlots',0);
                [level2, Intensity_Segment ] = thresholdLocalSeg(currentImg,'Otsu',IntPatch_Size,IntPace_Size,Intlowerbound,'showPlots',0);
                current_seg = and(Intensity_Segment,SteerabelRes_Segment);
                current_model=[];
                
            case 'st_only'
                %                 [level1, SteerabelRes_Segment ] = thresholdLocalSeg(MAX_st_res,'Otsu',StPatch_Size,StPace_Size,Stlowerbound,0,Whole_movie_stat_cell{iChannel}.otsu_ST);
                [level1, SteerabelRes_Segment ] = thresholdLocalSeg(MAX_st_res,'Otsu',StPatch_Size,StPace_Size,Stlowerbound,'showPlots',0);
                current_seg = SteerabelRes_Segment;
                Intensity_Segment = current_seg;
                SteerabelRes_Segment = current_seg;
                current_model=[];
                
            case 'st_nms_two'
                %                 [level1, SteerabelRes_Segment ] = thresholdLocalSeg(MAX_st_res,'Otsu',StPatch_Size,StPace_Size,Stlowerbound*0.7,0,Whole_movie_stat_cell{iChannel}.otsu_ST);
                %                 [level2, NMS_Segment ] = thresholdLocalSeg(nms,'Rosin',StPatch_Size,StPace_Size,Stlowerbound*1.3,0,Whole_movie_stat_cell{iChannel}.otsu_NMS);
                [level1, SteerabelRes_Segment ] = thresholdLocalSeg(MAX_st_res,'Otsu',StPatch_Size,StPace_Size,Stlowerbound*0.7,'showPlots',0);
                [level2, NMS_Segment ] = thresholdLocalSeg(nms,'Rosin',StPatch_Size,StPace_Size,Stlowerbound*1.3,showPlots,0);
                current_seg = imdilateWithScale(NMS_Segment,scaleMap,BaseSteerableFilterSigma.*(2.^((1:Levelsofsteerablefilters)-1)))...
                    .*SteerabelRes_Segment;
                
                Intensity_Segment = current_seg;
                SteerabelRes_Segment = current_seg;
                current_model=[];
                
            case 'st_nms_only'
                %                 [level2, NMS_Segment ] = thresholdLocalSeg(nms,'Rosin',StPatch_Size,StPace_Size,Stlowerbound,0,Whole_movie_stat_cell{iChannel}.otsu_NMS);
                [level2, NMS_Segment ] = thresholdLocalSeg(nms,'Rosin',StPatch_Size,StPace_Size,Stlowerbound,'showPlots',0);
                current_seg = NMS_Segment;
                Intensity_Segment = current_seg;
                SteerabelRes_Segment = current_seg;
                Min_area = 6;
                current_model=[];
                
            case 'geo_based_training'
                
                if(~isempty(funParams.F_classifier{iChannel}))
                    load(funParams.F_classifier{iChannel});
                else
                    F_classifer_train_this_channel=[];
                end
                
                [level2, NMS_Segment,current_model ] = ...
                    geoBasedNmsSeg_withtraining(nms,currentImg, F_classifer_train_this_channel,1,...
                    MaskCell,iFrame,FilamentSegmentationChannelOutputDir,funParams);
                current_seg = NMS_Segment;
                Intensity_Segment = current_seg;
                SteerabelRes_Segment = current_seg;
                Min_area = 6;
                
            case 'reserved_for_test'
                
                % Liya: test for the running of the comparison
                current_seg_canny_cell=cell(1,1);
                display('Canny Test:');
                tic
                %                 for PercentOfPixelsNotEdges = 0.8: 0.1: 0.95
                %                     for ThresholdRatio = 0.8 : 0.1: 0.95
                for iP = 1 : 5
                    for iT = 1 : 5
                        
                        PercentOfPixelsNotEdges = iP/20+0.70;
                        ThresholdRatio = iT/20+0.70;
                        
                        [lowThresh, highThresh, current_seg]...
                            = proximityBasedNmsSeg(MAX_st_res,...
                            orienation_map,funParams,...
                            PercentOfPixelsNotEdges,ThresholdRatio);
                        current_seg_canny_cell{iP,iT} = current_seg;
                        
                        
                    end
                end
                toc
                
                % Assume no training is done for the classifier, so use the
                % linear plane classifier with the input parameters.
                
                %                 if(~isempty(funParams.F_classifier{iChannel}))
                %                     load(funParams.F_classifier{iChannel});
                %                 else
                F_classifer_train_this_channel=[];
                %                 end
                
                display(['Geo based GM Frame',num2str(iFrame),':']);
                tic
                [level2, NMS_Segment,current_model ] = ...
                    geoBasedNmsSeg_withGM(nms,currentImg, F_classifer_train_this_channel,1,...
                    MaskCell,iFrame,FilamentSegmentationChannelOutputDir,funParams,iChannel);
                toc
                
                current_seg = NMS_Segment;
                Intensity_Segment = current_seg;
                SteerabelRes_Segment = current_seg;
                Min_area = 6;
                
                
                %% for screen no graph matching, no curvature/intensity control
            case 'geo_based_no_GM'
                F_classifer_train_this_channel=[];
                %                 end
                
                display(['Geo based NO-GM, Frame',num2str(iFrame),':']);
                tic
                [level2, NMS_Segment,current_model ] = ...
                    geoBasedNmsSeg_withoutGM(nms,currentImg, F_classifer_train_this_channel,1,...
                    MaskCell,iFrame,FilamentSegmentationChannelOutputDir,funParams,iChannel);
                toc
                
                current_seg = NMS_Segment;
                Intensity_Segment = current_seg;
                SteerabelRes_Segment = current_seg;
                Min_area = 6;
                
            case 'geo_based_GM'
                
                % Assume no training is done for the classifier, so use the
                % linear plane classifier with the input parameters.
                
                %                 if(~isempty(funParams.F_classifier{iChannel}))
                %                     load(funParams.F_classifier{iChannel});
                %                 else
                F_classifer_train_this_channel=[];
                %                 end
                
                display(['Geo based GM Frame',num2str(iFrame),':']);
                tic
                [level2, NMS_Segment,current_model ] = ...
                    geoBasedNmsSeg_withGM(nms,currentImg, F_classifer_train_this_channel,1,...
                    MaskCell,iFrame,FilamentSegmentationChannelOutputDir,funParams,iChannel);
                toc
                
                current_seg = NMS_Segment;
                Intensity_Segment = current_seg;
                SteerabelRes_Segment = current_seg;
                Min_area = 6;
                
            case 'int_only'
                %                 [level2, Intensity_Segment ] = thresholdLocalSeg(currentImg,'Otsu',IntPatch_Size,IntPace_Size,Intlowerbound,0);
                [level2, Intensity_Segment ] = thresholdLocalSeg(currentImg,'Otsu',IntPatch_Size,IntPace_Size,Intlowerbound,'showPlots',0);
                
                current_seg = Intensity_Segment;
                SteerabelRes_Segment = current_seg;
                current_model=[];
                
            case 'canny_method'                
                
                % Liya: test for the running of the comparison
                current_seg_canny_cell=cell(1,1);
                display('Canny Test:');
%                 tic
%                 %                 for PercentOfPixelsNotEdges = 0.8: 0.1: 0.95
%                 %                     for ThresholdRatio = 0.8 : 0.1: 0.95
%                 for iP = 1 : 5
%                     for iT = 1 : 5
%                         
%                         PercentOfPixelsNotEdges = iP/20+0.70;
%                         ThresholdRatio = iT/20+0.70;
%                         
%                         [lowThresh, highThresh, current_seg]...
%                             = proximityBasedNmsSeg(MAX_st_res,...
%                             orienation_map,funParams,...
%                             PercentOfPixelsNotEdges,ThresholdRatio);
%                         current_seg_canny_cell{iP,iT} = current_seg;
%                         
%                         
%                     end
%                 end
%                 toc
                
                % get the percentage threshold from funParam
                HigherThresdhold = funParams.CannyHigherThreshold(iChannel)/100;
                LowerThresdhold = funParams.CannyLowerThreshold(iChannel)/100;
                
                tic
                [lowThresh, highThresh, current_seg]...
                    = proximityBasedNmsSeg(MAX_st_res,orienation_map,funParams,HigherThresdhold,LowerThresdhold);
                toc
                
                level2 = highThresh;
                Intensity_Segment = current_seg;
                SteerabelRes_Segment = current_seg;
                Min_area = 6;
                current_model=[];
                
            otherwise
                warning('Use the default of union');
                %                    [level1, SteerabelRes_Segment ] = thresholdLocalSeg(MAX_st_res,'Otsu',StPatch_Size,StPace_Size,Stlowerbound,0,Whole_movie_stat_cell{iChannel}.otsu_ST);
                %                 [level2, Intensity_Segment ] = thresholdLocalSeg(currentImg,'Otsu',IntPatch_Size,IntPace_Size,Intlowerbound,0,Whole_movie_stat_cell{iChannel}.otsu_INT);
                [level1, SteerabelRes_Segment ] = thresholdLocalSeg(MAX_st_res,'Otsu',StPatch_Size,StPace_Size,Stlowerbound,'showPlots',0);
                [level2, Intensity_Segment ] = thresholdLocalSeg(currentImg,'Otsu',IntPatch_Size,IntPace_Size,Intlowerbound,'showPlots',0);
                % The segmentation is set as the union of two segmentation.
                current_seg = or(Intensity_Segment,SteerabelRes_Segment);
                current_model=[];
        end
        
        MaskCell=MaskCell>0;
        current_seg = current_seg.*MaskCell;
        
        
        %%
        % A smoothing done only at the steerable filtering results, if only intensity only, then the same
        orienation_map_filtered = OrientationSmooth(orienation_map, SteerabelRes_Segment);
        
        %%
        % % Voting of the orientation field for the non-steerable filter
        % % segmented places.
        
        % % the voting is not in use
        % OrientationVoted = OrientationVote(orienation_map,SteerabelRes_Segment,3,45);
        OrientationVoted= orienation_map_filtered;
        
        intensity_addon = current_seg - SteerabelRes_Segment ==1;
        if (~isempty(max(max(intensity_addon))>0))
            orienation_map_filtered(find(intensity_addon>0)) = OrientationVoted(find(intensity_addon>0));
        end
        
        if(~strcmp(Combine_Way,'geo_based'))
            % if the segmentation is not done with geo_based method, do
            % some geometry based checking on the results
            
            %% Deleting the small isolated dots
            labelMask = bwlabel(current_seg);
            ob_prop = regionprops(labelMask,'Area','MajorAxisLength','Eccentricity','MinorAxisLength');
            
            obAreas = [ob_prop.Area];
            obLongaxis = [ob_prop.MajorAxisLength];
            obShortaxis = [ob_prop.MinorAxisLength];
            obEccentricity = [ob_prop.Eccentricity];
            ratio  = obShortaxis./obLongaxis;
            % for now the parameters are set here, without adaptiveness
            for i_area = 1 : length(obAreas)
                if obAreas(i_area) < 100
                    angle_area{i_area} = orienation_map_filtered(find(labelMask==i_area));
                    [h_area, bin] = hist(angle_area{i_area},-pi/2:5/180*pi:pi/2);
                    ind_t = find(h_area==max(h_area));
                    temp = mod((angle_area{i_area} - bin(ind_t(1)) + pi/2), pi) - pi/2;
                    if std(temp)>0.75 && max(h_area)<0.2*length(angle_area{i_area}) && ratio(i_area) >0.5 && obLongaxis(i_area)<20
                        labelMask(find(labelMask==i_area))=0;
                    end
                end
            end
            
            current_seg = labelMask > 0;
            
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
                    
                    if obAreas(i_area) <Min_area || obLongaxis(i_area) <Min_longaxis
                        labelMask(labelMask==i_area) = 0;
                    end
                    
                end
            end
            
            current_seg = labelMask > 0;
        end
        
        %
        %         [ind_a,ind_b] = find(current_seg>0);
        %
        %         cone_bins = cell(size(current_seg,1), size(current_seg,2));
        %         cone_weight_bins = cell(size(current_seg,1), size(current_seg,2));
        %
        %         for si = 1 : length(ind_a)
        %             pixel_angle = round((-orienation_map(ind_a(si), ind_b(si))+pi/2)*180/pi);
        %             weight_res = MAX_st_res(ind_a(si), ind_b(si));
        %             if pixel_angle ==0
        %                 pixel_angle = 180;
        %             end
        %             try
        %                 [ind_c,ind_d] = find(cone_mask{pixel_angle}>0);
        %                 for p_i = 1 : length(ind_c)
        %                     cone_bins{ind_c(p_i)+ind_a(si)-cone_size-1, ind_d(p_i)+ind_b(si)-cone_size-1} ...
        %                         = [cone_bins{ind_c(p_i)+ind_a(si)-cone_size-1, ind_d(p_i)+ind_b(si)-cone_size-1} pixel_angle];
        %                     cone_weight_bins{ind_c(p_i)+ind_a(si)-cone_size-1, ind_d(p_i)+ind_b(si)-cone_size-1} ...
        %                         = [cone_bins{ind_c(p_i)+ind_a(si)-cone_size-1, ind_d(p_i)+ind_b(si)-cone_size-1} weight_res];
        %                 end
        %             end
        %
        %         end
        %
        %         for p_i = 1 : size(cone_bins,1)
        %             for p_j = 1 : size(cone_bins,2)
        %                 length_cone_bin(p_i,p_j) = length(cone_bins{p_i,p_j});
        %                  weight_bin(p_i,p_j) = m(cone_weight_bins{p_i,p_j});
        %                 if(~isempty(cone_bins{p_i,p_j}))
        %                 h = hist(cone_bins{p_i,p_j},0:10:180);
        %                 centernumber_cone_bin(p_i,p_j) = max(h);
        %                 else
        %                      centernumber_cone_bin(p_i,p_j) =0;
        %                 end
        %             end
        %         end
        
        
        %% For heat presentation of the segmented filaments
        
        
        for sub_i = 1 : Sub_Sample_Num
            if iFrame + sub_i-1 <= nFrame
                imwrite(current_seg, ...
                    [FilamentSegmentationChannelOutputDir,'/segment_binary_',...
                    filename_short_strs{iFrame+ sub_i-1},'.tif']);
                imwrite(orienation_map_filtered.*double(current_seg), ...
                    [OrientationOutputDir,'/segment_orientation_',...
                    filename_short_strs{iFrame+ sub_i-1},'.tif']);
            end
        end
        
        
        if(~isempty(current_model))
            [Vif_digital_model,Vif_orientation_model,VIF_XX,VIF_YY,VIF_OO] ...
                = filament_model_to_digital_with_orientation(current_model);
            
            OO_flip = pi-VIF_OO;
            
            OO_flip(OO_flip<-pi/2)=OO_flip(OO_flip<-pi/2)+pi;
            OO_flip(OO_flip>pi/2)=OO_flip(OO_flip>pi/2)-pi;
            OO_flip(OO_flip<-pi/2)=OO_flip(OO_flip<-pi/2)+pi;
            OO_flip(OO_flip>pi/2)=OO_flip(OO_flip>pi/2)-pi;
            OO_flip(OO_flip<-pi/2)=OO_flip(OO_flip<-pi/2)+pi;
            OO_flip(OO_flip>pi/2)=OO_flip(OO_flip>pi/2)-pi;
            OO_flip(OO_flip<-pi/2)=OO_flip(OO_flip<-pi/2)+pi;
            OO_flip(OO_flip>pi/2)=OO_flip(OO_flip>pi/2)-pi;
            
            orienation_map_filtered(sub2ind(size(currentImg), VIF_YY,VIF_XX))=OO_flip;
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
        
        enhanced_im_r = 255-currentImg;
        enhanced_im_g = 255-currentImg;
        enhanced_im_b = 255-currentImg;
        
        enhanced_im_r(find(current_seg>0))=255*R_seg_orient_heat_map(find(current_seg>0));
        enhanced_im_g(find(current_seg>0))=255*G_seg_orient_heat_map(find(current_seg>0));
        enhanced_im_b(find(current_seg>0))=255*B_seg_orient_heat_map(find(current_seg>0));
        
        RGB_seg_orient_heat_map(:,:,1 ) = enhanced_im_r;
        RGB_seg_orient_heat_map(:,:,2 ) = enhanced_im_g;
        RGB_seg_orient_heat_map(:,:,3 ) = enhanced_im_b;
        
        for sub_i = 1 : Sub_Sample_Num
            if iFrame + sub_i-1 <= nFrame
                imwrite(RGB_seg_orient_heat_map, ...
                    [HeatEnhOutputDir,'/white_segment_heat_',...
                    filename_short_strs{iFrame+ sub_i-1},'.tif']);
            end
        end
        
        RGB_seg_orient_heat_map_nms=[];
        
        if(strcmp(Combine_Way,'st_nms_two'))
            RGB_seg_orient_heat_map_nms = RGB_seg_orient_heat_map*0;
            
            enhanced_im_r = currentImg;
            enhanced_im_g = currentImg;
            enhanced_im_b = currentImg;
            
            enhanced_im_r(find(NMS_Segment>0))=255*R_seg_orient_heat_map(find(NMS_Segment>0));
            enhanced_im_g(find(NMS_Segment>0))=255*G_seg_orient_heat_map(find(NMS_Segment>0));
            enhanced_im_b(find(NMS_Segment>0))=255*B_seg_orient_heat_map(find(NMS_Segment>0));
            
            RGB_seg_orient_heat_map_nms(:,:,1 ) = enhanced_im_r;
            RGB_seg_orient_heat_map_nms(:,:,2 ) = enhanced_im_g;
            RGB_seg_orient_heat_map_nms(:,:,3 ) = enhanced_im_b;
            
            
            for sub_i = 1 : Sub_Sample_Num
                if iFrame + sub_i-1 <= nFrame
                    imwrite(RGB_seg_orient_heat_map_nms, ...
                        [HeatEnhOutputDir,'/NMS_Segment_heat_',...
                        filename_short_strs{iFrame+ sub_i-1},'.tif']);
                end
            end            
            RGB_seg_orient_heat_map = RGB_seg_orient_heat_map_nms;
        end
        
        current_seg_orientation = current_seg.*orienation_map_filtered;
        current_seg_orientation(find(current_seg==0)) = nan;
        end_points_map = bwmorph(current_seg,'endpoints');
        
        tip_orientation = double(end_points_map).*double(orienation_map_filtered);
        tip_int = double(end_points_map).*double(currentImg);
        tip_NMS = double(end_points_map).*double(nms);
        tip_orientation(find(end_points_map==0)) = nan;
        tip_int(find(end_points_map==0)) = nan;
        tip_NMS(find(end_points_map==0)) = nan;
                
        
        % see if canny test array exist
        if(~exist('current_seg_canny_cell','var'))
            current_seg_canny_cell=[];
        end
        
        %% Save segmentation results
        for sub_i = 1 : Sub_Sample_Num
            if iFrame + sub_i-1 <= nFrame
                save([DataOutputDir,'/steerable_vote_', ...
                    filename_short_strs{iFrame+ sub_i-1},'.mat'],...
                    'currentImg','orienation_map_filtered','OrientationVoted','orienation_map','RGB_seg_orient_heat_map','RGB_seg_orient_heat_map_nms', ...
                    'MAX_st_res', 'current_seg','Intensity_Segment','SteerabelRes_Segment','NMS_Segment', ...
                    'current_model', 'RGB_seg_orient_heat_map','current_seg_orientation','tip_orientation',...
                    'tip_int','tip_NMS',...
                    'current_seg_canny_cell');
                
            end
        end
        
        
        %% %tif stack cost too much memory, comment these
        %         if( save_tif_flag==1)
        % %             current_seg = (imread([FilamentSegmentationChannelOutputDir,'/segment_binary_',filename_short_strs{iFrame},'.tif']))>0;
        % %             RGB_seg_orient_heat_map = imread([HeatEnhOutputDir,'/segment_heat_',filename_short_strs{iFrame},'.tif']);
        % %
        %             tif_stack_binary_seg_image_data(:,:,iFrame_index) = uint8(current_seg*255);
        %             tif_stack_RGB_heat_image_data(:,:,:,iFrame_index) = uint8(RGB_seg_orient_heat_map);
        %
        %         end
        %%
    end
    %% For Gelfand Lab, save results as tif stack file
    if( save_tif_flag==1)
        options.comp = false;
        options.ask = false;
        options.message = true;
        options.append = false;
        
        % Save the multi-frame RGB color image
        options.color = true;
        %         saveastiff(tif_stack_RGB_heat_image_data, [FilamentSegmentationProcessOutputDir,'/channel_',num2str(iChannel),'_seg_heat.tif'], options);
        options.color = false;
        %         saveastiff(tif_stack_binary_seg_image_data, [FilamentSegmentationProcessOutputDir,'/channel_',num2str(iChannel),'_seg_binary.tif'], options);
        
    end
    
end



%% For Gelfand Lab, outgrowth calculation
if(VIF_Outgrowth_Flag==1)
    VIF_outgrowth_measurement(movieData);
end




