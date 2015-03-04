function movieData = filament_segmentation_wholemoviestat_save(movieData, paramsIn, wholemovie_input_filename, varargin)
% Function to take care of the whole movie stat before the filament segmentation part is done

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

% Created 2015.2. by Liya Ding, Matlab R2012b

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
saveallresults_movie = funParams.savestepfigures;

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
    FilamentSegmentationChannelOutputDir = [FilamentSegmentationProcessOutputDir,filesep,'Channel',num2str(iChannel)];
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

if indexCellSegProcess == 0 && (Cell_Mask_ind(1) == 1 || Cell_Mask_ind(1) == 3 || Cell_Mask_ind(1) == 4 || Cell_Mask_ind(1) == 6)
    msgbox('Please run segmentation and refinement first.')
    return;
end
Rerun_WholeMovie = 1;
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
        
        % check if the existing whole movie file include the currently
        % selected channel
        flag_complete = zeros(numel(movieData.channels_),1);
        for iChannel = selected_channels
            if(numel(Whole_movie_stat_cell)<iChannel)
                flag_complete(iChannel)=0;
            else
                if(isempty(Whole_movie_stat_cell{iChannel}))
                    flag_complete(iChannel)=0;
                end
            end
        end
        
        %if some channels are missing, rerun the whole_movie_stat
        if(min(flag_complete)==0)        
            % this version of "addon" whole_movie_stat_function
            % accept what ever was in the mat file
            % and do for the missing channel(this previous channels will be
            % kept even if it is not selected in current setting
            Whole_movie_stat_cell = whole_movie_stat_function_addon(movieData);
            save([FilamentSegmentationProcessOutputDir, filesep, 'whole_movie_stat.mat'],'Whole_movie_stat_cell');
        end
    else
        % this version of whole_movie_stat_function calculate for
        % currently selected channels, disregarding whether any thing
        % already existed
        Whole_movie_stat_cell = whole_movie_stat_function(movieData);
        save([FilamentSegmentationProcessOutputDir, filesep, 'whole_movie_stat.mat'],'Whole_movie_stat_cell');
    end
    
    funParams.Whole_movie_stat_cell = Whole_movie_stat_cell;
    
end
