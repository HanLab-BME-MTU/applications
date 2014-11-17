function  F_classifer  = load_classifier_trained(movieData)
% nms_classifier_train trains an classifier for segments filaments from input image(nms) based on the geometrical features of the curves/lines in the image and user input
% Input:
%    MD:                                the MD for the image project
% Output:
%    F_classifer:                       trained classifier
%

% Liya Ding
% 2013.03



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


nFrame = movieData.nFrames_;

for iChannel = selected_channels
    Channel_FilesNames = movieData.channels_(iChannel).getImageFileNames(1:movieData.nFrames_);
    
    % Make output directory for the steerable filtered images
    FilamentSegmentationChannelOutputDir =  movieData.processes_{indexFilamentSegmentationProcess}.outFilePaths_{iChannel};
    if (~exist(FilamentSegmentationChannelOutputDir,'dir'))
        mkdir(FilamentSegmentationChannelOutputDir);
    end
    
    load([FilamentSegmentationChannelOutputDir,'/F_classifer_channel.mat'],'F_classifer_train_this_channel');
    F_classifer{iChannel} = F_classifer_train_this_channel;
end



