
% for now MT Channe 1, VIF/Actin/Adhesion... Channel 2

if(~exist('indexMTChannel','var'))
    indexMTChannel = 1;
end

if(~exist('indexVIFChannel','var'))
    indexVIFChannel = min(2,length(MD.channels_));
end

% display_msg_flag = 0;

% Find the package of Filament Analysis
nFrame = MD.nFrames_;

nPackage = length(MD.packages_);

indexTrackingPackage = 0;
for i = 1 : nPackage
    if(strcmp(MD.packages_{i}.getName,'Tracking')==1)
        indexTrackingPackage = i;
        break;
    end
end

% if(indexTrackingPackage==0)
%     disp('Need to run tracking package first.')
% %     return;
% end

indexFilamentPackage = 0;
for i = 1 : nPackage
    if(strcmp(MD.packages_{i}.getName,'FilamentAnalysis')==1)
        indexFilamentPackage = i;
        break;
    end
end

if(indexFilamentPackage==0)
    if(display_msg_flag>0)
        disp('No filament analysis package.')
    end
    %     return;
end

% Find the process of plusTiptracking.
nProcesses = length(MD.processes_);

indexMTTrackingProcess = 0;
for i = 1 : nProcesses
    if(strcmp(MD.processes_{i}.getName,'Tracking')==1)
        indexMTTrackingProcess = i;
        break;
    end
end

% if indexMTTrackingProcess==0
%     disp('Please run plusTipTracker package first.')
% %     return;
% end


indexMTPostTrackingProcess = 0;
for i = 1 : nProcesses
    if(strcmp(MD.processes_{i}.getName,'Microtubule dynamics classification')==1)
        indexMTPostTrackingProcess = i;
        break;
    end
end

% if indexMTPostTrackingProcess==0
%     disp('Please run plusTipTracker postprocessing first.')
% %     return;
% end

%%
indexSteerabeleProcess = 0;
for i = 1 : nProcesses
    if(strcmp(MD.processes_{i}.getName,'Steerable filtering')==1)
        indexSteerabeleProcess = i;
        break;
    end
end

if indexSteerabeleProcess==0
    if(display_msg_flag>0)
        
        disp('No steerable filtering process.');
    end
    %     return;
end

indexFilamentSegmentationProcess = 0;
for i = 1 : nProcesses
    if(strcmp(MD.processes_{i}.getName,'Filament Segmentation')==1)
        indexFilamentSegmentationProcess = i;
        break;
    end
end

if indexFilamentSegmentationProcess==0
    if(display_msg_flag>0)
        
        disp('Filament Segmentation parameters not set.')
    end
    %     return;
end

indexFlattenProcess = 0;
for i = 1 : nProcesses
    if(strcmp(MD.processes_{i}.getName,'Image Flatten')==1)
        indexFlattenProcess = i;
        break;
    end
end

if indexFlattenProcess==0
    if(display_msg_flag>0)        
        display('The setting shows you want to use flattened image for steerable filtering. Please set parameters for Image Flatten and run.')
    end
    %     return;
end

indexCellSegProcess = 0;
for i = 1 : nProcesses
    if(strcmp(MD.processes_{i}.getName,'Thresholding')==1)
        indexCellSegProcess = i;
        break;
    end
end

if indexCellSegProcess == 0
    if(display_msg_flag>0)        
        disp('Please run segmentation and refinement first.')
    end
    %     return;
end

indexCellRefineProcess = 0;
for i = 1 : nProcesses
    if(strcmp(MD.processes_{i}.getName,'Mask Refinement')==1)
        indexCellRefineProcess = i;
        break;
    end
end

if indexCellRefineProcess == 0
    if(display_msg_flag>0)        
        disp('Please run segmentation and refinement first.')
    end
    %     return;
end

% Get the sub sample rate on the VIF channel
if(indexSteerabeleProcess>0)
    Sub_Sample_Num_st = MD.processes_{indexSteerabeleProcess}.funParams_.Sub_Sample_Num;    
else    
    Sub_Sample_Num_st=1;    
end

if(indexFilamentSegmentationProcess>0)
    Sub_Sample_Num_filament = MD.processes_{indexFilamentSegmentationProcess}.funParams_.Sub_Sample_Num;
else
    Sub_Sample_Num_filament=1;
end


% In case they are not set to the same, take the largest common divider, if not 1
Sub_Sample_Num = gcd(Sub_Sample_Num_st,Sub_Sample_Num_filament);


% Get frame number from the title of the image, this not neccesarily
% the same as iFrame due to some shorting problem of the channel
if(isempty(MD.channels_(indexVIFChannel).getImageFileNames))
    tic
    for iFrame = 1 : nFrame
        filename_one = (MD.channels_(indexVIFChannel).getImageFileNames(iFrame));
        channel_filename{1,iFrame} = filename_one{1};        
    end
    toc
    filename_short_strs = uncommon_str_takeout(channel_filename);
else
    filename_short_strs = uncommon_str_takeout(MD.channels_(indexVIFChannel).getImageFileNames);
end

Frames_to_Seg = 1 : Sub_Sample_Num_filament : nFrame;
Frames_results_correspondence = im2col(repmat(Frames_to_Seg, [Sub_Sample_Num_filament,1]),[1 1]);
Frames_results_correspondence = Frames_results_correspondence(1:nFrame);

if indexFilamentPackage > 0 && indexFilamentSegmentationProcess > 0
    % Set the directories
    FilamentSegmentationProcessOutputDir  = [MD.packages_{indexFilamentPackage}.outputDirectory_, filesep 'FilamentSegmentation'];
    FilamentSegmentationChannelOutputDir =  MD.processes_{indexFilamentSegmentationProcess}.outFilePaths_{indexVIFChannel};
    HeatOutputDir = [FilamentSegmentationChannelOutputDir,'/HeatOutput'];
    HeatEnhOutputDir = [HeatOutputDir,'/Enh'];
    DataOutputDir = [FilamentSegmentationChannelOutputDir,'/DataOutput'];
end
