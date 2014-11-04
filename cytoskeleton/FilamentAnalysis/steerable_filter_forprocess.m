function movieData = steerable_filter_forprocess(movieData, paramsIn, varargin)
% Created 07 2012 by Liya Ding, Matlab R2011b

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

% Find the process of segmentation mask refinement.
nProcesses = length(movieData.processes_);

indexSteerabeleProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Steerable filtering')==1)
        indexSteerabeleProcess = i;
        break;
    end
end

if indexSteerabeleProcess==0
    msgbox('Please set parameters for steerable filtering.')
    return;
end


% with no input funparam, use the one the process has on its own
if nargin < 2
    paramsIn = [];
    funParams = movieData.processes_{indexSteerabeleProcess}.funParams_;
else
    funParams = paramsIn;    
end

selected_channels = funParams.ChannelIndex;
BaseSteerableFilterSigma = funParams.BaseSteerableFilterSigma;
Levelsofsteerablefilters = funParams.Levelsofsteerablefilters;
ImageFlattenFlag = funParams.ImageFlattenFlag;
Sub_Sample_Num = funParams.Sub_Sample_Num;

%% Output Directories

% default steerable filter process output dir
ImageSteerableFilterProcessOutputDir = [movieData.outputDirectory_, filesep 'SteerableFiltering'];

% if there is filamentanalysispackage
if (indexFilamentPackage>0)
    % and a directory is defined for this package
    if (~isempty(movieData.packages_{indexFilamentPackage}.outputDirectory_))
        % and this directory exists
        if (exist(movieData.packages_{indexFilamentPackage}.outputDirectory_,'dir'))
            ImageSteerableFilterProcessOutputDir  = [movieData.packages_{indexFilamentPackage}.outputDirectory_, filesep 'SteerableFiltering'];
        end
    end
end

if (~exist(ImageSteerableFilterProcessOutputDir,'dir'))
    mkdir(ImageSteerableFilterProcessOutputDir);
end


for iChannel = selected_channels
    ImageSteerableFilterChannelOutputDir = [ImageSteerableFilterProcessOutputDir,'/Channel',num2str(iChannel)];
    if (~exist(ImageSteerableFilterChannelOutputDir,'dir'))
        mkdir(ImageSteerableFilterChannelOutputDir);
    end
    
    
    output_dir_content = dir(fullfile([ImageSteerableFilterChannelOutputDir,filesep,'*.*']));
    
    %if there are files in this dir, clear them
    if(length(output_dir_content)>2)
        delete([ImageSteerableFilterChannelOutputDir,filesep,'*.*']);
        if(exist([ImageSteerableFilterChannelOutputDir,filesep,'NMS'],'dir'))
            rmdir([ImageSteerableFilterChannelOutputDir,filesep,'NMS'], 's');
        end
    end
    
    movieData.processes_{indexSteerabeleProcess}.setOutImagePath(iChannel,ImageSteerableFilterChannelOutputDir);
end

%%
indexFlattenProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Image Flatten')==1)
        indexFlattenProcess = i;
        break;
    end
end


if indexFlattenProcess==0 && ImageFlattenFlag == 2
    display('The setting shows you want to use flattened image for steerable filtering. Please set parameters for Image Flatten and run.')
    return;
end

nFrame = movieData.nFrames_;

% RES_cell = cell(1,nFrame);
% Image_cell = cell(1,nFrame);

Frames_to_Seg = 1:Sub_Sample_Num:nFrame;
Frames_results_correspondence = im2col(repmat(Frames_to_Seg, [Sub_Sample_Num,1]),[1 1]);
Frames_results_correspondence = Frames_results_correspondence(1:nFrame);

for iChannel = selected_channels
    % Get frame number from the title of the image, this not neccesarily
    % the same as iFrame due to some shorting problem of the channel
    Channel_FilesNames = movieData.channels_(iChannel).getImageFileNames(1:movieData.nFrames_);
    
    filename_short_strs = uncommon_str_takeout(Channel_FilesNames);
    
    % Make output directory for the steerable filtered images
    ImageSteerableFilterChannelOutputDir = movieData.processes_{indexSteerabeleProcess}.outFilePaths_{iChannel};
    
    if (~exist(ImageSteerableFilterChannelOutputDir,'dir'))
        mkdir(ImageSteerableFilterChannelOutputDir);
    end
    
    
    display('======================================');
    
    display(['Current movie: as in ',movieData.outputDirectory_]);
    
    display(['Start steerable filtering in Channel ',num2str(iChannel)]);
    
    if (~exist([ImageSteerableFilterChannelOutputDir,filesep,'NMS'],'dir'))
        mkdir(ImageSteerableFilterChannelOutputDir,'NMS');
    end
    
    for iFrame_subsample = 1 : length(Frames_to_Seg)
        iFrame = Frames_to_Seg(iFrame_subsample);
        disp(['Frame: ',num2str(iFrame)]);
        
        % Read in the intensity image.
        if indexFlattenProcess > 0
            currentImg = imread([movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel}, filesep, 'flatten_',filename_short_strs{iFrame},'.tif']);
        else
            currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        end
        currentImg = double(currentImg);
        
        levels_sizes = 2.^((1:Levelsofsteerablefilters)-1);
        
        % Steerable filtering using four scales one doubling the previous one.
        % function multiscaleSteerableDetector will automatically merge the results
%         [MAX_st_res, orienation_map, nms, scaleMap] = multiscaleSteerableDetector(currentImg, 4, BaseSteerableFilterSigma.*levels_sizes);
        
        [MAX_st_res, orienation_map, nms, scaleMap] = multiscaleSteerableDetector(currentImg, 4, BaseSteerableFilterSigma.*levels_sizes);
        
        for sub_i = 1 : Sub_Sample_Num
            if iFrame + sub_i-1 <= nFrame
                imwrite((MAX_st_res)/(max(max(MAX_st_res))), ...
                    [ImageSteerableFilterChannelOutputDir,'/MAX_st_res_', ...
                    filename_short_strs{iFrame + sub_i-1},'.tif']);
               imwrite((nms)/(max(max(nms))), ...
                    [ImageSteerableFilterChannelOutputDir,'/NMS/NMS_', ...
                    filename_short_strs{iFrame + sub_i-1},'.tif']);
            end
        end
        
        save([ImageSteerableFilterChannelOutputDir,'/steerable_',filename_short_strs{iFrame},'.mat'],...
            'orienation_map', 'MAX_st_res','nms','scaleMap');
    end
end
