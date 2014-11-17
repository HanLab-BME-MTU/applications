function Whole_movie_stat_cell = whole_movie_stat_function(movieData)

%%
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

% Output Directories

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


%%
funParams.Whole_movie_stat_cell = cell(1,max(funParams.ChannelIndex));

Whole_movie_stat.mean_ST = 0;
Whole_movie_stat.std_ST = 0;
Whole_movie_stat.mode_ST = 0;
Whole_movie_stat.otsu_ST = 0;
Whole_movie_stat.otsu_mode_ST = 0;

Whole_movie_stat.mean_NMS = 0;
Whole_movie_stat.std_NMS = 0;
Whole_movie_stat.mode_NMS = 0;
Whole_movie_stat.otsu_NMS = 0;
Whole_movie_stat.otsu_mode_NMS = 0;

Whole_movie_stat.mean_INT = 0;
Whole_movie_stat.std_INT = 0;
Whole_movie_stat.mode_INT = 0;
Whole_movie_stat.otsu_INT = 0;
Whole_movie_stat.otsu_mode_INT = 0;

for iCh = 1 :max(funParams.ChannelIndex)
    funParams.Whole_movie_stat_cell{iCh} = Whole_movie_stat;
end


%%
nFrame = movieData.nFrames_;

for iChannel = selected_channels
    
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
    
    Frames_to_Seg = 1:Sub_Sample_Num:nFrame;
    Frames_results_correspondence = im2col(repmat(Frames_to_Seg, [Sub_Sample_Num,1]),[1 1]);
    Frames_results_correspondence = Frames_results_correspondence(1:nFrame);
    
    INT_pool = [];
    
    whole_modie_selectframe_pace = max(1,round(length(Frames_to_Seg)/30));
    
    for iFrame_index = 1 : whole_modie_selectframe_pace:length(Frames_to_Seg)
        iFrame = Frames_to_Seg(iFrame_index);
        
        % Read in the intensity image.
        if indexFlattenProcess > 0 && ImageFlattenFlag==2
            currentImg = imread([movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel}, filesep, 'flatten_',filename_short_strs{iFrame},'.tif']);
        else
            currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        end
        
        currentImg = currentImg(1:3:end,1:3:end);
        
        INT_pool = [INT_pool; currentImg(:)];
    end
    
    INT_pool =  double(INT_pool);
    
    % find the mode of the intensity of the curves/lines
    [hist_n,bin] = hist((INT_pool),50);
    ind_mode = find(hist_n==max(hist_n));
    mode_INT = bin(ind_mode(1));
    
    Whole_movie_stat.mean_INT = mean(INT_pool);
    Whole_movie_stat.std_INT = std(INT_pool);
    Whole_movie_stat.mode_INT = mode_INT;
    Whole_movie_stat.otsu_INT = thresholdOtsu(INT_pool);
    Whole_movie_stat.otsu_mode_INT = thresholdOtsu(INT_pool(find(INT_pool>mode_INT)));
    
    INT_pool = [];
    
    
    ST_pool = [];
    % this line in commandation for shortest version of filename
    filename_shortshort_strs = all_uncommon_str_takeout(Channel_FilesNames{1});
    
    for iFrame_index = 1 : whole_modie_selectframe_pace: length(Frames_to_Seg)
        iFrame = Frames_to_Seg(iFrame_index);
        
        
        try
            load([SteerableChannelOutputDir, filesep, 'steerable_',...
                filename_short_strs{iFrame},'.mat']);
        catch
            % in the case of only having the short-old version
            load([SteerableChannelOutputDir, filesep, 'steerable_',...
                filename_shortshort_strs{iFrame},'.mat']);
        end
        
        
        MAX_st_res = MAX_st_res(1:3:end,1:3:end);
        MAX_st_res = MAX_st_res(MAX_st_res>0);
        
        ST_pool = [ST_pool; MAX_st_res(:)];
        
    end
    
    [hist_n,bin] = hist(ST_pool,50);
    ind_mode = find(hist_n==max(hist_n));
    mode_ST = bin(ind_mode(1));
    
    Whole_movie_stat.mean_ST = mean(ST_pool);
    Whole_movie_stat.std_ST = std(ST_pool);
    Whole_movie_stat.mode_ST = mode_ST;
    Whole_movie_stat.otsu_ST = thresholdOtsu(ST_pool);
    Whole_movie_stat.otsu_mode_ST = thresholdOtsu(ST_pool(find(ST_pool>mode_ST)));
    
    try
        Whole_movie_stat.rosin_ST = thresholdRosin(ST_pool);
    catch
        Whole_movie_stat.rosin_ST = Whole_movie_stat.otsu_ST;
    end
    
    try
        Whole_movie_stat.rosin_mode_ST = thresholdRosin(ST_pool(find(ST_pool>mode_ST)));
    catch
        Whole_movie_stat.rosin_mode_ST = Whole_movie_stat.otsu_mode_ST;
    end
    
    
    ST_pool = [];
    
    NMS_pool = [];
    
    Length_pool=[];
    
    for iFrame_index = 1 : whole_modie_selectframe_pace:length(Frames_to_Seg)
        iFrame = Frames_to_Seg(iFrame_index);
                
        try
            load([SteerableChannelOutputDir, filesep, 'steerable_',...
                filename_short_strs{iFrame},'.mat']);
        catch
            % in the case of only having the short-old version
            load([SteerableChannelOutputDir, filesep, 'steerable_',...
                filename_shortshort_strs{iFrame},'.mat']);
        end
        
        imageNMS= nms;
        
        [hist_n,bin] = hist(nms(find(nms>0)),200);
        ind_mode = find(hist_n==max(hist_n));
        mode_nms = bin(ind_mode(1));
        % And find the Otsu threshold for the intensity
        T_otsu = thresholdOtsu(nms(find(nms>mode_nms)));
        T_otsu_start =  max(mode_nms/3, (-abs(T_otsu - mode_nms)*0.1+mode_nms));
        if(isempty(T_otsu_start))
           T_otsu_start = mode_nms;
        end
        % first, get almost all the curves/lines, by using a low threshold
        imageMask = imageNMS > T_otsu_start;
        
        % further thin it, since the nms version of steerable filtering is not real skeleton
        bw_out = bwmorph(imageMask,'thin','inf');
        
        % Find the branching points
        nms_seg_brancing = bwmorph(bw_out,'branchpoints');
        
        % Delete these branching points for now
        nms_seg_no_brancing = bw_out - nms_seg_brancing;
        
        % again Find the branching points
        nms_seg_brancing = bwmorph(nms_seg_no_brancing,'branchpoints');
        
        % Delete these branching points for now
        nms_seg_no_brancing = nms_seg_no_brancing - nms_seg_brancing;
        
        % Label all isolated lines(curves)
        labelMask = bwlabel(nms_seg_no_brancing);
        
        % Get properties for each of curve
        ob_prop = regionprops(labelMask,'Area','MajorAxisLength','Eccentricity','MinorAxisLength','Centroid');
        
        obAreas = [ob_prop.Area];
        
        nms = nms(1:3:end,1:3:end);
        
        nms = nms(~isnan(nms));
        nms = nms(nms>0);
        
        NMS_pool = [NMS_pool; nms(:)];
        
        Length_pool = [Length_pool obAreas];
        
    end
        
    [hist_n,bin] = hist(NMS_pool,50);
    ind_mode = find(hist_n==max(hist_n));
    mode_NMS = bin(ind_mode(1));
    
    Whole_movie_stat.mean_NMS = mean(NMS_pool);
    Whole_movie_stat.std_NMS = std(NMS_pool);
    Whole_movie_stat.mode_NMS = mode_NMS;
    Whole_movie_stat.otsu_NMS = thresholdOtsu(NMS_pool);
    Whole_movie_stat.otsu_mode_NMS = thresholdOtsu(NMS_pool(find(NMS_pool>mode_NMS)));
    
    try
      Whole_movie_stat.rosin_NMS = thresholdRosin(NMS_pool);
    catch
        Whole_movie_stat.rosin_NMS = Whole_movie_stat.otsu_NMS;
    end
    
    try
    Whole_movie_stat.rosin_mode_NMS = thresholdRosin(NMS_pool(find(NMS_pool>mode_NMS)));
    catch
        Whole_movie_stat.rosin_mode_NMS = Whole_movie_stat.otsu_mode_NMS;
    end
    
    Length_pool = Length_pool(Length_pool>1);
    
    
    [hist_n,bin] = hist(Length_pool,50);
    ind_mode = find(hist_n==max(hist_n));
    mode_Length = bin(ind_mode(1));
    
    Whole_movie_stat.mean_Length = mean(Length_pool);
    Whole_movie_stat.std_Length = std(Length_pool);
    Whole_movie_stat.mode_Length = mode_Length;
    Whole_movie_stat.rosin_Length = thresholdRosin(Length_pool);
    Whole_movie_stat.rosin_mode_Length = thresholdRosin(Length_pool(find(Length_pool>mode_Length)));
    
    
    NMS_pool = [];
    
    Whole_movie_stat_cell{iChannel} = Whole_movie_stat;
    
end

