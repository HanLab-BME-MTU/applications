function  F_classifer_train_output  = nms_classifier_train_moviedata(movieData,channelIndex)
% nms_classifier_train trains an classifier for segments filaments from input image(nms) based on the geometrical features of the curves/lines in the image and user input
% Input:
%    MD:                                the MD for the image project
%    channelIndex:                      the channels to run through ( at this point the process might be not ready yet)
% Output:
%    F_classifer:                       trained classifiers for selected channels, one for each
%

% Liya Ding
% 2013.03


save_tif_flag=1;
F_classifer_train_output= cell(1,1);

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


% if there is no st results, but the combine way is not int only, then st
% is missing, error here.
if indexSteerabeleProcess==0 && (~strcmp(Combine_Way,'int_only'))
    display('Please run steerable filtering first.')
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

% if indexCellSegProcess == 0 && Cell_Mask_ind == 1
%     msgbox('Please run segmentation and refinement first.')
%     return;
% end
%

nFrame = movieData.nFrames_;


F_classifer_train_output = cell(1,max(channelIndex));

for iChannel = channelIndex
    
    Channel_FilesNames = movieData.channels_(iChannel).getImageFileNames(1:movieData.nFrames_);
    
    filename_short_strs = uncommon_str_takeout(Channel_FilesNames);
    
    % Make output directory for the steerable filtered images
    FilamentSegmentationChannelOutputDir =  [FilamentSegmentationProcessOutputDir, '/Channel',num2str(iChannel)];
    if (~exist(FilamentSegmentationChannelOutputDir,'dir'))
        mkdir(FilamentSegmentationChannelOutputDir);
    end
    
    
    % If steerable filter process is run
    if indexSteerabeleProcess>0
        SteerableChannelOutputDir = movieData.processes_{indexSteerabeleProcess}.outFilePaths_{iChannel};
    end
    
    iFrame_index = 1;
    iFrame = 1;
    
    % Read in the intensity image.
    if indexFlattenProcess > 0 && ImageFlattenFlag==2
        currentImg = imread([movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel}, filesep, 'flatten_',filename_short_strs{iFrame},'.tif']);
    else
        currentImg = movieData.channels_(iChannel).loadImage(iFrame);
    end
    currentImg = double(currentImg);
    
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
    
    imageNMS = nms;
    
    
    imageInt = currentImg;
    
    T_Rosin_otsu = thresholdRosin(imfilter(imageInt,fspecial('gaussian',11,2)));
    MaskCell = imageInt>T_Rosin_otsu*6/3;
    MaskCell = imdilate(MaskCell,fspecial('disk', 71)>0);
    
    [hist_n,bin] = hist(imageNMS(find(imageNMS>0)),200);
    ind_mode = find(hist_n==max(hist_n));
    mode_nms = bin(ind_mode(1));
    % And find the Otsu threshold for the intensity
    T_otsu = thresholdOtsu(imageNMS(find(imageNMS>mode_nms)));
    T_otsu_start =  abs(T_otsu - mode_nms)*0.2 + mode_nms;
    
    imageNMS = imageNMS.*MaskCell;
    imageInt = imageInt.*MaskCell;
    
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
    %     ob_prop = regionprops(labelMask,'Area','MinorAxisLength','MajorAxisLength','Centroid');
    ob_prop = regionprops(labelMask,'Area','Centroid');
    
    % Redefine variable for easy of notation
    obAreas = [ob_prop.Area];
    
    nLine = length(obAreas);
    
    % Some feature for later consideration
    %     obLongaxis = [ob_prop.MajorAxisLength];
    %     obShortaxis = [ob_prop.MinorAxisLength];
    obCentroid = zeros(2, length(obAreas));
    obCentroid(:) = [ob_prop.Centroid];
    % The ratio of short vs long axis
    %     ratio  = obShortaxis./obLongaxis;
    
    
    feature_MeanInt = nan(nLine,1);
    feature_MeanNMS = nan(nLine,1);
    feature_MeanX = zeros(nLine,1);
    feature_MeanY = zeros(nLine,1);
    feature_Length = obAreas';
    feature_Curvature = nan(nLine,1);
    
    feature_MeanX = (obCentroid(1,:))';
    feature_MeanY = (obCentroid(2,:))';
    
    % for the features, only include those curves/lines longer than 4 pixels
    ind_long = find(feature_Length>4);
    
    ordered_points = cell(1,1);
    smoothed_ordered_points = cell(1,1);
    
    tic
    % get the mean intensity of the curves
    for i_area = ind_long'
        [all_y_i, all_x_i] = find(labelMask == i_area);
        NMS = imageNMS(sub2ind(size(bw_out), round(all_y_i),round(all_x_i)));
        feature_MeanNMS(i_area) = mean(NMS);
        INT = imageInt(sub2ind(size(bw_out), round(all_y_i),round(all_x_i)));
        feature_MeanInt(i_area) = mean(INT);
        
        feature_MeanX(i_area) = mean(all_x_i);
        feature_MeanY(i_area) = mean(all_y_i);
        
    end
    display('Time spend in look for mean int:');
    toc
    % figure; plot3(feature_Length,feature_MeanInt,feature_Curvature,'.');
    
    % find the mode of the intensity of the curves/lines
    [hist_n,bin] = hist(feature_MeanNMS,200);
    ind_mode = find(hist_n==max(hist_n));
    mode_nms = bin(ind_mode(1));
    % And find the Otsu threshold for the intensity
    hotsu = thresholdOtsu(feature_MeanNMS(find(feature_MeanNMS>mode_nms)));
    
    % the up one
    % Set the slanted classification line cutoff as twice of the Otsu with
    % respect to the mode
    T_xie_int_up =  abs(hotsu - mode_nms)*2 + mode_nms;
    
    % And the length as Otsu threshold
    T_xie_length_up = 2*max(thresholdOtsu(feature_Length),thresholdRosin(feature_Length));
    
    % Make a classification function as whether it is above the line
    F_classifer_up = @(nms,length,curvature,int) (((T_xie_int_up + (T_xie_int_up/T_xie_length_up)*(-length) )<nms));
    
    % the down one
    T_xie_int_down =  abs(hotsu - mode_nms)*1.5 + mode_nms;
    T_xie_length_down = 1.5*max(thresholdOtsu(feature_Length),thresholdRosin(feature_Length));
    
    F_classifer_down = @(nms,length,curvature,int)  (((T_xie_int_down + (T_xie_int_down/T_xie_length_down)*(-length) )<nms));
    
    % the points in between the two classifier is the not sure ones that
    % requires annotation
    not_sure_ind = find(F_classifer_up(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)==0 & F_classifer_down(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)>0);
    good_ind = find(F_classifer_up(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)>0);
    
    F_classifer_notuse =  @(nms,length,curvature,int)  (((T_xie_int_down/2 + (T_xie_int_down/T_xie_length_down)*(-length) )<nms));
    bad_ind = find(F_classifer_down(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)==0 & F_classifer_notuse(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)>0);
    
    
    h1 = figure(1); hold off;
    imagescc(currentImg); hold on;
    scrsz = get(0,'ScreenSize');
    set(h1,'Position',scrsz);
    training_sample_number=funParams.training_sample_number;
    
    MaskCell = keep_largest_area(MaskCell);
    S_mask = regionprops(MaskCell,'MajorAxisLength');
    
    radius_for_kd = (S_mask.MajorAxisLength)/20;
    
    tic
    [idx, dist] = KDTreeBallQuery([feature_MeanX feature_MeanY], [feature_MeanX feature_MeanY], radius_for_kd);
    display('Find Neighbor');
    toc
    
    tic
    for fi = 1 : length(feature_MeanX)
        neighbor_number(fi) = length(idx{fi});
    end
    display('Set number of Neighbor');
    toc
    
    W_neighbor = 1./neighbor_number;
    training_ind = randsample(not_sure_ind,min(length(not_sure_ind),training_sample_number),true,W_neighbor(not_sure_ind));
    training_ind_unique = unique(training_ind);
    
    if(length(training_ind_unique)<training_sample_number)
        not_sure_ind_again = setdiff(not_sure_ind,training_ind_unique);
        training_ind_again = randsample(not_sure_ind_again,min(length(not_sure_ind_again),training_sample_number-length(training_ind_unique)),true,W_neighbor(not_sure_ind_again));
        training_ind = [training_ind_unique;training_ind_again]';
    end
    
    
    display('Random Sample');
    
    training_bad_ind = [];
    training_good_ind = [];
    
    good_bad_label = [];
    i_ind = 1;
    
    for i_ind = 1 : length(training_ind)
        i_ind
        iLine = training_ind(i_ind);
        [y,x] = find(labelMask==iLine);
        h1=figure(1);hold on;
        plot(x,y,'y.');
    end
    
    i_ind = 1;
    
    for i_mark = 1 : 2*length(training_ind)
        i_ind
        iLine = training_ind(i_ind);
        [y,x] = find(labelMask==iLine);
        h1=figure(1);title(['Frame: ',num2str(i_mark)]);
        plot(x,y,'g.');
        saveas(h1,[FilamentSegmentationChannelOutputDir,'/train_g_',num2str(i_mark),'.jpg']);
        
        AA = [max(1, round(mean(x))-50), min(size(imageNMS,2), round(mean(x))+50), ...
            max(1, round(mean(y))-50), min(size(imageNMS,1), round(mean(y))+50)];
        axis(AA);
        ch = getkey();
        
        if(strcmp(ch,'')==1)
            ch=0;
        end
        
        if(ch==28)
            i_ind = max(1, i_ind - 1);
        else
            if(ch==32)
                good_bad_label(i_ind) = 1;
                plot(x,y,'r.');
                i_ind = i_ind + 1;
            else
                good_bad_label(i_ind) = 2;
                plot(x,y,'b.');
                i_ind = i_ind + 1;
            end
            saveas(h1,[FilamentSegmentationChannelOutputDir,'/train_rb_',num2str(i_mark),'.jpg']);
            
            
            if(ch==27 || i_ind > length(training_ind))
                good_bad_label = good_bad_label(1:length(training_ind));
                break;
            end
        end
    end
    
    % if everything is good, then lower the cut
    if(mean(good_bad_label)>1.9)
        
        % the up one
        % Set the slanted classification line cutoff as twice of the Otsu with
        % respect to the mode
        T_xie_int_up =  abs(hotsu - mode_nms)*1.5 + mode_nms;
        
        % And the length as Otsu threshold
        T_xie_length_up = 1.5*max(thresholdOtsu(feature_Length),thresholdRosin(feature_Length));
        
        % Make a classification function as whether it is above the line
        F_classifer_up = @(nms,length,curvature,int)  (((T_xie_int_up + (T_xie_int_up/T_xie_length_up)*(-length) )<nms));
        
        
        % the down one
        T_xie_int_down =  abs(hotsu - mode_nms)*1 + mode_nms;
        T_xie_length_down = 1.2*max(thresholdOtsu(feature_Length),thresholdRosin(feature_Length));
        
        F_classifer_down = @(nms,length,curvature,int)  (((T_xie_int_down + (T_xie_int_down/T_xie_length_down)*(-length) )<nms));
        
        % the points in between the two classifier is the not sure ones that
        % requires annotation
        not_sure_ind = find(F_classifer_up(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)==0 & F_classifer_down(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)>0);
        good_ind = find(F_classifer_up(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)>0);
        
        F_classifer_notuse =  @(nms,length,curvature,int)  (((T_xie_int_down/2 + (T_xie_int_down/T_xie_length_down)*(-length) )<nms));
        bad_ind = find(F_classifer_down(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)==0 ...
            & F_classifer_notuse(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)>0);
        
        
        h1 = figure(1); hold off;
        imagescc(currentImg); hold on;
        scrsz = get(0,'ScreenSize');
        set(h1,'Position',scrsz);
        
        
        training_ind = randsample(not_sure_ind,min(length(not_sure_ind),training_sample_number),true,W_neighbor(not_sure_ind));
        training_ind_unique = unique(training_ind);
        
        if(length(training_ind_unique)<training_sample_number)
            not_sure_ind_again = setdiff(not_sure_ind,training_ind_unique);
            training_ind_again = randsample(not_sure_ind_again,min(length(not_sure_ind_again),training_sample_number-length(training_ind_unique)),true,W_neighbor(not_sure_ind_again));
            training_ind = [training_ind_unique;training_ind_again]';
        end
        
        training_bad_ind = [];
        training_good_ind = [];
        
        good_bad_label = [];
        i_ind = 1;
        
        for i_mark = 1 : 2*length(training_ind)
            i_ind
            iLine = training_ind(i_ind);
            [y,x] = find(labelMask==iLine);
            h1=figure(1);title(['Frame: ',num2str(i_mark)]);
            plot(x,y,'g.');
            saveas(h1,[FilamentSegmentationChannelOutputDir,'/train_g_',num2str(i_mark),'.jpg']);
            
            AA = [max(1, round(mean(x))-50), min(size(imageNMS,2), round(mean(x))+50), ...
                max(1, round(mean(y))-50), min(size(imageNMS,1), round(mean(y))+50)];
            axis(AA);
            ch = getkey();
            if(ch==28)
                i_ind = max(1, i_ind - 1);
            else
                if(ch==32)
                    good_bad_label(i_ind) = 1;
                    plot(x,y,'r.');
                    i_ind = i_ind + 1;
                else
                    good_bad_label(i_ind) = 2;
                    plot(x,y,'b.');
                    i_ind = i_ind + 1;
                end
                saveas(h1,[FilamentSegmentationChannelOutputDir,'/train_rb_',num2str(i_mark),'.jpg']);
                
                if(ch==27 || i_ind > length(training_ind))
                    good_bad_label = good_bad_label(1:length(training_ind));
                    break;
                end
                
            end
        end
        
        
    end
    
    
    % if everything is bad, then higher the cut
    if(mean(good_bad_label)<1.1)
        
        % the up one
        % Set the slanted classification line cutoff as twice of the Otsu with
        % respect to the mode
        T_xie_int_up =  abs(hotsu - mode_nms)*3 + mode_nms;
        
        % And the length as Otsu threshold
        T_xie_length_up = 3*max(thresholdOtsu(feature_Length),thresholdRosin(feature_Length));
        
        % Make a classification function as whether it is above the line
        F_classifer_up = @(nms,length,curvature,int)  (((T_xie_int_up + (T_xie_int_up/T_xie_length_up)*(-length) )<nms));
        
        
        % the down one
        T_xie_int_down =  abs(hotsu - mode_nms)*2 + mode_nms;
        T_xie_length_down = 2*max(thresholdOtsu(feature_Length),thresholdRosin(feature_Length));
        
        F_classifer_down = @(nms,length,curvature,int)  (((T_xie_int_down + (T_xie_int_down/T_xie_length_down)*(-length) )<nms));
        
        % the points in between the two classifier is the not sure ones that
        % requires annotation
        not_sure_ind = find(F_classifer_up(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)==0 & F_classifer_down(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)>0);
        good_ind = find(F_classifer_up(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)>0);
        
        F_classifer_notuse =  @(nms,length,curvature,int)  (((T_xie_int_down/2 + (T_xie_int_down/T_xie_length_down)*(-length) )<nms));
        bad_ind = find(F_classifer_down(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)==0 & F_classifer_notuse(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)>0);
        
        
        h1 = figure(1); hold off;
        imagescc(currentImg); hold on;
        scrsz = get(0,'ScreenSize');
        set(h1,'Position',scrsz);
        
        training_ind = randsample(not_sure_ind,min(length(not_sure_ind),training_sample_number),true,W_neighbor(not_sure_ind));
        training_ind_unique = unique(training_ind);
        
        if(length(training_ind_unique)<training_sample_number)
            not_sure_ind_again = setdiff(not_sure_ind,training_ind_unique);
            training_ind_again = randsample(not_sure_ind_again,min(length(not_sure_ind_again),training_sample_number-length(training_ind_unique)),true,W_neighbor(not_sure_ind_again));
            training_ind = [training_ind_unique;training_ind_again]';
        end
        
        training_bad_ind = [];
        training_good_ind = [];
        
        good_bad_label = [];
        i_ind = 1;
        
        for i_mark = 1 : 2*length(training_ind)
            i_ind
            iLine = training_ind(i_ind);
            [y,x] = find(labelMask==iLine);
            h1=figure(1);
            plot(x,y,'g.');
            title(['Frame: ',num2str(i_mark)]);
            saveas(h1,[FilamentSegmentationChannelOutputDir,'/train_g_',num2str(i_mark),'.jpg']);
            
            AA = [max(1, round(mean(x))-50), min(size(imageNMS,2), round(mean(x))+50), ...
                max(1, round(mean(y))-50), min(size(imageNMS,1), round(mean(y))+50)];
            axis(AA);
            ch = getkey();
            if(ch==28)
                i_ind = max(1, i_ind - 1);
            else
                i_ind = i_ind + 1;
                if(ch==32)
                    good_bad_label(i_ind) = 1;
                    plot(x,y,'r.');
                else
                    good_bad_label(i_ind) = 2;
                    plot(x,y,'b.');
                end
                saveas(h1,[FilamentSegmentationChannelOutputDir,'/train_rb_',num2str(i_mark),'.jpg']);
                
                if(ch==27 || i_ind > length(training_ind))
                    i_ind = i_ind-1;
                    good_bad_label = good_bad_label(1:length(training_ind));
                    break;
                end
            end
        end
    end
    
    
    figure(1); axis auto;
    
    sample_good_ind = randsample(good_ind,min(length(good_ind),200),true,W_neighbor(good_ind));
    sample_bad_ind = randsample(bad_ind,min(length(bad_ind),200),true,W_neighbor(bad_ind));
    
    if size(sample_good_ind,1)==1
        sample_good_ind = sample_good_ind';
    end
    
    if size(sample_bad_ind,1)==1
        sample_bad_ind = sample_bad_ind';
    end
    
    training_good_ind = sample_good_ind;
    training_bad_ind = sample_bad_ind;
    
    if(~isempty(find(good_bad_label(1:end)==1)))
        train_good_ind_from_training = training_ind(find(good_bad_label(1:end)==1));
        if size(train_good_ind_from_training,1)==1
            train_good_ind_from_training = train_good_ind_from_training';
        end
        training_good_ind = [sample_good_ind; train_good_ind_from_training];
    end
    
    if(~isempty(find(good_bad_label(1:end)==2)))
        train_bad_ind_from_training = training_ind(find(good_bad_label(1:end)==2));
        if size(train_bad_ind_from_training,1)==1
            train_bad_ind_from_training = train_bad_ind_from_training';
        end
        training_bad_ind = [sample_bad_ind; train_bad_ind_from_training];
    end
    
    for i_area = [training_good_ind' training_bad_ind' ]
        [all_y_i, all_x_i] = find(labelMask == i_area);
        % get the curvature for training samples only
        
        bw_i = zeros(size(bw_out));
        bw_i(sub2ind(size(bw_i), round(all_y_i),round(all_x_i)))=1;
        end_points_i = bwmorph(bw_i,'endpoints');
        [y_i, x_i]=find(end_points_i);
        
        if isempty(x_i)
            % if there is no end point, then it is a enclosed circle
            [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 1000, all_x_i(1),all_y_i(1));
        else
            [y_i, x_i]=find(end_points_i);
            [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 1000, x_i(1),y_i(1));
        end
        
        ordered_points{i_area} = [line_i_x, line_i_y];
        
        line_smooth_H = fspecial('gaussian',5,1.5);
        
        line_i_x = (imfilter(line_i_x, line_smooth_H, 'replicate', 'same'));
        line_i_y = (imfilter(line_i_y, line_smooth_H, 'replicate', 'same'));
        
        smoothed_ordered_points{i_area} = [line_i_x, line_i_y];
        
        Vertices = [line_i_x' line_i_y'];
        Lines=[(1:size(Vertices,1)-1)' (2:size(Vertices,1))'];
        k=LineCurvature2D(Vertices,Lines);
        
        feature_Curvature(i_area) = mean(k);
    end
    
    train_length_good = feature_Length(training_good_ind);
    train_length_bad = feature_Length(training_bad_ind);
    
    train_int_good = feature_MeanInt(training_good_ind);
    train_int_bad = feature_MeanInt(training_bad_ind);
    
    train_nms_good = feature_MeanNMS(training_good_ind);
    train_nms_bad = feature_MeanNMS(training_bad_ind);
    
    train_cur_good = feature_Curvature(training_good_ind);
    train_cur_bad = feature_Curvature(training_bad_ind);
    
    feature_good = [train_nms_good train_length_good train_cur_good train_int_good  ];
    feature_bad = [train_nms_bad train_length_bad train_cur_bad train_int_bad  ];
    
    label_good = ones(size(train_length_good));
    label_bad = -ones(size(train_length_bad));
    
    Classifier_Type_ind=funParams.Classifier_Type_ind;
    if(funParams.Classifier_Type_ind==2)
        % if user intended for SVM classifier
        %         Linear Kernel
        model_polynomial  = libsvmtrain([label_good; label_bad], [feature_good; feature_bad], '-t 1');
        [predict_label_L, accuracy_L, dec_values_L] = libsvmpredict([label_good; label_bad], [feature_good; feature_bad], model_polynomial);
        
        F_classifer_train_this_channel =  @(nms,length,curvature,int) (libsvmpredict(ones(size(nms)), [nms length curvature int], model_polynomial)>0);
        
        F_classifer_train_output{iChannel} = [FilamentSegmentationChannelOutputDir,'/F_classifer_channel.mat'];
        
        save([FilamentSegmentationChannelOutputDir,'/F_classifer_channel.mat'],'F_classifer_train_this_channel','Classifier_Type_ind');
        
    end
    
    % feature_training = [feature_Length(training_ind) feature_MeanInt(training_ind)];
    % label_training = good_bad_label;
    
    % mean_good = mean(feature_good);
    % mean_bad = mean(feature_bad);
    % v = corr_LDA([feature_good; feature_bad], 2, [length(label_good); length(label_bad)], 0.7);
    %
    % v*feature_good
    %
    % v*feature_bad
    %
    % F_classifer=[];
    %
    % save('F_classifer.mat','F_classifer');
    
    if(funParams.Classifier_Type_ind==1)
        % if the user want the linear classifier
        T_xie_int_mid = (T_xie_int_up + T_xie_int_down)/2;
        T_xie_length_mid = (T_xie_length_up +  T_xie_length_down)/2;
        
        F_classifer_mid = @(nms,length,curvature,int) (((T_xie_int_mid + (T_xie_int_mid/T_xie_length_mid)*(-length) )<nms));
        
        train_mat = [];
        for T_xie_int_grid = T_xie_int_down : (T_xie_int_up - T_xie_int_down)/10 : T_xie_int_up
            for T_xie_length_grid = T_xie_length_down : (T_xie_length_up - T_xie_length_down)/10 : T_xie_length_up
                
                F_classifer_train = @(nms,length,curvature,int) (((T_xie_int_grid + (T_xie_int_grid/T_xie_length_grid)*(-length) -nms )));
                train_mat = [train_mat; T_xie_int_grid T_xie_length_grid ...
                    (-sum(F_classifer_train(feature_MeanNMS(training_good_ind),...
                    feature_Length((training_good_ind))))...
                    +sum(F_classifer_train(feature_MeanNMS(training_bad_ind),...
                    feature_Length((training_bad_ind)))))];
            end
        end
        
        ind = find(train_mat(:,3)==max(train_mat(:,3)));
        
        T_xie_int_train = train_mat(ind(1), 1);
        T_xie_length_train = train_mat(ind(1), 2);
        clear train_mat training_good_ind training_bad_ind feature_Length feature_MeanNMS feature_good  feature_bad;
        clear imageNMS nms imageInt currentImg;
        clear F_classifer_down F_classifer_mid F_classifer_notuse F_classifer_train F_classifer_up;
        
        F_classifer_train_this_channel = @(nms,length,curvature,int) (((T_xie_int_train + (T_xie_int_train/T_xie_length_train)*(-length) )<nms));
        
        % due to some problem of the size of the function, with wired
        % difference of being here or outside the framework, now change to save
        % only the mat file name in the MD file.
        
    end
    
    F_classifer_train_output{iChannel} = [FilamentSegmentationChannelOutputDir,'/F_classifer_channel.mat'];
    
    save([FilamentSegmentationChannelOutputDir,'/F_classifer_channel.mat'],'F_classifer_train_this_channel','Classifier_Type_ind');
    
end

% save([FilamentSegmentationProcessOutputDir,'/F_classifer_process.mat'],'F_classifer_train_output');

