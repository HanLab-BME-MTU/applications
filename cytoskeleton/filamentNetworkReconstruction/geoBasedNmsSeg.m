function  [T_otsu, current_all_matching_bw, current_model ]  = geoBasedNmsSeg(imageNMS, imageInt, classifier_trained, graph_matching_flag,MaskCell,iFrame,FilamentSegmentationChannelOutputDir,funParams)
% geoBasedNmsSeg segments filaments from input image(nms) based on the geometrical features of the curves/lines in the image
% Input:
%    ImageIn:                           the input image, typically the non maximum supress version of the steerable filtering output
%    classifier_trained:                the trained or provided classifier of the curvee, if not provided, use empirical function
% Output:
%    T_otsu:                            the threshold defined by Otsu method for intensity of the input image, just as a format thing here.
%    current_all_matching_bw:           the segmented results, this serves as the starting point of the graphic matching
%

% Liya Ding
% 2013.02



% to show or not to show the figures; not affecting the saving part
 if(funParams.nofiguredisruption ==1 )
     set_visible='off';
 else
     set_visible='on';
 end

 % save figures or not, show messages or not
 SaveFigures = funParams.savestepfigures;
 ShowDetailMessages = funParams.savestepfigures;


% define a mask to reduce the number of curvelet to classify
if(~isempty(MaskCell) || mean(double(MaskCell))==1)
    
%     MaskCell =  imdilate(MaskCell,fspecial('disk', 3)>0);
else
    T_Rosin_otsu = thresholdRosin(imfilter(imageInt,fspecial('gaussian',11,2)));
    MaskCell = imageInt>T_Rosin_otsu*6/3;
    if(sum(sum(MaskCell))<size(MaskCell,1)*size(MaskCell,2)*0.05)
        MaskCell = imageInt>T_Rosin_otsu*5/3;
    end
    if(sum(sum(MaskCell))<size(MaskCell,1)*size(MaskCell,2)*0.05)
        MaskCell = imageInt>T_Rosin_otsu*4/3;
    end
    if(sum(sum(MaskCell))<size(MaskCell,1)*size(MaskCell,2)*0.05)
        MaskCell = imageInt>T_Rosin_otsu*3/3;
    end
    MaskCell = imdilate(MaskCell,fspecial('disk', 71)>0);
end

% MaskCell = ones(size(MaskCell));

% figure(1); hold off;
% subplot(121);
% imagescc(imageInt);
% 
% % imageNMS = imageNMS.*MaskCell;
% % imageInt = imageInt.*MaskCell;
% subplot(122);
% imagescc(MaskCell);

% pause;
%
% the threshold defined by Otsu method


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
ob_prop = regionprops(labelMask,'Area','MajorAxisLength','Eccentricity','MinorAxisLength','Centroid');

% Redefine variable for easy of notation
obAreas = [ob_prop.Area];

nLine = length(obAreas);

% Some feature for later consideration
obLongaxis = [ob_prop.MajorAxisLength];
obShortaxis = [ob_prop.MinorAxisLength];
obEccentricity = [ob_prop.Eccentricity];
obCentroid = zeros(2, length(obAreas));
obCentroid(:) = [ob_prop.Centroid];
% The ratio of short vs long axis
ratio  = obShortaxis./obLongaxis;


feature_MeanInt = nan(nLine,1);
feature_MeanNMS = nan(nLine,1);
feature_Length = obAreas';
smoothed_ordered_points = cell(1,1);
feature_Curvature = nan(nLine,1);

% for the features, only include those curves/lines longer than 4 pixels
ind_long = find(feature_Length>5);

% get the mean intensity of the curves
for i_area = ind_long'
    [all_y_i, all_x_i] = find(labelMask == i_area);
    NMS = imageNMS(sub2ind(size(bw_out), round(all_y_i),round(all_x_i)));
    feature_MeanNMS(i_area) = mean(NMS);
    INT = imageInt(sub2ind(size(bw_out), round(all_y_i),round(all_x_i)));
    feature_MeanInt(i_area) = mean(INT);
    % this version with the curvature measure, to save time.
    
    %     bw_i = zeros(size(bw_out));
    %     bw_i(sub2ind(size(bw_i), round(all_y_i),round(all_x_i)))=1;
    %     end_points_i = bwmorph(bw_i,'endpoints');
    %     [y_i, x_i]=find(end_points_i);
    %
    %     if isempty(x_i)
    %         % if there is no end point, then it is a enclosed circle
    %         [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 1000, all_x_i(1),all_y_i(1));
    %     else
    %         [y_i, x_i]=find(end_points_i);
    %         [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 1000, x_i(1),y_i(1));
    %     end
    %
    %     ordered_points{i_area} = [line_i_x, line_i_y];
    %
    %     line_smooth_H = fspecial('gaussian',5,1.5);
    %
    %     line_i_x = (imfilter(line_i_x, line_smooth_H, 'replicate', 'same'));
    %     line_i_y = (imfilter(line_i_y, line_smooth_H, 'replicate', 'same'));
    %      smoothed_ordered_points{i_area} = [line_i_x, line_i_y];
    %
    %     Vertices = [line_i_x' line_i_y'];
    %     Lines=[(1:size(Vertices,1)-1)' (2:size(Vertices,1))'];
    %     k=LineCurvature2D(Vertices,Lines);
    %
    %     feature_Curvature(i_area) = mean(k);
end

% figure; plot3(feature_Length,feature_MeanInt,feature_InvCurvature,'.');

if(isempty(classifier_trained))
    % if there is no input classifier, build one, with empirical setting
    
    % find the mode of the intensity of the curves/lines
    [hist_n,bin] = hist(feature_MeanNMS,200);
    ind_mode = find(hist_n==max(hist_n));
    mode_nms = bin(ind_mode(1));
    % And find the Otsu threshold for the intensity
    hotsu = thresholdOtsu(feature_MeanNMS(find(feature_MeanNMS>mode_nms)));
    
    % Set the slanted classification line cutoff as twice of the Otsu with
    % respect to the mode
    T_xie_int =  abs(hotsu - mode_nms)*1.8+ mode_nms;
    
    % And the length as Otsu threshold
    T_xie_length = 1.8*max(thresholdOtsu(feature_Length),thresholdRosin(feature_Length));
   
    % Make a classification function as whether it is above the line
    T_xie_int_train = T_xie_int;
    T_xie_length_train = T_xie_length;
    
    F_classifer = @(int,length) (((T_xie_int_train + (T_xie_int_train/T_xie_length_train)*(-length) )<int));
    
else
    % when there is an input classifer, use the input one
    F_classifer = classifier_trained;
end


master_flassier = F_classifer;

% Select those above the line as our good ones as starting point of graph
% matching
Good_ind = find(F_classifer(feature_MeanNMS, feature_Length)>0);
Bad_ind = find(F_classifer(feature_MeanNMS, feature_Length)==0);
if(ShowDetailMessages==1)
    display(['Length and NMS only: Number of good ones: ',num2str(length(Good_ind)),', number of bad ones: ',num2str(length(Bad_ind))]);
end

% get the mean intensity of the curves
for i_area = Good_ind'
    [all_y_i, all_x_i] = find(labelMask == i_area);
    % only for those good ones, check the curvature
    
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

feature_all = [];

feature_all.feature_MeanNMS = feature_MeanNMS;
feature_all.feature_Length = feature_Length;
feature_all.feature_MeanInt = feature_MeanInt;
feature_all.feature_Curvature = feature_Curvature;


Good_ind = find(F_classifer(feature_MeanNMS, feature_Length)>0 & feature_Curvature<0.051);
Bad_ind = find(F_classifer(feature_MeanNMS, feature_Length)==0 | feature_Curvature>=0.051);

if(ShowDetailMessages==1)
    display(['Length,NMSand Curvature: Number of good ones: ',num2str(length(Good_ind)),', number of bad ones: ',num2str(length(Bad_ind))]);
end
    
%     
%     h12 =  figure(12);hold off;
%     plot(feature_Length(Good_ind),feature_MeanNMS(Good_ind),'b.'); hold on;
%     plot(feature_Length(Bad_ind),feature_MeanNMS(Bad_ind),'r.');
%     title('Classifier Plane with matched data round 0');
%     saveas(h12,['./GEO_frame_',num2str(iFrame),'_round0_trained_plane.tif']);

% plot the output image with these good ones
current_all_seg_bw = zeros(size(labelMask));
current_model = [];
model_ind = [];

for i_E = 1 : length(Good_ind)
    current_good_bw = labelMask==Good_ind(i_E);
    current_all_seg_bw = or(current_all_seg_bw, current_good_bw);
    [y_i, x_i] = find(labelMask==Good_ind(i_E));
    
    current_model{i_E} = [x_i y_i];
    model_ind{i_E} = Good_ind(i_E);
end

current_all_matching_bw = current_all_seg_bw;

% restore parameters from saved trained functions
for int_test = max(feature_MeanNMS):-0.1:0
    if F_classifer(int_test,0)==0;
        T_xie_int_train = int_test;
        break;
    end
end

for length_test = max(feature_Length):-0.1:0
    if F_classifer(0,length_test)==0;
        T_xie_length_train = length_test;
        break;
    end
end

% graph_matching_flag=1;

% Now if user intended for graph matching part, do it
if(graph_matching_flag==1)
    
    h3 = figure(3); set(h3,'Visible',set_visible);
    imagescc(imageInt);hold on;
    [y,x]=find(current_all_seg_bw>0);
    plot(x,y,'.','MarkerSize',3);
    
    if(~exist([FilamentSegmentationChannelOutputDir,'/GEO'],'dir'))
        mkdir(FilamentSegmentationChannelOutputDir,'GEO');
    end
    
    if(SaveFigures==1)        
        saveas(h3,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round0_all_match_bw.tif']);
    end
    
    confidency_interval = 0.7;
%     close(h1);
        iIteration=0;
    [current_model,current_matching_bw,model_ind]...
    = graph_matching_linking_once([], current_all_seg_bw, confidency_interval,imageInt, ...
                                Good_ind,ind_long, [],feature_all,labelMask,master_flassier,iIteration,funParams);
    if(SaveFigures==1)
        imwrite(double(current_all_seg_bw*3/4),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round1_begin.tif']);
        imwrite(double(current_matching_bw/2),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round1_end.tif']);
        h1=figure(1);set(h1,'Visible',set_visible);saveas(h1,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round1_match_color.tif']);
        h3=figure(3);set(h3,'Visible',set_visible);saveas(h3,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round1_all_match_bw.tif']);
    end

    retrain_flag = 1;
     
    T_xie_int_train = T_xie_int_train*0.9;
    T_xie_length_train = T_xie_length_train*0.9;
    F_classifer = @(int,length) (((T_xie_int_train + (T_xie_int_train/T_xie_length_train)*(-length) )<int));
    Good_ind = find(F_classifer(feature_MeanNMS, feature_Length)>0);
    Bad_ind = find(F_classifer(feature_MeanNMS, feature_Length)==0);
    
    h12 =  figure(12); set(h12,'Visible',set_visible);hold off;
    plot(feature_Length(Good_ind),feature_MeanNMS(Good_ind),'b.'); hold on;
    plot(feature_Length(Bad_ind),feature_MeanNMS(Bad_ind),'r.');
    
    title('Classifier Plane with matched data round 1');
    
    if(SaveFigures==1)
        saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round1_lowered_plane.tif']);
    end
    
    % plot the output image with these good ones
    current_all_seg_bw = zeros(size(labelMask));
    for i_E = 1 : length(Good_ind)
        current_good_bw = labelMask==Good_ind(i_E);
        current_all_seg_bw = or(current_all_seg_bw, current_good_bw);
    end
    
     if(SaveFigures==1)   
        imwrite(double(current_all_seg_bw*3/4),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round2_begin.tif']);
     end
     
     if(exist('h1','var'))
     close(h1);
     end
    iIteration=1;
    [current_model,current_matching_bw, model_ind]...
            = graph_matching_linking_once(current_model, current_all_seg_bw, confidency_interval,imageInt, ...
                                Good_ind,ind_long, model_ind,feature_all,labelMask,master_flassier,iIteration,funParams);

    if(SaveFigures==1)        
        imwrite(double(current_matching_bw/2),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round2_end.tif']);
        h1=figure(1);set(h1,'Visible',set_visible);saveas(h1,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round2_match_color.tif']);
        h3=figure(3);set(h3,'Visible',set_visible);saveas(h3,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round2_all_match_bw.tif']);
    end
    
    good_bw = nms_seg_no_brancing.*current_matching_bw;
    
    Matched_ind=[];
    UnMatched_ind=[];
    
    for i_area = ind_long'
        
        [all_y_i, all_x_i] = find(labelMask == i_area);
        if_matched = mean(good_bw(sub2ind(size(bw_out), round(all_y_i),round(all_x_i))));
        if(if_matched>0.1)
            Matched_ind = [Matched_ind; i_area];
        else
            UnMatched_ind = [UnMatched_ind; i_area];
        end
        
    end

%     retrain_flag=0;
%     
%     for iIteration = 2 : 3
%         graph_matching_iteration_in_segmentation;
%     end
    current_all_matching_bw = current_matching_bw>0;
end

