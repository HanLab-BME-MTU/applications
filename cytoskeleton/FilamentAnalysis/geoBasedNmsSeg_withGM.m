function  [T_otsu, current_all_matching_bw, current_model ]  = geoBasedNmsSeg_withGM(imageNMS, imageInt, classifier_trained, graph_matching_flag,MaskCell,iFrame,FilamentSegmentationChannelOutputDir,funParams,iChannel)
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

CoefAlpha = funParams.CoefAlpha;
LengthThreshold = funParams.LengthThreshold;
IternationNumber = funParams.IternationNumber;
CurvatureThreshold = funParams.CurvatureThreshold;

%% % these are for new GUI allowing different setting of different channels
if(length(CoefAlpha)>=iChannel)
    CoefAlpha = CoefAlpha(iChannel);
end
if(length(LengthThreshold)>=iChannel)
    LengthThreshold = LengthThreshold(iChannel);
end
if(length(IternationNumber)>=iChannel)
    IternationNumber = IternationNumber(iChannel);
end
if(length(CurvatureThreshold)>=iChannel)
    CurvatureThreshold = CurvatureThreshold(iChannel);
end
%%

Good_ind_cell=cell(1,4);
model_ind_cell =cell(1,10);


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

% figure(1);set(h1,'Visible',set_visible); hold off;
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
T_otsu_start =  max(mode_nms/2, (-abs(T_otsu - mode_nms)*0.1+mode_nms));

MaskCell=MaskCell>0;

T_otsu_start=0.4;

MaskCell_nohole = MaskCell;
MaskCell_nohole = bwmorph(MaskCell,'fill');
MaskCell_nohole(140:170,140:170)=1;

imageNMS_nohole = imageNMS.*MaskCell_nohole;
imageInt_nohole = imageInt.*MaskCell_nohole;

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
ind_long = find(feature_Length>=LengthThreshold);


h14 = figure(14); set(h14,'Visible',set_visible); hold off;
display_labelMask = 0*labelMask;

for iiii = ind_long'
    display_labelMask(labelMask==iiii)=iiii+100;
end

imagesc(display_labelMask);
hold on;
for iiii = ind_long'
    %     obCentroid(1,iiii)
    text(obCentroid(1,iiii),obCentroid(2,iiii),num2str(iiii),'color','r');
end
if(SaveFigures==1)
    
    if(  ~exist([FilamentSegmentationChannelOutputDir,'/GEO'],'dir'))
        mkdir([FilamentSegmentationChannelOutputDir,'/GEO']);
    end
    
    
    saveas(h14,[FilamentSegmentationChannelOutputDir,'/GEO/numbers_filament_f',num2str(iFrame),'.fig']);
    saveas(h14,[FilamentSegmentationChannelOutputDir,'/GEO/numbers_filament_f',num2str(iFrame),'.tif']);
    saveas(h14,[FilamentSegmentationChannelOutputDir,'/GEO/numbers_filament_f',num2str(iFrame),'.eps']);
end






% get the mean intensity of the curves
for i_area = ind_long'
    [all_y_i, all_x_i] = find(labelMask == i_area);
    NMS = imageNMS(sub2ind(size(bw_out), round(all_y_i),round(all_x_i)));
    %     feature_MeanNMS(i_area) = mean(NMS);
    %     feature_MeanNMS(i_area) = max(NMS);
    feature_MeanNMS(i_area) = median(NMS);
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
    try
        hRosin = thresholdRosin(feature_MeanNMS(find(feature_MeanNMS>mode_nms)));
        hRosin_or_otsu = hRosin;
    catch
        hotsu = thresholdOtsu(feature_MeanNMS(find(feature_MeanNMS>mode_nms)));
        hRosin_or_otsu = hotsu;
    end
    
    % Set the slanted classification line cutoff as twice of the Otsu with
    % respect to the mode
    T_xie_int_this =  abs(hRosin_or_otsu - mode_nms)*CoefAlpha+ mode_nms;
    
    whole_movie_otsu =  funParams.Whole_movie_stat_cell{iChannel}.otsu_mode_NMS;
    whole_movie_mode =  funParams.Whole_movie_stat_cell{iChannel}.mode_NMS;
    
    whole_T_xie_int =  abs(whole_movie_otsu - whole_movie_mode)*CoefAlpha+ whole_movie_mode;
    
    weight_whole = (-0.5)*funParams.Whole_movie_ind(iChannel) + 1.5;
    weight_this =  1 - weight_whole;
    
    T_xie_int = weight_whole*whole_T_xie_int + weight_this*T_xie_int_this;
    
    
    % And the length as Otsu threshold
    %     T_xie_length =CoefAlpha*max(thresholdOtsu(feature_Length),thresholdRosin(feature_Length));
    T_xie_length = CoefAlpha*(thresholdRosin(feature_Length));
    
    % Make a classification function as whether it is above the line
    T_xie_int_train = T_xie_int;
    T_xie_length_train = weight_whole*(CoefAlpha*(funParams.Whole_movie_stat_cell{iChannel}.rosin_mode_Length)) + weight_this*T_xie_length;
    
    
    %     thresholdOtsu(feature_Length)
    %     thresholdRosin(feature_Length)
    %
    %     T_xie_length_train
    
    F_classifer = @(nms,length,int,curv) (((T_xie_int_train + (T_xie_int_train/T_xie_length_train)*(-length) )<nms));
    
else
    % when there is an input classifer, use the input one
    F_classifer = classifier_trained;
end


master_flassier = F_classifer;

% Select those above the line as our good ones as starting point of graph
% matching

Good_ind = find(F_classifer(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)>0);
Bad_ind = find(F_classifer(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)==0);

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
    
    feature_Curvature(i_area) = mean(abs(k));
    
end

feature_all = [];

feature_all.feature_MeanNMS = feature_MeanNMS;
feature_all.feature_Length = feature_Length;
feature_all.feature_MeanInt = feature_MeanInt;
feature_all.feature_Curvature = feature_Curvature;

IntThreshold = thresholdRosin(feature_all.feature_MeanInt)/2;

Good_ind = find(F_classifer(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)>0 & feature_Curvature<CurvatureThreshold & feature_Length>=LengthThreshold);
Bad_ind = find(F_classifer(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)==0 | feature_Curvature>=CurvatureThreshold | feature_Length<LengthThreshold);

F2_Classifier = @(nms,length,int,curv) (curv<=CurvatureThreshold & int> IntThreshold);
    
master_flassier = @(nms,length,int,curv)(F2_Classifier(nms,length,int,curv) >0 & F_classifer(nms,length,int,curv) >0 );

% CurvatureThreshold=0.2;
changed_ind_good_K = find(F_classifer(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)>0 ...
    & feature_Curvature >= CurvatureThreshold & feature_Length>LengthThreshold);

changed_ind_good_I = find(F_classifer(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)>0 ...
    & feature_MeanInt <= IntThreshold*5/5 & feature_Length>LengthThreshold);

if(SaveFigures==1)
    plot_detailed_high_low_confidence_fragments;    
end

Original_set_Good_ind = Good_ind;
Original_set_Bad_ind = Bad_ind;

final_set_good = Good_ind;
final_set_rescue = [];
final_set_connection = cell(1,1);
connect_count = 0;

h12 = figure(12);set(h12,'Visible',set_visible);hold off;
plot(feature_Length(Original_set_Bad_ind),feature_MeanNMS(Original_set_Bad_ind),'r.');hold on;
plot(feature_Length(Original_set_Good_ind),feature_MeanNMS(Original_set_Good_ind),'b.');
for iiii = ind_long'
    text(feature_Length(iiii),feature_MeanNMS(iiii),num2str(iiii));
end
if(SaveFigures==1)
    saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/number_feature_all_', ...
        num2str(iFrame),'.fig']);
    saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/number_feature_all_', ...
        num2str(iFrame),'.tif']);
    saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/number_feature_all_', ...
        num2str(iFrame),'.eps']);
end

if(ShowDetailMessages==1)
    display(['Length,NMSand Curvature: Number of good ones: ',num2str(length(Good_ind)),', number of bad ones: ',num2str(length(Bad_ind))]);
end



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
    if F_classifer(int_test,0,0,0)==0;
        T_xie_int_train = int_test;
        break;
    end
end

for length_test = max(feature_Length):-0.1:0
    if F_classifer(0,length_test,0,0)==0;
        T_xie_length_train = length_test;
        break;
    end
end
%
% graph_matching_flag=0;

Good_ind_cell{1} = Good_ind;
Bad_ind_cell{1} = Bad_ind;



% Now if user intended for graph matching part, do it
if(graph_matching_flag==1)
    tomatch_ind_cell{1} = Good_ind;
    
    h3 = figure(3); set(h3,'Visible',set_visible);
    imagescc(imageInt);hold on;
    [y,x]=find(current_all_seg_bw>0);
    plot(x,y,'.','MarkerSize',3);
    
    if(~exist([FilamentSegmentationChannelOutputDir,'/GEO'],'dir'))
        mkdir(FilamentSegmentationChannelOutputDir,'GEO');
    end
    
    if(SaveFigures==1)
        saveas(h3,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round0_all_match_bw.tif']);
        saveas(h3,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round0_all_match_bw.fig']);
        saveas(h3,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round0_all_match_bw.eps']);
    end
    
    
    confidency_interval = 0.9;
    %     close(h1);
    iIteration=0;
    [current_model,current_matching_bw,model_ind,transformer, final_set_rescue, final_set_connection, connect_count,bw_rgb_three_color,bw_whitergb_three_color]...
        = graph_matching_linking_once([], current_all_seg_bw, confidency_interval,imageInt, ...
        Good_ind,ind_long, [],feature_all,labelMask,master_flassier,iIteration,funParams,final_set_rescue, final_set_connection, connect_count,Original_set_Good_ind);
    
    model_ind_cell{1} = model_ind;
    current_model_cell{1} = current_model;
    Matched_ind=[];
    UnMatched_ind=[];
    good_bw = nms_seg_no_brancing.*current_matching_bw;
    
    Matched_ind = model_ind(find(~isnan(model_ind)));
    Matched_ind_cell{1} = Matched_ind;
    UnMatched_ind_cell{1} = setdiff(ind_long,Matched_ind);
    
     imwrite(double(bw_whitergb_three_color),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_whitethreecolor_pixel.tif']);
     
    
    if(SaveFigures==1)
        imwrite(double(current_all_seg_bw*3/4),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_begin.tif']);
        imwrite(double(current_matching_bw/2),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_end.tif']);
        imwrite(double(bw_rgb_three_color),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_threecolor_pixel.tif']);
       h1=figure(1);set(h1,'Visible',set_visible);
        saveas(h1,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_match_color.tif']);
        saveas(h1,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_match_color.fig']);
        print(h1,'-depsc',[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_match_color.eps']);
        h3=figure(3);set(h3,'Visible',set_visible);
        saveas(h3,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_all_match_bw.tif']);
        saveas(h3,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_all_match_bw.fig']);
        print(h3,'-depsc',[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_all_match_bw.eps']);
        h11=figure(11);set(h11,'Visible',set_visible);
        saveas(h11,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_threecolor.tif']);
        saveas(h11,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_threecolor.fig']);
        print(h11,'-depsc',[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_threecolor.eps']);
        h21=figure(21);set(h21,'Visible',set_visible);
        saveas(h21,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_whitethreecolor.tif']);
        saveas(h21,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_whitethreecolor.fig']);
        print(h21,'-depsc',[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_whitethreecolor.eps']);
        
        h12 = figure(12);set(h12,'Visible',set_visible);hold off;
        plot(feature_Length(Original_set_Bad_ind),feature_MeanNMS(Original_set_Bad_ind),'r.');hold on;
        plot(feature_Length(Matched_ind),feature_MeanNMS(Matched_ind),'g.');
        plot(feature_Length(Original_set_Good_ind),feature_MeanNMS(Original_set_Good_ind),'b.');
        
        title(['Classifier Plane with matched data after round ',num2str(iIteration)]);
        print(h12,'-depsc2',[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_afterround',num2str(iIteration),'_match_plane.eps']);
        print(h12,'-dtiff',[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_afterround',num2str(iIteration),'_match_plane.tiff']);
        saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_afterround',num2str(iIteration),'_match_plane.fig']);
        
        for iiii = ind_long'
            text(feature_Length(iiii),feature_MeanNMS(iiii),num2str(iiii));
        end
        
        saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/number_feature_',num2str(iFrame),num2str(iIteration),'.fig']);
        
        h13 = figure(13);set(h13,'Visible',set_visible);
        RGB_seg_orient_heat_map = heat_display_filament_from_model(imageInt, current_model);
        imshow(RGB_seg_orient_heat_map);
        title(['Heat display for Segmentation after round ',num2str(iIteration)]);
        imwrite(RGB_seg_orient_heat_map,[FilamentSegmentationChannelOutputDir,'/GEO/heat_frame_',num2str(iFrame),'_afterround',num2str(iIteration),'.tif']);
        imwrite(RGB_seg_orient_heat_map,[FilamentSegmentationChannelOutputDir,'/GEO/heat_frame_',num2str(iFrame),'_afterround',num2str(iIteration),'.bmp']);
        
    end
    
    
    
    if(IternationNumber>1)
        retrain_flag = 1;
        
        T_xie_int_train = T_xie_int_train*0.95;
        T_xie_length_train = T_xie_length_train*0.95;
        F_classifer = @(nms,length,curv,int) (((T_xie_int_train + (T_xie_int_train/T_xie_length_train)*(-length) )<nms));
        
        Good_ind = find(F_classifer(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)>0 & feature_Curvature<CurvatureThreshold & feature_Length>=LengthThreshold);
        Bad_ind = find(F_classifer(feature_MeanNMS, feature_Length,feature_Curvature,feature_MeanInt)==0 | feature_Curvature>=CurvatureThreshold | feature_Length<LengthThreshold);
        
        Good_ind_cell{2} = Good_ind;
        Bad_ind_cell{2} = Bad_ind;
        
        
        if(SaveFigures==1)
            h12 =  figure(12); set(h12,'Visible',set_visible);hold off;
            plot(feature_Length(Good_ind),feature_MeanNMS(Good_ind),'b.'); hold on;
            plot(feature_Length(Bad_ind),feature_MeanNMS(Bad_ind),'r.');
            
            title('Classifier Plane with matched data round 1');
            saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round1_lowered_plane.tif']);
            saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round1_lowered_plane.fig']);
            saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round1_lowered_plane.eps']);
        end
        
        % plot the output image with these good ones
        current_all_seg_bw = zeros(size(labelMask));
        for i_E = 1 : length(Good_ind)
            current_good_bw = labelMask==Good_ind(i_E);
            current_all_seg_bw = or(current_all_seg_bw, current_good_bw);
        end
        
        if(SaveFigures==1)
            imwrite(double(current_all_seg_bw*3/4),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round2_begin.tif']);
            imwrite(double(bw_rgb_three_color),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_threecolor_pixel.tif']);
            h12 = figure(12);set(h12,'Visible',set_visible);hold off;
            plot(feature_Length(Bad_ind),feature_MeanNMS(Bad_ind),'r.');hold on;
            plot(feature_Length(Good_ind),feature_MeanNMS(Good_ind),'b.');
            
            title(['Classifier Plane with matched data before round ',num2str(iIteration)]);
            saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_beforeround',num2str(iIteration),'_match_plane.tif']);
            saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_beforeround',num2str(iIteration),'_match_plane.eps']);
            saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_beforeround',num2str(iIteration),'_match_plane.fig']);
        end
        
        if(exist('h1','var'))
            close(h1);
        end
        iIteration=1;
        [current_model,current_matching_bw, model_ind,transformer,final_set_rescue, final_set_connection, connect_count,bw_rgb_three_color,bw_whitergb_three_color]...
            = graph_matching_linking_once(current_model, current_all_seg_bw, confidency_interval,imageInt, ...
            Good_ind,ind_long, model_ind,feature_all,labelMask,master_flassier,iIteration,funParams,final_set_rescue, final_set_connection, connect_count,Original_set_Good_ind);
        
        imwrite(double(bw_whitergb_three_color),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_whitethreecolor_pixel.tif']);
     
        if(SaveFigures==1)
            imwrite(double(current_matching_bw/2),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_end.tif']);
            imwrite(double(bw_rgb_three_color),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_threecolor_pixel.tif']);
            imwrite(double(bw_whitergb_three_color),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_whitethreecolor_pixel.tif']);
            h1=figure(1);set(h1,'Visible',set_visible);
            saveas(h1,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_match_color.tif']);
            saveas(h1,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_match_color.fig']);
            print(h1,'-depsc',[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_match_color.eps']);
            h3=figure(3);set(h3,'Visible',set_visible);
            saveas(h3,[FilamentSegmentationChannelOutputDir,'/GEO/output_frame_',num2str(iFrame),'_round',num2str(iIteration),'_all_match_bw.tif']);
            saveas(h3,[FilamentSegmentationChannelOutputDir,'/GEO/output_frame_',num2str(iFrame),'_round',num2str(iIteration),'_all_match_bw.fig']);
            print(h3,'-depsc',[FilamentSegmentationChannelOutputDir,'/GEO/output_frame_',num2str(iFrame),'_round',num2str(iIteration),'_all_match_bw.eps']);
            h11=figure(11);set(h11,'Visible',set_visible);
            saveas(h11,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_threecolor.tif']);
            saveas(h11,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_threecolor.fig']);
            print(h11,'-depsc',[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_threecolor.eps']);
            h21=figure(21);set(h21,'Visible',set_visible);
            saveas(h21,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_whitethreecolor.tif']);
            saveas(h21,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_whitethreecolor.fig']);
            print(h21,'-depsc',[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration),'_whitethreecolor.eps']);
        end
        
        good_bw = nms_seg_no_brancing.*current_matching_bw;
        
        model_ind_cell{2} = model_ind;
        current_model_cell{2} = current_model;
        
        Matched_ind = model_ind(find(~isnan(model_ind)));
        
        Matched_ind_cell{2} = Matched_ind;
        UnMatched_ind_cell{2} = setdiff(ind_long,Matched_ind);
        
        original_Matched_ind = Matched_ind_cell{1};
        
        if(SaveFigures==1)
            
            h12 = figure(12);set(h12,'Visible',set_visible);hold off;
            plot(feature_Length(Original_set_Bad_ind),feature_MeanNMS(Original_set_Bad_ind),'r.');hold on;
            plot(feature_Length(Matched_ind),feature_MeanNMS(Matched_ind),'g.');
            plot(feature_Length(setdiff(original_Matched_ind,Original_set_Good_ind)), ...
                feature_MeanNMS(setdiff(original_Matched_ind,Original_set_Good_ind)),'m*');
            plot(feature_Length(Original_set_Good_ind),feature_MeanNMS(Original_set_Good_ind),'b.');
            title(['Classifier Plane with matched data after round ',num2str(iIteration)]);
            print(h12,'-depsc2',[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_afterround',num2str(iIteration),'_match_plane.eps']);
            print(h12,'-dtiff',[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_afterround',num2str(iIteration),'_match_plane.tiff']);
            saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_afterround',num2str(iIteration),'_match_plane.eps']);
            
            for iiii = setdiff(ind_long,UnMatched_ind_cell{iIteration+1})
                text(feature_Length(iiii),feature_MeanNMS(iiii),num2str(iiii));
            end
            
            %            saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/number_feature_',num2str(iFrame),'_round',num2str(iIteration),'.fig']);
            saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/number_feature_',num2str(iFrame),num2str(iIteration),'.fig']);
            saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/number_feature_',num2str(iFrame),num2str(iIteration),'.tif']);
            saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/number_feature_',num2str(iFrame),num2str(iIteration),'.eps']);
            
            h13 = figure(13);set(h13,'Visible',set_visible);
            RGB_seg_orient_heat_map = heat_display_filament_from_model(imageInt, current_model);
            imshow(RGB_seg_orient_heat_map);
            title(['Heat display for Segmentation after round ',num2str(iIteration)]);
            imwrite(RGB_seg_orient_heat_map,[FilamentSegmentationChannelOutputDir,'/GEO/heat_frame_',num2str(iFrame),'_afterround',num2str(iIteration),'.tif']);
            imwrite(RGB_seg_orient_heat_map,[FilamentSegmentationChannelOutputDir,'/GEO/heat_frame_',num2str(iFrame),'_afterround',num2str(iIteration),'.bmp']);
            
            
            
        end
        
    end
    
    retrain_flag=0;
    
    if(IternationNumber>2)
        for iIteration = 2 : IternationNumber
            graph_matching_iteration_in_segmentation;
        end
    end
    
    current_all_matching_bw = current_matching_bw>0;
end
