function [output_model,current_matching_bw,new_model_ind,transformer,final_set_rescue, final_set_connection, connect_count,bw_rgb_three_color,bw_whitergb_three_color] = ...
    graph_matching_linking_once(current_model, bw_to_be_matched,  confidency_interval,currentImg,...
    Good_ind,ind_long, model_ind,feature_all,labelMask,master_flassier,...
    iIteration,funParams,final_set_rescue, final_set_connection, connect_count,Original_set_Good_ind)
% Graph matching one iteration


new_model_ind=[];
output_model=current_model;
transformer=[];
bw_rgb_three_color=[];
bw_whitergb_three_color=[];

%%
% to show or not to show the figures; not affecting the saving part
if(funParams.nofiguredisruption ==1 )
    set_visible='off';
else
    set_visible='on';
end

% save figures or not, show messages or not
SaveFigures = funParams.savestepfigures;
ShowDetailMessages = funParams.savestepfigures;

%%
orig_current_model = current_model;
count_linking=0;
tomatch_ind = [];
Matched_ind=[];
UnMatched_ind=[];

% plot the output image with these good ones
current_matching_bw = zeros(size(bw_to_be_matched));

% if there is already a current model, build the image for matched lines
if(~isempty(current_model))
    for i_E = 1 : length(current_model)
        y_i = current_model{i_E}(:,2);
        x_i = current_model{i_E}(:,1);
        current_matching_bw(sub2ind(size(current_matching_bw), round(y_i),round(x_i)))=1;
        %         figure(2);imagesc(current_matching_bw);title(num2str(i_E));
        %         pause;
    end
    count_linking = length(current_model);
    
    
    % the input Good_ind includes those already been matched ones, so get rid of them
    for i_area = Good_ind'
        [all_y_i, all_x_i] = find(labelMask == i_area);
        if_matched = mean(current_matching_bw(sub2ind(size(currentImg), round(all_y_i),round(all_x_i))));
        if(if_matched>0.1)
            Matched_ind = [Matched_ind; i_area];
        else
            UnMatched_ind = [UnMatched_ind; i_area];
        end
    end
    tomatch_ind = UnMatched_ind;
else
    % if there is no current model, then this is the first matching, try
    % match every good ones.
    tomatch_ind = Good_ind;
    
end

if(ShowDetailMessages==1)
    display(['iIteration ',num2str(iIteration),' start:']);
    display(['Matching: Number of good ones: ',num2str(length(Good_ind)),', number of to match ones: ',num2str(length(tomatch_ind))]);
end
% get the new unmatched bw image, if there is a model, this is the part
% that is not matched yet, if there is not current model, this is all the
% lines in the to be matched bw image.
% new_unmatched_bw = bw_to_be_matched - current_matching_bw;

% Get properties for each of curve
ob_prop = regionprops(labelMask,'Area','MajorAxisLength','Eccentricity','MinorAxisLength','Centroid');

% Redefine variable for easy of notation
obAreas = [ob_prop.Area];
obLongaxis = [ob_prop.MajorAxisLength];
obShortaxis = [ob_prop.MinorAxisLength];
obEccentricity = [ob_prop.Eccentricity];

obCentroid = zeros(2, length(obAreas));
obCentroid(:) = [ob_prop.Centroid];

% The ratio of short vs long axis
ratio  = obShortaxis./obLongaxis;

nms_seg_no_brancing  = zeros(size(labelMask));

for i_ind = 1 : length(tomatch_ind)
    nms_seg_no_brancing(labelMask==tomatch_ind(i_ind))=1;
end

%%
% Use kdtree to limit pairs
% tic
% Input data
end_points_all = bwmorph(nms_seg_no_brancing,'endpoints');
ordered_end_points=[];

ordered_curve_xy = cell(2*(length(tomatch_ind)+length(current_model)),1);
ordered_smooth_curve_xy = cell(2*(length(tomatch_ind)+length(current_model)),1);
end_angles = nan(2*(length(tomatch_ind)+length(current_model)),1);
line_smooth_H = fspecial('gaussian',5,2);
circle_flag = zeros(2*(length(tomatch_ind)+length(current_model)),1);


model_ind_all = nan(length(tomatch_ind)+size(model_ind,1),20);


% for every new curve to be matched
for i_area = 1 : length(tomatch_ind)
    % the index of curvelet in every part
    model_ind_all(i_area,1)  = tomatch_ind(i_area);
    
    bw_i = labelMask == tomatch_ind(i_area);
    end_points_this = bwmorph(bw_i,'endpoints');
    [y,x] = find(end_points_this);
    if(isempty(x))
        ordered_end_points = [ordered_end_points; [-100, -100; -100, -100];];
        circle_flag(2*i_area-1:2*i_area)=1;
    else
        ordered_end_points = [ordered_end_points; [x, y];];
        [line_i_x, line_i_y] = line_following_with_limit(bw_i, sum(sum(bw_i)), x(1),y(1));
        line_i_x = line_i_x';
        line_i_y = line_i_y';
        
        ordered_curve_xy{2*i_area-1,1} = [line_i_x line_i_y];
        ordered_curve_xy{2*i_area,1} = flipud([line_i_x line_i_y]);
        
        line_i_x = (imfilter(line_i_x, line_smooth_H, 'replicate', 'same'));
        line_i_y = (imfilter(line_i_y, line_smooth_H, 'replicate', 'same'));
        
        ordered_smooth_curve_xy{2*i_area-1,1} = [line_i_x line_i_y];
        ordered_smooth_curve_xy{2*i_area,1} = flipud([line_i_x line_i_y]);
        
        length_line = length(line_i_x);
        
        end_angles(2*i_area-1) = atan2(line_i_x(1)-line_i_x(min(3,length_line)),line_i_y(1)-line_i_y(min(3,length_line)));
        end_angles(2*i_area) = atan2(line_i_x(end)-line_i_x(max(1,length_line-2)),line_i_y(end)-line_i_y(max(1,length_line-2)));
    end
end

if(~isempty(current_model))
    % for every matched curve to from previous matching
    for i_E = 1 : length(current_model)
        % the index of curvelet in current_model as obtained in previous iterations
        model_ind_all(length(tomatch_ind)+i_E,1:length(model_ind(i_E,:)))  =  model_ind(i_E,:);
        
        y_i = current_model{i_E}(:,2);
        x_i = current_model{i_E}(:,1);
        end_points_this = [x_i(1) y_i(1); x_i(end) y_i(end)];
        ordered_end_points = [ordered_end_points; end_points_this;];
        ordered_curve_xy{2*(length(tomatch_ind)+i_E)-1,1} = [x_i, y_i];
        ordered_curve_xy{2*(length(tomatch_ind)+i_E),1} = flipud([x_i, y_i]);
        
        line_i_x = (imfilter(x_i, line_smooth_H, 'replicate', 'same'));
        line_i_y = (imfilter(y_i, line_smooth_H, 'replicate', 'same'));
        
        ordered_smooth_curve_xy{2*(length(tomatch_ind)+i_E)-1,1} =  [line_i_x, line_i_y];
        ordered_smooth_curve_xy{2*(length(tomatch_ind)+i_E),1} = flipud( [line_i_x, line_i_y]);
        
        length_line = length(line_i_x);
        
        end_angles(2*(length(tomatch_ind)+i_E)-1) = atan2(line_i_x(1)-line_i_x(min(3,length_line)),line_i_y(1)-line_i_y(min(3,length_line)));
        end_angles(2*(length(tomatch_ind)+i_E)) = atan2(line_i_x(end)-line_i_x(max(1,length_line-2)),line_i_y(end)-line_i_y(max(1,length_line-2)));
    end
end

%%
% temp fix problem of ind matrix big
for ti = 1: size(model_ind_all,1)
    ind_line = model_ind_all(ti,:);
    
    ind_line= ind_line(~isnan(ind_line));
    ind_line =ind_line(ind_line~=0);
    
    
    model_ind_all(ti,1:length(ind_line)) = ind_line;
    model_ind_all(ti,length(ind_line)+1:end) = nan;
    
end

model_ind_all = model_ind_all(:,1:20);
%%

% Build tree (Memory automatically freed)
kd = KDTree(ordered_end_points);
% toc

if(ShowDetailMessages==1)
    display('finish build kdtree and angle calculation');
end

radii = 40;

[idx, dist] = KDTreeBallQuery(ordered_end_points, ordered_end_points, radii);


%%
% Initialize the graph
% tic
E = [];
W = [];
angle_array=[];
connection_xy_cell=cell(1,1);
%

% for now the parameters are set here, without adaptiveness
% for the loop number,  include everything in the unmatched bw and the
% current model lines( could be connected lines)
num_loop=0;

for i_end = 1 : size(idx,1)
    if(circle_flag(i_end)==0 && ~isempty(ordered_curve_xy{i_end}))
        
        ax_i = ordered_curve_xy{i_end}(:,1);
        ay_i = ordered_curve_xy{i_end}(:,2);
        line_i_x =  ordered_smooth_curve_xy{i_end}(:,1);
        line_i_y =  ordered_smooth_curve_xy{i_end}(:,2);
        
        angle_i = end_angles(i_end,1);
        x_i = ordered_end_points(i_end,1);
        y_i = ordered_end_points(i_end,2);
        
        %   j end cannot be i it self and cannot be smaller, since it i big - j samll will
        %   be calculated on the other pair, which is i small- j big
        
        for j_end =  setdiff(idx{i_end,1}',[1:2*round(i_end/2)])
            %         try
            num_loop = num_loop+1;
            %             display(['loop num: ',num2str(num_loop)]);
            
            ax_j = ordered_curve_xy{j_end}(:,1);
            ay_j = ordered_curve_xy{j_end}(:,2);
            line_j_x =  ordered_smooth_curve_xy{j_end}(:,1);
            line_j_y =  ordered_smooth_curve_xy{j_end}(:,2);
            
            x_j = ordered_end_points(j_end,1);
            y_j = ordered_end_points(j_end,2);
            angle_j = end_angles(j_end,1)+pi;
            
            connect_xy = [line_j_x(1)-line_i_x(1) line_j_y(1)-line_i_y(1)];
            connect_n = connect_xy/norm(connect_xy);
            connect_length = norm(connect_xy);
            
            connect_line_fine = 0:0.2:connect_length;
            
            connect_x = line_i_x(1)+ connect_n(1).*connect_line_fine;
            connect_y = line_i_y(1)+ connect_n(2).*connect_line_fine;
            
            
            
            % flip the i part due to the order is starting the end point
            % toward i, but both connecting and j is leaving i
            
            angle_c = atan2(connect_x(end)-connect_x(1),connect_y(end)-connect_y(1));
            
            distance_ij = sqrt((ax_i(1)-ax_j(1)).^2+(ay_i(1)-ay_j(1)).^2);
            
            angle_ij = angle_i - angle_j;
            if angle_ij >pi
                angle_ij = angle_ij -2*pi;
            end
            if angle_ij < -pi
                angle_ij = angle_ij + 2*pi;
            end
            
            
            dist_T = min( 3, min(length(ax_i),length(ax_j))/5);
            
            % this is to exclude the angles larger than certain value, for now, include everything
            if abs(angle_ij) > pi/(3) && distance_ij >dist_T
                continue;
            end
            
            
            angle_ic = angle_i - angle_c;
            if angle_ic >pi
                angle_ic = angle_ic -2*pi;
            end
            if angle_ic < -pi
                angle_ic = angle_ic + 2*pi;
            end
            if abs(angle_ic) > pi/(3) && distance_ij >dist_T
                continue;
            end
            %
            
            angle_cj = angle_c - angle_j;
            % normalize the angle between -pi to pi
            if angle_cj >pi
                angle_cj = angle_cj -2*pi;
            end
            if angle_cj < -pi
                angle_cj = angle_cj + 2*pi;
            end
            if abs(angle_cj) > pi/(3) && distance_ij >dist_T
                continue;
            end
            %
            %
            
            if distance_ij > min(min(length(ax_i),length(ax_j)),30)
                continue;
            end
            
            
            extend_i_n = [-cos(angle_i+pi/2) sin(angle_i+pi/2)];
            extend_j_n = [-cos(angle_j-pi/2) sin(angle_j-pi/2)];
            line_fine = -(0.5)*distance_ij:0.1:(1.5)*distance_ij;
            
            extend_i_x = ax_i(1)+ extend_i_n(1).*line_fine;
            extend_i_y = ay_i(1)+ extend_i_n(2).*line_fine;
            extend_j_x = ax_j(1)+ extend_j_n(1).*line_fine;
            extend_j_y = ay_j(1)+ extend_j_n(2).*line_fine;
            
            [ind_i_j, distance_i_p_j] = distance_two_curves([ax_i(1) ay_i(1)], [extend_j_x' extend_j_y']);
            [ind_j_i, distance_j_p_i] = distance_two_curves([ax_j(1) ay_j(1)], [extend_i_x' extend_i_y']);
            
            sigma_dis = min(15, min(length(ax_i),length(ax_j))/2);
            
            sigma_ij = pi/6;
            sigma_ic = pi/6;
            sigma_cj = pi/6;
            
            angle_array = [angle_array; angle_ij angle_ic angle_cj];
            
            angle_w_1 = (distance_i_p_j/distance_ij).*(angle_cj.^2)/2/(sigma_cj.^2);
            angle_w_2 = (distance_j_p_i/distance_ij).*(angle_ic.^2)/2/(sigma_ic.^2);
            dis_w = (1/2)*(distance_ij.^2)/2/(sigma_dis.^2);
            
            % the threshold of the weight (over the exp(-...))
            T_w = 2*(1/2).*((pi/3).^2)/2/(sigma_cj.^2) + (1/2)*(20.^2)/2/(sigma_dis.^2) ;
            
            % if the weight is too small, skip it;
            % keep it if the exp(-...) term is smaller than the T_w threshold
            if (angle_w_1 + angle_w_2 + dis_w) < T_w*2
                W = [W; exp(-angle_w_1)...
                    *exp(-angle_w_2)...
                    *exp(-dis_w)...
                    ];
                
                E = [E; i_end j_end];
                
                connection_xy_cell{size(E,1)} = [connect_x;connect_y];
            end
            
            %             h11=figure(11); hold off; plot(ax_i,ay_i,'.');
            %             axis equal
            %             hold on;
            %             plot(ax_j,ay_j,'r.');
            %             plot(ax_j(1),ay_j(1),'r*');
            %             plot(ax_i(1),ay_i(1),'*');
            %
            %             plot(line_i_x,line_i_y,'-');
            %             plot(line_j_x,line_j_y,'r-');
            %
            %             plot(extend_i_x,extend_i_y,':');
            %             plot(extend_j_x,extend_j_y,'r:');
            %
            %             plot([ax_i(1) extend_j_x(ind_i_j(2))],[ay_i(1) extend_j_y(ind_i_j(2))],'m');
            %             plot([ax_j(1) extend_i_x(ind_j_i(2))],[ay_j(1) extend_i_y(ind_j_i(2))],'g');
            %
            %
            %             legend({['angle-ic ',num2str(angle_ic,'%.3f')], ['angle-cj ',num2str(angle_cj,'%.3f')], ...
            %                 ['angle-ij ',num2str(angle_ij,'%.3f')],...
            %                 ['pd i to j ',num2str(distance_i_p_j,'%.3f')],...
            %                 ['pd j to i ',num2str(distance_j_p_i,'%.3f')],...
            %                 ['angel i w ', num2str(angle_w_1,'%.3f')],...
            %                 ['angel j w ', num2str(angle_w_2,'%.3f')],...
            %                 ['dist w ', num2str(dis_w,'%.3f')]...
            %                 });
            %
            %             saveas(h11, ['./endij/end_i_',num2str(i_end),'_end_j_',num2str(j_end),'.tif']);
            
            %                       pause;
            
            
        end
    end
end
if length(E)<4
    return;
end
% toc

if(ShowDetailMessages==1)
    display('finish build graph');
end

M = maxGreedyMatching(E,W);
%%

matched_ind = find(M>0);
Weight_Matched = W(matched_ind);

% find the weight at the confidency interval
weight_ci = prctile(Weight_Matched,(1-confidency_interval)*100);

% the match with weight lower than the ci wants are eliminated
M(matched_ind( find( Weight_Matched < weight_ci) ) ) = 0;

%%

h1=figure(1); set(h1,'Visible',set_visible); hold off; imagescc(currentImg); hold on;

current_model = cell(1,1);
count_linking = 0;
new_current_all_matching_bw = zeros(size(bw_to_be_matched));

for i_E = 1 : length(W)
    if M(i_E)==1
        
        x_i = ordered_end_points(E(i_E,1),1);
        y_i = ordered_end_points(E(i_E,1),2);
        x_j = ordered_end_points(E(i_E,2),1);
        y_j = ordered_end_points(E(i_E,2),2);
        
        ax_i = ordered_curve_xy{E(i_E,1)}(:,1);
        ay_i = ordered_curve_xy{E(i_E,1)}(:,2);
        ax_j = ordered_curve_xy{E(i_E,2)}(:,1);
        ay_j = ordered_curve_xy{E(i_E,2)}(:,2);
        
        
        connect_xy = connection_xy_cell{i_E};
        end_1_ind = E(i_E,1);
        end_2_ind = E(i_E,2);
        
        
        plot([ x_i x_j], [y_i y_j],'y-','LineWidth',3);
        hold on;
        plot(ax_i, ay_i,'b.');
        plot(ax_j, ay_j,'r.');
    end
end


build_open_model =cell(1,1);
copy_lines = ordered_curve_xy;
ind_matched_model = nan(size(ordered_curve_xy,1),1);
iM =0;
copy_M = M;

transformer = [];


model_ind_all(model_ind_all==0)=nan;


if(iIteration==3)
    stophere=1;
end

%%
iM=0;
while(sum(M)>0)
    
    for i_E = 1 : length(M)
        if M(i_E)==1
            %%
            % find all the indices first
            % here for index
            % i_E the pair of end points
            % all the i_area, j_area, i_end, j_end is for index in the tomatched ind
            % all model_ind_all is for original gobal index for the curves
            
            i_end = E(i_E,1);
            j_end = E(i_E,2);
            
            part_i = ordered_curve_xy{i_end};
            part_j = ordered_curve_xy{j_end};
            part_c =  connection_xy_cell{i_E}';
            
            i_area = round(E(i_E,1)/2);
            j_area = round(E(i_E,2)/2);
            
            % put into final rescue set and the connections
            
            try
                % if failed, mean the i_area is not in tomatch_ind, so it is already matched
                i_curve_index = tomatch_ind(i_area);
            catch
                i_curve_index = nan;
            end
            
            try
                % if failed, mean the i_area is not in tomatch_ind, so it is already matched
                j_curve_index = tomatch_ind(j_area);
            catch
                j_curve_index = nan;
            end
            
            final_set_rescue = [final_set_rescue; i_curve_index j_curve_index];
            
            
            connect_count = connect_count+1;
            final_set_connection{connect_count} = connection_xy_cell{i_E};
            
            
            if(mod(i_end,2)==1)
                i_other_end = i_end+1;
            else
                i_other_end = i_end-1;
            end
            
            if(mod(j_end,2)==1)
                j_other_end = j_end+1;
            else
                j_other_end = j_end-1;
            end
            
            
            %%
            % see if the ind has been matched already matched to something else
            % find the index first, assumming it is going to be added
            if (isnan(ind_matched_model(E(i_E,1)))...
                    && isnan(ind_matched_model(E(i_E,2))) )
                iM=iM+1;
                current_iM = iM;
            else
                if (~isnan(ind_matched_model(E(i_E,1))))
                    % this matching, one end has been linked to a matching E(i_E,1)
                    current_iM = ind_matched_model(E(i_E,1));
                else
                    current_iM = ind_matched_model(E(i_E,2));
                end
            end
            %%
            
            % if both of the endpoints come from curvelet in the new pool,
            % and this is not the first time matching
            if ( E(i_E,1)< length(tomatch_ind)*2 && E(i_E,2)< length(tomatch_ind)*2 && ~isempty(orig_current_model) ...
                    && isnan(model_ind_all(i_area,2)) && isnan(model_ind_all(j_area,2)))
                % the last two requirements of the if, means, they are not
                % linked to something else, like master network, or new
                % masters(in case of three new(new in candidate) segments
                % become one filament, then the first pair need to have
                % requirement. If the first pair doesn't qualify for
                % requirement below, but later pairs does, then it doesn't
                % matter, since the missed segments will be picked up in
                % the next iteration
                
                % then, the two curvelet need to become a "master" curvelet
                % to qualify as new one add to current model
                new_meannms = mean([feature_all.feature_MeanNMS(tomatch_ind(i_area)) feature_all.feature_MeanNMS(tomatch_ind(j_area))]);
                new_length = sum([feature_all.feature_Length(tomatch_ind(i_area)) feature_all.feature_Length(tomatch_ind(j_area))]);
                new_meanint = mean([feature_all.feature_MeanInt(tomatch_ind(i_area)) feature_all.feature_MeanInt(tomatch_ind(j_area))]);
                new_curv = sum([feature_all.feature_Curvature(tomatch_ind(i_area)) feature_all.feature_Curvature(tomatch_ind(j_area))]);
                %                 if this qualifies, add it, if not, skip it.
                
                line_xy = [flipud(part_i);part_c;part_j];
                line_i_x = line_xy(:,1);
                line_i_y = line_xy(:,2);
                
                line_i_x = (imfilter(line_i_x, [1 3 1]'/5, 'replicate', 'same'));
                line_i_y = (imfilter(line_i_y, [1 3 1]'/5, 'replicate', 'same'));
                
                Vertices = [line_i_x'; line_i_y']';
                Lines=[(1:size(Vertices,1)-1)' (2:size(Vertices,1))'];
                k=LineCurvature2D(Vertices,Lines);
                new_curv=mean(abs(k));
                
                % if these two qualify to be new feature
                if(master_flassier(new_meannms,new_length,new_curv,new_meanint)>0)
                    
                    % add to the system
                    
                    ind_matched_model(E(i_E,1))=current_iM;
                    ind_matched_model(E(i_E,2))=current_iM;
                    
                    ordered_curve_xy{i_end} = flipud([flipud(part_i);part_c;part_j]);
                    ordered_curve_xy{i_other_end} = ([flipud(part_i);part_c;part_j]);
                    
                    ordered_curve_xy{j_other_end} = flipud([flipud(part_i);part_c;part_j]);
                    ordered_curve_xy{j_end} = ([flipud(part_i);part_c;part_j]);
                    
                    ind_matched_model(i_other_end)=current_iM;
                    ind_matched_model(j_other_end)=current_iM;
                    
                    % add to the model and (big) index info
                    
                    build_open_model{current_iM,1} = [flipud(part_i);part_c;part_j];
                    temp_indices = [model_ind_all(i_area,:) model_ind_all(j_area,:)];
                    temp_indices = temp_indices(find(isnan(temp_indices)==0));
                    temp_indices = temp_indices(find(temp_indices~=0));
                    new_model_ind(current_iM,1:length(temp_indices))= temp_indices;
                    if length(temp_indices)>50
                        stophere=1;
                    end
                    %                     figure(21);imagescc(model_ind_all');axis on;
                    %                     figure(22);bar(temp_indices);title(['i area:',num2str(i_area),', j area:',num2str(j_area)]);
                    
                    
                    model_ind_all(i_area,1:length(temp_indices))= temp_indices;
                    model_ind_all(j_area,1:length(temp_indices))= temp_indices;
                    
                    transformer = [transformer; tomatch_ind(i_area) tomatch_ind(j_area)];
                    
                    % put into final rescue and connecting lines
                    try
                        % if failed, mean the i_area is not in tomatch_ind, so it is already matched
                        i_curve_index = tomatch_ind(i_area);
                    catch
                        i_curve_index = nan;
                    end
                    
                    try
                        % if failed, mean the i_area is not in tomatch_ind, so it is already matched
                        j_curve_index = tomatch_ind(j_area);
                    catch
                        j_curve_index = nan;
                    end
                    
                    final_set_rescue = [final_set_rescue; i_curve_index j_curve_index];
                    
                    connect_count = connect_count+1;
                    final_set_connection{connect_count} = connection_xy_cell{i_E};
                    
                else
                    % if this is not successful, then remove this from the
                    % matched index, since the if requirement made it sure
                    % to be a new iM, so can remove this here.
                    iM=iM-1;
                    current_iM = iM;
                end
            else
                
                % if this is not the case, add this right on
                
                % add to the system
                
                ind_matched_model(E(i_E,1))=current_iM;
                ind_matched_model(E(i_E,2))=current_iM;
                
                ordered_curve_xy{i_end} = flipud([flipud(part_i);part_c;part_j]);
                ordered_curve_xy{i_other_end} = ([flipud(part_i);part_c;part_j]);
                
                ordered_curve_xy{j_other_end} = flipud([flipud(part_i);part_c;part_j]);
                ordered_curve_xy{j_end} = ([flipud(part_i);part_c;part_j]);
                
                ind_matched_model(i_other_end)=current_iM;
                ind_matched_model(j_other_end)=current_iM;
                
                % add to the model and (big) index info
                
                build_open_model{current_iM,1} = [flipud(part_i);part_c;part_j];
                temp_indices = [model_ind_all(i_area,:) model_ind_all(j_area,:)];
                temp_indices = temp_indices(find(isnan(temp_indices)==0));
                temp_indices = temp_indices(find(temp_indices~=0));
                new_model_ind(current_iM,1:length(temp_indices))= temp_indices;
                if length(temp_indices)>50
                    stophere=1;
                end
                %                 figure(21);imagescc(model_ind_all');axis on;
                %                 figure(22);bar(temp_indices);title(['i area:',num2str(i_area),', j area:',num2str(j_area)]);
                
                model_ind_all(i_area,1:length(temp_indices))= temp_indices;
                model_ind_all(j_area,1:length(temp_indices))= temp_indices;
                
                
            end
            
            M(i_E)=0;
            
        end
    end
end

current_model = build_open_model;

if(iIteration==3)
    stophere=1;
end
%%
% check how the new current model is updated from the original model
count_linking=size(build_open_model,1);
current_good_bw = zeros(size(new_current_all_matching_bw));
if(~isempty(orig_current_model))
    for i_E = 1 : length(orig_current_model)
        y_i = orig_current_model{i_E}(:,2);
        x_i = orig_current_model{i_E}(:,1);
        
        % if this good line is not linked to something else
        if(mean(new_current_all_matching_bw(sub2ind(size(new_current_all_matching_bw), round(y_i),round(x_i))))<1)
            
            count_linking = count_linking + 1;
            current_model{count_linking} = [x_i, y_i];
            new_model_ind(count_linking,1:length(model_ind(i_E,:)))= model_ind(i_E,:);
        end
    end
else
    %%
    % in the case when original model is empty, we know this is the first
    % matching with the input bw containing only the good curves, so add
    % the rest of the unmatched curve to the model.
    matched_ind = E(find(copy_M>0),:);
    
    for i_area = 1 : length(tomatch_ind)
        
        ind = find(matched_ind==2*i_area-1 | matched_ind==2*i_area);
        if(isempty(ind) && ~isempty(ordered_curve_xy{2*i_area-1,1}))
            count_linking = count_linking+1;
            current_model{count_linking} = ordered_curve_xy{2*i_area-1,1};
            new_model_ind(count_linking,1)= tomatch_ind(i_area);
        end
    end
end

new_model_ind(new_model_ind==0)=nan;

%% plot the output image with these good ones
current_matching_bw = zeros(size(labelMask));

% if there is already a current model, build the image for matched lines
if(~isempty(current_model))
    for i_E = 1 : length(current_model)
        y_i = current_model{i_E}(:,2);
        x_i = current_model{i_E}(:,1);
        ind = find(~isnan(x_i));
        x_i = x_i(ind);
        y_i = y_i(ind);
        current_model{i_E} = [x_i y_i];
        
        if(~isempty(x_i))
            if(i_E>size(build_open_model,1) )
                current_matching_bw(sub2ind(size(current_matching_bw), round(y_i),round(x_i)))=1.5;
            else
                current_matching_bw(sub2ind(size(current_matching_bw), round(y_i),round(x_i)))=2;
            end
        end
    end
end

h2=figure(2);set(h2,'Visible',set_visible);imagescc(current_matching_bw);

output_model = current_model;

h3=figure(3);set(h3,'Visible',set_visible);imagescc(max(max(currentImg))-currentImg);hold on;
[y,x]=find(current_matching_bw>0);

plot(x,y,'.','MarkerSize',3);

% h1 = figure(2);print(h1,'-dtiff','new_all.tif');

end_point_stop =1 ;

bw_rgb_three_color = nan(size(currentImg,1),size(currentImg,2),3);
R = currentImg*0;
G = currentImg*0;
B = currentImg*0;

bw_whitergb_three_color = nan(size(currentImg,1),size(currentImg,2),3);
whiteR = currentImg*0+1;
whiteG = currentImg*0+1;
whiteB = currentImg*0+1;

h11 = figure(11);
set(h11,'Visible',set_visible);imagescc(currentImg*0); caxis([0 1])

h21 = figure(21);
set(h21,'Visible',set_visible);imagescc(currentImg*0+1); caxis([0 1])

% Original_set_Good_ind=Good_ind;
for ind = Original_set_Good_ind'
    [ind_y,ind_x] = find(labelMask==ind);
    h11=figure(11);set(h11,'Visible',set_visible);
    hold on;
    plot(ind_x,ind_y,'.');
    
    R(sub2ind(size(currentImg),round(ind_y),round(ind_x)))=0;
    G(sub2ind(size(currentImg),round(ind_y),round(ind_x)))=0;
    B(sub2ind(size(currentImg),round(ind_y),round(ind_x)))=1;
end


for ind = Original_set_Good_ind'
    [ind_y,ind_x] = find(labelMask==ind);
    h21=figure(21);set(h21,'Visible',set_visible);
    hold on;
    plot(ind_x,ind_y,'.');
    
    whiteR(sub2ind(size(currentImg),round(ind_y),round(ind_x)))=0;
    whiteG(sub2ind(size(currentImg),round(ind_y),round(ind_x)))=0;
    whiteB(sub2ind(size(currentImg),round(ind_y),round(ind_x)))=1;
end

for ind = setdiff(final_set_rescue(:),Original_set_Good_ind)'
    [ind_y,ind_x] = find(labelMask==ind);
    h11=figure(11);hold on;set(h11,'Visible',set_visible);
    plot(ind_x,ind_y,'r.');
    
    R(sub2ind(size(currentImg),round(ind_y),round(ind_x)))=1;
    G(sub2ind(size(currentImg),round(ind_y),round(ind_x)))=0;
    B(sub2ind(size(currentImg),round(ind_y),round(ind_x)))=0;
    
end

for ind = setdiff(final_set_rescue(:),Original_set_Good_ind)'
    [ind_y,ind_x] = find(labelMask==ind);
    h21=figure(21);hold on;set(h21,'Visible',set_visible);
    hold on;
    plot(ind_x,ind_y,'r.');
    
    whiteR(sub2ind(size(currentImg),round(ind_y),round(ind_x)))=1;
    whiteG(sub2ind(size(currentImg),round(ind_y),round(ind_x)))=0;
    whiteB(sub2ind(size(currentImg),round(ind_y),round(ind_x)))=0;
    
end

for ind = 1:connect_count
    XY = final_set_connection{ind};
    h11=figure(11);hold on;set(h11,'Visible',set_visible);
    hold on;
    plot(XY(1,:),XY(2,:),'g.');
    
    R(sub2ind(size(currentImg),round(XY(2,:)),round(XY(1,:))))=0;
    G(sub2ind(size(currentImg),round(XY(2,:)),round(XY(1,:))))=0.7;
    B(sub2ind(size(currentImg),round(XY(2,:)),round(XY(1,:))))=0;
    
    
end

for ind = 1:connect_count
    XY = final_set_connection{ind};
    h21=figure(21);hold on;set(h21,'Visible',set_visible);
    hold on;
    plot(XY(1,:),XY(2,:),'g.');
    
    whiteR(sub2ind(size(currentImg),round(XY(2,:)),round(XY(1,:))))=0;
    whiteG(sub2ind(size(currentImg),round(XY(2,:)),round(XY(1,:))))=0.7;
    whiteB(sub2ind(size(currentImg),round(XY(2,:)),round(XY(1,:))))=0;
    
    
end

bw_rgb_three_color(:,:,1)=R;
bw_rgb_three_color(:,:,2)=G;
bw_rgb_three_color(:,:,3)=B;
bw_whitergb_three_color(:,:,1)=whiteR;
bw_whitergb_three_color(:,:,2)=whiteG;
bw_whitergb_three_color(:,:,3)=whiteB;

