function  [direction_array,angle_array,path, time_setps, profile]  = track_angles_path(x_tr, y_tr, SteerableOutputDir)

nFrame = length(DirectionMap);

[inda,indb]=find(~isnan(x_tr));
start_ind = min(indb);
end_ind = max(indb);

x_tr = x_tr(start_ind:end_ind);
y_tr = y_tr(start_ind:end_ind);

track_length = length(x_tr);

time_setps = start_ind : 0.05:end_ind;

x_tr(find(isnan(x_tr))) = 0;
path_x = interp1(start_ind : 1:end_ind, x_tr,time_setps,'pchip');
y_tr(find(isnan(y_tr))) = 0;
path_y = interp1(start_ind : 1:end_ind, y_tr,time_setps,'pchip');

path = [path_x;path_y];

direction_array = [path_x(2:end) - path_x(1:end-1); ...
    path_y(2:end) - path_y(1:end-1);];

angle_array = atan2(direction_array(1,:),direction_array(2,:));

for i = 2 : length(angle_array)
    diff = angle_array(i)-angle_array(i-1);
    
    if(abs(diff)>5.8)
        diff_pis = round(diff/pi)*pi;
        angle_array(i:end) = angle_array(i:end) - diff_pis;
    end    
end
   
profile = [];
path = round(path);
profile(1,:) = time_setps;

previous_point_frame = 0;
for i = 1 : size(direction_array,2)
   currentFrame = round(time_setps(i));
   if(currentFrame~=previous_point_frame)
    load([SteerableOutputDir,'/steerable_',num2str(currentFrame),'.mat'],'currentImg','ST_map','MAX_res','current_seg','Intensity_Segment','SteerabelRes_Segment');       
   end
   
   VIF_direction = [ cos(ST_map(path(1,i),path(2,i))); sin(ST_map(path(1,i),path(2,i)))];
   MT_direction = direction_array(:,i);
   MT_direction = MT_direction/norm(MT_direction);
   
   profile(2,i) = acos(dot(VIF_direction,MT_direction));
   profile(3,i) =  current_seg(path(1,i),path(2,i));
   profile(4,i) =  MAX_res(path(1,i),path(2,i));
   profile(5,i) =  currentImg(path(1,i),path(2,i));
   
   previous_point_frame = currentFrame;
end
