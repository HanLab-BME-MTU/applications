function [digital_model,orientation_model,XX,YY,OO, II] ...
    = filament_model_to_digital_with_orientation(current_model)

model_length = length(current_model);

digital_model = cell(1,model_length);
orientation_model = cell(1,model_length);

line_smooth_H = fspecial('gaussian',5,2);
% initialize the array output, II is for index
XX=[];YY=[];OO=[];II=[];

for iFila = 1 : model_length
    try
    line_i_x = current_model{iFila}(:,1);
    line_i_y = current_model{iFila}(:,2);
    
    line_i_x = (imfilter(line_i_x, line_smooth_H, 'replicate', 'same'));
    line_i_y = (imfilter(line_i_y, line_smooth_H, 'replicate', 'same'));
    if(length(line_i_x)==1)
        continue;
    end
    
    if(length(line_i_x)==2)
     angles = atan2(-line_i_x(1)+line_i_x(2),-line_i_y(1)+line_i_y(2));
    else
    angles = atan2(-line_i_x(1:end-2)+line_i_x(3:end),-line_i_y(1:end-2)+line_i_y(3:end));
    end
    
    angles = [angles(1); angles; angles(end)];
    
    digital_x = round(line_i_x);
    digital_y = round(line_i_y);
    
    sort_x = digital_x(1);
    sort_y = digital_y(1);
    sort_ang = [];
    angle_pool = angles(1);
    
    for j = 2 : length(digital_x)
        if digital_x(j)==sort_x(end) && digital_y(j)==sort_y(end)
            %same point as previous, put the angle to the pool
            angle_pool = [angle_pool; angles(j)];
        else
            % new digital point? add to the sorted x and put the angle of
            % last group to the 
            sort_x = [sort_x;digital_x(j)];
            sort_y = [sort_y;digital_y(j)];     
            sort_ang = [sort_ang; mean(angle_pool)];
            angle_pool = angles(j);
        end
    end
     sort_ang = [sort_ang; mean(angle_pool)];
    digital_model{iFila} = [sort_x sort_y];
    orientation_model{iFila} = sort_ang;
    
    II = [II;iFila*(sort_x*0+1)];
    XX = [XX;sort_x;];
    YY = [YY;sort_y;];
    OO = [OO;sort_ang;];
    end
end