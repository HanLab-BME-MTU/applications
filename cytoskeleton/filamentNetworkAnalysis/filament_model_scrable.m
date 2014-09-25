function [scrable_digital_model,scrable_orientation_model,scrable_VIF_current_seg,scrable_VIF_current_orientation,...
    scrable_XX,scrable_YY,scrable_OO, scrable_II] ...
    = filament_model_scrable(current_model,img_size,T_sigma,O_sigma)


if(nargin<2)
    img_size = [];
end

if(nargin<3)
    T_sigma = 10;
end

if(nargin<4)
    O_sigma = pi/2;
end

model_length = length(current_model);

scrable_digital_model = cell(1,model_length);
scrable_orientation_model = cell(1,model_length);

line_smooth_H = fspecial('gaussian',5,2);
% initialize the array output, II is for index
scrable_XX=[];scrable_YY=[];scrable_OO=[];scrable_II=[];

for iFila = 1 : model_length
    try
    line_i_x = current_model{iFila}(:,1);
    line_i_y = current_model{iFila}(:,2);
    
    mean_line_i_x = mean(line_i_x);
    mean_line_i_y = mean(line_i_y);
       
    translation_rand = randn(2,1)*T_sigma;
    orientation_rand = randn(1,1)*O_sigma;
    
    scrable_line_i_x = (line_i_x-mean_line_i_x)*cos(orientation_rand) ...
        + (line_i_y-mean_line_i_y)*sin(orientation_rand) ...
        + mean_line_i_x +translation_rand(1);
    scrable_line_i_y = -(line_i_x-mean_line_i_x)*sin(orientation_rand) ...
        + (line_i_y-mean_line_i_y)*cos(orientation_rand) ...
        + mean_line_i_y +translation_rand(2);
    
    
    
    line_i_x = (imfilter(scrable_line_i_x, line_smooth_H, 'replicate', 'same'));
    line_i_y = (imfilter(scrable_line_i_y, line_smooth_H, 'replicate', 'same'));
  
    if(length(line_i_x)==1)
        continue;
    end
    
    if(length(line_i_x)==2)
     angles = atan2(-line_i_x(1)+line_i_x(2),-line_i_y(1)+line_i_y(2));
    else
    angles = atan2(-line_i_x(1:end-2)+line_i_x(3:end),-line_i_y(1:end-2)+line_i_y(3:end));
    end
    
    angles = [angles(1); angles; angles(end)];
    
    digital_x = round(scrable_line_i_x);
    digital_y = round(scrable_line_i_y);
    
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
    scrable_digital_model{iFila} = [sort_x sort_y];
    scrable_orientation_model{iFila} = sort_ang;
    
    scrable_II = [scrable_II;iFila*(sort_x*0+1)];
    scrable_XX = [scrable_XX;sort_x;];
    scrable_YY = [scrable_YY;sort_y;];
    scrable_OO = [scrable_OO;sort_ang;];
    end
end

if(isempty(img_size))
img_size(1) = round(max(scrable_YY));
img_size(2) = round(max(scrable_XX));
end

scrable_OO(scrable_OO>pi/2) = scrable_OO(scrable_OO>pi/2) -pi;
scrable_OO(scrable_OO>pi/2) = scrable_OO(scrable_OO>pi/2) -pi;
scrable_OO(scrable_OO>pi/2) = scrable_OO(scrable_OO>pi/2) -pi;

scrable_OO(scrable_OO<-pi/2) = scrable_OO(scrable_OO<-pi/2) +pi;
scrable_OO(scrable_OO<-pi/2) = scrable_OO(scrable_OO<-pi/2) +pi;
scrable_OO(scrable_OO<-pi/2) = scrable_OO(scrable_OO<-pi/2) +pi;


scrable_VIF_current_seg = zeros(img_size);
scrable_VIF_current_orientation = nan(img_size);

new_scrable_XX = scrable_XX(find(scrable_XX>0 & scrable_YY>0 & scrable_XX<=img_size(2) & scrable_YY<=img_size(1)));
new_scrable_YY = scrable_YY(find(scrable_XX>0 & scrable_YY>0 & scrable_XX<=img_size(2) & scrable_YY<=img_size(1)));
new_scrable_OO = scrable_OO(find(scrable_XX>0 & scrable_YY>0 & scrable_XX<=img_size(2) & scrable_YY<=img_size(1)));

scrable_VIF_current_seg(sub2ind(img_size, round(new_scrable_YY),round(new_scrable_XX)))=1; 

scrable_VIF_current_orientation(sub2ind(img_size, round(new_scrable_YY),round(new_scrable_XX)))=new_scrable_OO; 

