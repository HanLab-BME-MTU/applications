function [scrable_digital_model,scrable_orientation_model,scrable_VIF_current_seg,scrable_VIF_current_orientation,...
    scrable_XX,scrable_YY,scrable_OO, scrable_II] ...
    = filament_model_scrable(current_model,img_size,T_sigma,O_sigma,CellROI)


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

%some new part due to CellROI
count_new_part = model_length;
i_iterate_array=[];

for iFila = 1 : model_length
%     try
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
        
        in_ROI_array=[];
        
        for i_point = 1 : numel(digital_x)
            try
                 in_ROI_array(i_point) = CellROI(digital_y, digital_x);
            catch
                 in_ROI_array(i_point) = 0;      
            end
        end
        
        
        if(sum(1-in_ROI_array))>0
            [L_out,numComponents_out]=bwlabel(1-in_ROI_array);
            [L_in,numComponents_in]=bwlabel(in_ROI_array);
            
            for i_in = 1 : numComponents_in
                ind_a = min(find(L_in==i_in));
                ind_b = max(find(L_in==i_in));
                digital_x_section = round(scrable_line_i_x(ind_a:ind_b));
                digital_y_section = round(scrable_line_i_y(ind_a:ind_b));
                angles_section = angles(ind_a:ind_b);
                
                if i_in ==1
                    iFila_section = iFila;
                else
                    count_new_part = count_new_part+1;
                    iFila_section = count_new_part;
                end
                
                % put this one in
                scrabled_one_fila;
            end
            
             for i_out = 1 : numComponents_out
                ind_a = min(find(L_out==i_out));
                ind_b = max(find(L_out==i_out));
                digital_x_section = round(scrable_line_i_x(ind_a:ind_b));
                digital_y_section = round(scrable_line_i_y(ind_a:ind_b));
                angles_section = angles(ind_a:ind_b);
                
                for i_iterate = 1 :10000
                    try
                     translation_rand_A = randn(2,1)*T_sigma;
                     digital_x_section_new = round(digital_x_section +translation_rand_A(1));
                     digital_y_section_new = round(digital_y_section +translation_rand_A(1));
                     
                     new_in_ROI_array = CellROI(sub2ind(img_size, digital_y_section_new, digital_x_section_new));
                     if(sum(1-new_in_ROI_array))==0
                          digital_x_section = digital_x_section_new;
                          digital_y_section = digital_y_section_new;  
                          
                         break;
                     end
                    end
                end     
                
                if(i_iterate==10000)
                    
                    
                    digital_x_section = round(scrable_line_i_x(ind_a:round((ind_b-ind_a)/2+ind_a)));
                    digital_y_section = round(scrable_line_i_y(ind_a:round((ind_b-ind_a)/2+ind_a)));
                    angles_section = angles(ind_a:round((ind_b-ind_a)/2+ind_a));
                    
                    for i_iterate = 1 :10000
                        try
                            translation_rand_A = randn(2,1)*T_sigma;
                            digital_x_section_new = round(digital_x_section +translation_rand_A(1));
                            digital_y_section_new = round(digital_y_section +translation_rand_A(1));
                            
                            new_in_ROI_array = CellROI(sub2ind(img_size, digital_y_section_new, digital_x_section_new));
                            if(sum(1-new_in_ROI_array))==0
                                digital_x_section = digital_x_section_new;
                                digital_y_section = digital_y_section_new;
                                
                                break;
                            end
                        end
                    end
                    
                
                
                    i_iterate_array(count_new_part)=i_iterate;
                    
                    count_new_part = count_new_part+1;
                    iFila_section = count_new_part;
                    
                    % put this one in
                    scrabled_one_fila;
                    
                    
                    
                    
                    digital_x_section = round(scrable_line_i_x(round((ind_b-ind_a)/2+ind_a))+1:ind_b);
                    digital_y_section = round(scrable_line_i_y(round((ind_b-ind_a)/2+ind_a))+1:ind_b);
                    angles_section = angles(round((ind_b-ind_a)/2+ind_a)+1:ind_b);
                    
                    for i_iterate = 1 :10000
                        try
                            translation_rand_A = randn(2,1)*T_sigma;
                            digital_x_section_new = round(digital_x_section +translation_rand_A(1));
                            digital_y_section_new = round(digital_y_section +translation_rand_A(1));
                            
                            new_in_ROI_array = CellROI(sub2ind(img_size, digital_y_section_new, digital_x_section_new));
                            if(sum(1-new_in_ROI_array))==0
                                digital_x_section = digital_x_section_new;
                                digital_y_section = digital_y_section_new;
                                
                                break;
                            end
                        end
                    end
                    
                
                
                i_iterate_array(count_new_part)=i_iterate;
                
                count_new_part = count_new_part+1;
                iFila_section = count_new_part;                
                
                % put this one in
                scrabled_one_fila;
                
                
                else
                    
                count_new_part = count_new_part+1;
                iFila_section = count_new_part;                
                
                % put this one in
                scrabled_one_fila;  
        end
            
        
            ind_a = 1;
            ind_b = numel(scrable_line_i_x);
            digital_x_section = round(scrable_line_i_x(ind_a:ind_b));
            digital_y_section = round(scrable_line_i_y(ind_a:ind_b));
            angles_section = angles(ind_a:ind_b);
            iFila_section = iFila;
          
            % put this one in
            scrabled_one_fila;
        end
        
%     end
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

