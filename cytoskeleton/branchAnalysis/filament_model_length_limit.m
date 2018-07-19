function [length_limited_model,length_limited_seg] = filament_model_length_limit(current_model,current_seg,min_filament_length)

length_limited_model=[];
length_limited_seg = nan(size(current_seg));

model_length = length(current_model);
new_iFila=0;
XX=[];YY=[];

for iFila = 1 : model_length
    try
        line_i_x = current_model{iFila}(:,1);
        line_i_y = current_model{iFila}(:,2);
        
        if(length(line_i_x)<min_filament_length)
            continue;
        end
        
        new_iFila = new_iFila+1;
        
        length_limited_model{new_iFila} = current_model{iFila};        
        
        XX = [XX;line_i_x;];
        YY = [YY;line_i_y;];
    end
end

length_limited_seg(sub2ind(size(current_seg), round(YY),round(XX)))=1; 



