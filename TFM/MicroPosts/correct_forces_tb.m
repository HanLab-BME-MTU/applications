% CAL 2/5/03
%
% Correct forces by subtracting average background noise

step_factor=1;

% Subtract background noise from post deflections

delta_x_micron=delta_x_micron-avg_grid_err_x/step_factor;
delta_y_micron=delta_y_micron-avg_grid_err_y/step_factor;
 
delta_x=delta_x_micron/convert;
delta_y=delta_y_micron/convert;


for i=1:row;
    for j=1:col;
        
        if or(border_posts_x(i,j) ~= 0, border_posts_y(i,j) ~= 0);        
            border_posts_x(i,j) = border_posts_x(i,j)-avg_grid_err_x;
            border_posts_y(i,j) = border_posts_y(i,j)-avg_grid_err_y;
        end
        
        if or(cell_posts_x(i,j) ~=0, cell_posts_y(i,j) ~=0);
            cell_posts_x(i,j)=cell_posts_x(i,j)-avg_grid_err_x;
            cell_posts_y(i,j)=cell_posts_y(i,j)-avg_grid_err_y;
        end
        
        if or(exterior_posts_x(i,j) ~=0, exterior_posts_y(i,j) ~=0);
            exterior_posts_x(i,j)=exterior_posts_x(i,j)-avg_grid_err_x;
            exterior_posts_y(i,j)=exterior_posts_y(i,j)-avg_grid_err_y;
        end
        
        if or(interior_posts_x(i,j) ~=0, interior_posts_y(i,j) ~=0);
            interior_posts_x(i,j)=interior_posts_x(i,j)-avg_grid_err_x;
            interior_posts_y(i,j)=interior_posts_y(i,j)-avg_grid_err_y;
        end
    end
end

force_x=cell_posts_x*force_convert;
force_y=cell_posts_y*force_convert;

force_ext_x=exterior_posts_x*force_convert;
force_ext_y=exterior_posts_y*force_convert;

force_int_x=interior_posts_x*force_convert;
force_int_y=interior_posts_y*force_convert;

force_unocc_x=border_posts_x*force_convert;
force_unocc_y=border_posts_y*force_convert;

% 
% % Re-establish post separation and conversions
% 
% %[cell_posts_x,cell_posts_y]=compare_clicks2(cell_x,cell_y,actual_x_bot,actual_y_bot,delta_x,delta_y);
% 
% [degen_posts_x,degen_posts_y]=compare_clicks2(degen_x,degen_y,actual_x_bot,actual_y_bot,delta_x,delta_y);
% 
% cell_num=0;
% 
% for i=1:row;
%     for j=1:col;
%         cell_posts_x(i,j)=0;
%         cell_posts_y(i,j)=0;
%     end
% end
% 
% for i=1:k-1;
%     roi_actin=poly2mask(boundaries_actin{b_actin(i)}(:,2),boundaries_actin{b_actin(i)}(:,1),size(y2_top,1),size(y2_top,2));
%     l=0;
%     for M = 1:80:(hw(1)*80+1-80);
%         l = l + 1;
%         j = 0;
%         for N = 1:80:(hw(2)*80+1-80);
%             j = j+1;
%             
%             if and(degen_posts_x(l,j)==0,degen_posts_y(l,j)==0);
%             
%                 A_act = roi_actin(M:M+79,N:N+79);
%             
%                 act_corr = immultiply(A_act,I6a(:,:,l,j));
%             
%                 if sum(sum(act_corr)) ~= 0;
%                     cell_posts_x(l,j)=delta_x(l,j);
%                     cell_posts_y(l,j)=delta_y(l,j);
%                     cell_num=cell_num+1;
%                     %else
%                     %cell_posts_x(j,l)=0;
%                     %cell_posts_y(j,l)=0;
%                 end
%             end
%         end 
%     end
% end
% 
% for i=1:row;
%     for j=1:col;
%         border_posts_x(i,j)=0;
%         border_posts_y(i,j)=0;
%     end
% end
% 
% for i=1:row;
%     for j=1:col;
%         if and(and(cell_posts_x(i,j)==0,degen_posts_x(i,j)==0),and(cell_posts_y(i,j)==0,degen_posts_x(i,j)==0));
%             border_posts_x(i,j)=delta_x(i,j);
%             border_posts_y(i,j)=delta_y(i,j);
%         end
%     end
% end

