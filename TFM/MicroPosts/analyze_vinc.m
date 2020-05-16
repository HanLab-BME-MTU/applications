warning off;

%vinc4a=imrotate(vinc4,rot_angle);
%vinc4b=imcrop(vinc4a,[p1 offset]);

% vinc2a=imrotate(vinc2,rot_angle);
% vinc2b=imcrop(vinc2a,[p1 offset]);

comp4v(:,:,1)=vinc2b(:,:);

comp4v(:,:,2)=f3b_top(:,:);

comp4v(:,:,3)=f3b_bot(:,:);

[boundaries_vinc,L_vinc]=bwboundaries(vinc4b,4,'noholes');
    
s_vinc=regionprops(L_vinc,'basic');

k_v=1;

for i=1:length(s_vinc);
    if(s_vinc(i).Area*convert^2 < 100);
        b=boundaries_vinc{i};
            if length(b)>(10);
%               post_vinc_h(i)=plot(b(:,2),b(:,1),'y');
                b_vinc(k_v)=i;
            k_v=k_v+1;
            end
    end
end


for i=1:row;
    for j=1:col;
        vinc2_area(i,j)=0;
        vinc2_Orientation(i,j)=0;
        vinc2_MajorAxisLength(i,j)=0;
        vinc2_MinorAxisLength(i,j)=0;
    end
end

roi_vinc=zeros(size(y2_top,1),size(y2_top,2));

for i=1:k_v-1;
    roi_vinc=roi_vinc+poly2mask(boundaries_vinc{b_vinc(i)}(:,2),boundaries_vinc{b_vinc(i)}(:,1),size(y2_top,1),size(y2_top,2));
end

box_size=80;

for l=1:row;
    for j=1:col;
            
            if and(degen_posts_x(l,j)==0,degen_posts_y(l,j)==0);
            
                if actual_x_top(l,j)-box_size/2 < 1;
                    vinc_lim1=1;
                else
                    vinc_lim1=round(actual_x_top(l,j)-box_size/2);
                end
                
                if actual_x_top(l,j)+box_size/2 > size(y2_top,1);
                    vinc_lim2=size(y2_top,1);
                else
                    vinc_lim2=round(actual_x_top(l,j)+box_size/2);
                end
                
                if actual_y_top(l,j)-box_size/2 < 1;
                    vinc_lim3=1;
                else
                    vinc_lim3=round(actual_y_top(l,j)-box_size/2);
                end
                
                if actual_y_top(l,j)+box_size/2 > size(y2_top,2);
                    vinc_lim4=size(y2_top,2);
                else
                    vinc_lim4=round(actual_y_top(l,j)+box_size/2);
                end
                
                box_ratio(l,j)=(vinc_lim2-vinc_lim1)*(vinc_lim4-vinc_lim3)/6400;
                
                A_vinc = roi_vinc(vinc_lim1:vinc_lim2,vinc_lim3:vinc_lim4);
                A_post=I7(vinc_lim1:vinc_lim2,vinc_lim3:vinc_lim4);
 %               vinc_corr = A_vinc*double(double(A_post));
                
                [boundaries_vinc2,L_vinc2]=bwboundaries(A_vinc,4,'noholes');
                s_vinc2=regionprops(L_vinc2,'Area','Orientation','MajorAxisLength','MinorAxisLength');
                vinc2_a=cat(1,s_vinc2.Area);
                vinc2_o=cat(1,s_vinc2.Orientation);
                vinc2_m1=cat(1,s_vinc2.MajorAxisLength);
                vinc2_m2=cat(1,s_vinc2.MinorAxisLength);
                
%                if vinc2_a < (80*80*0.50)
                    vinc2_area(l,j)=vinc2_area(l,j)+sum(vinc2_a);
%                end

                 if isnan(mean(vinc2_o));
                 else
                     vinc2_Orientation(l,j)=-mean(vinc2_o);
                 end
%                 
%                 if isnan(mean(vinc2_m1));
%                 else    
%                     vinc2_MajorAxisLength(l,j)=vinc2_MajorAxisLength(l,j)+mean(vinc2_m1);
%                 end
%                 
%                 if isnan(mean(vinc2_m2));                
%                 else
%                     vinc2_MinorAxisLength(l,j)=vinc2_MinorAxisLength(l,j)+mean(vinc2_m2);
%                 end
                
            end
            
        end 
    end

sum_vinc=0;
sum_ext_vinc=0;
unocc_post_vinc=0;

for i=1:row;
    for j=1:col;
        if and(cell_posts_x(i,j)==0,cell_posts_y(i,j)==0);
            unocc_post_vinc=unocc_post_vinc+vinc2_area(i,j);
        end
    end
end

avg_unocc_vinc=unocc_post_vinc/border_num;

for i=1:row;
    for j=1:col;
        
        if vinc2_area(i,j)-avg_unocc_vinc > 0;
            vinc2_area(i,j)=vinc2_area(i,j)-avg_unocc_vinc;
        else
            vinc2_area(i,j)=0;
        end
        
         if and(cell_posts_mag(i,j) ~= 0,vinc2_area(i,j) ~=0);                
             sum_vinc=sum_vinc+vinc2_area(i,j)/box_ratio(i,j);          
         end
        
        if exterior_posts_mag(i,j) ~=0;
            sum_ext_vinc=sum_ext_vinc+vinc2_area(i,j)/box_ratio(i,j);
        end
        
        
    end
end



average_vinc=sum_vinc/cell_num;
sum_int_vinc=sum_vinc-sum_ext_vinc;


k1=1;
k2=1;
k3=1;

force_exterior_dir_list_v=[];
vinc_exterior_dir_list=[];
force_interior_dir_list_v=[];
vinc_interior_dir_list=[];
force_unocc_dir_list_v=[];
vinc_unocc_dir_list=[];

force_exterior_list=[];
vinc_exterior_list=[];
force_interior_list=[];
vinc_interior_list=[];
force_unocc_list=[];
vinc_unocc_list=[];

for i=1:row;
    for j=1:col;
%         vinc2_x(i,j)=vinc2_MajorAxisLength(i,j)*-sin(vinc2_Orientation(i,j)/180*pi);
%         vinc2_y(i,j)=vinc2_MajorAxisLength(i,j)*cos(vinc2_Orientation(i,j)/180*pi);
%         vinc_corr(i,j)=vinc2_area(i,j)*force_mag(i,j)/(average_vinc*avg_post_force);
        
        if exterior_posts_mag(i,j) ~= 0;
            vinc_exterior_list(k1)=vinc2_area(i,j)*convert^2;
            force_exterior_list(k1)=force_mag(i,j);
            if and(vinc2_area(i,j)*convert^2 > 2.5,or(force_mag(i,j) > 5, border_posts_mag(i,j) > 5));
                vinc_exterior_dir_list(k1)=vinc2_Orientation(i,j);
                if and(force_dir(i,j) <= 90, force_dir(i,j) >= -90);
                    force_exterior_dir_list_v(k1)=force_dir(i,j);
                elseif force_dir(i,j) > 90;
                    force_exterior_dir_list_v(k1)=force_dir(i,j)-180;
                elseif force_dir(i,j) < -90;
                    force_exterior_dir_list_v(k1)=force_dir(i,j)+180;
                end
            end
            k1=k1+1;
        elseif cell_posts_mag(i,j) ~=0;
            vinc_interior_list(k2)=vinc2_area(i,j)*convert^2;
            force_interior_list(k2)=force_mag(i,j);
            if and(vinc2_area(i,j)*convert^2 > 2.5,or(force_mag(i,j) > 5, border_posts_mag(i,j) > 5));
                vinc_interior_dir_list(k2)=vinc2_Orientation(i,j);
                if and(force_dir(i,j) <= 90, force_dir(i,j) >= -90);
                    force_interior_dir_list_v(k2)=force_dir(i,j);
                elseif force_dir(i,j) > 90;
                    force_interior_dir_list_v(k2)=force_dir(i,j)-180;
                elseif force_dir(i,j) < -90;
                    force_interior_dir_list_v(k2)=force_dir(i,j)+180;
                end
            end
            k2=k2+1;
        else
            vinc_unocc_list(k3)=vinc2_area(i,j)*convert^2;
            force_unocc_list(k3)=force_unocc_mag(i,j);
            if and(vinc2_area(i,j)*convert^2 > 2.5,or(force_mag(i,j) > 5, border_posts_mag(i,j) > 5));
                vinc_unocc_dir_list(k3)=vinc2_Orientation(i,j);
                if and(border_posts_dir(i,j) <= 90, border_posts_dir(i,j) >= -90);
                    force_unocc_dir_list_v(k3)=border_posts_dir(i,j);
                elseif border_posts_dir(i,j) > 90;
                    force_unocc_dir_list_v(k3)=border_posts_dir(i,j)-180;
                elseif border_posts_dir(i,j) < -90;
                    force_unocc_dir_list_v(k3)=border_posts_dir(i,j)+180;
                end        
            end
            k3=k3+1;
        end
    end
end



