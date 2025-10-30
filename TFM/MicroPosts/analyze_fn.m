warning off;

%fib4a=imrotate(fib4,rot_angle);
%fib4b=imcrop(fib4a,[p1 offset]);

% fib2a=imrotate(fib2,rot_angle);
% fib2b=imcrop(fib2a,[p1 offset]);

comp4f(:,:,1)=fib2b(:,:);

comp4f(:,:,2)=f3b_top(:,:);

comp4f(:,:,3)=f3b_bot(:,:);

[boundaries_fib,L_fib]=bwboundaries(fib4b,4,'noholes');
    
s_fib=regionprops(L_fib,'basic');

k_f=1;

for i=1:length(s_fib);
    if(s_fib(i).Area*convert^2 < 100);
        b=boundaries_fib{i};
            if length(b)>(10);
%               post_fib_h(i)=plot(b(:,2),b(:,1),'y');
                b_fib(k_f)=i;
            k_f=k_f+1;
            end
    end
end


for i=1:row;
    for j=1:col;
        fib2_area(i,j)=0;
        fib2_Orientation(i,j)=0;
        fib2_MajorAxisLength(i,j)=0;
        fib2_MinorAxisLength(i,j)=0;
    end
end

roi_fn=zeros(size(y2_top,1),size(y2_top,2));

for i=1:k_f-1;
    roi_fn=roi_fn+poly2mask(boundaries_fib{b_fib(i)}(:,2),boundaries_fib{b_fib(i)}(:,1),size(y2_top,1),size(y2_top,2));
end

box_size=80;

for l=1:row;
    for j=1:col;
            
            if and(degen_posts_x(l,j)==0,degen_posts_y(l,j)==0);
            
                if actual_x_top(l,j)-box_size/2 < 1;
                    fib_lim1=1;
                else
                    fib_lim1=round(actual_x_top(l,j)-box_size/2);
                end
                
                if actual_x_top(l,j)+box_size/2 > size(y2_top,1);
                    fib_lim2=size(y2_top,1);
                else
                    fib_lim2=round(actual_x_top(l,j)+box_size/2);
                end
                
                if actual_y_top(l,j)-box_size/2 < 1;
                    fib_lim3=1;
                else
                    fib_lim3=round(actual_y_top(l,j)-box_size/2);
                end
                
                if actual_y_top(l,j)+box_size/2 > size(y2_top,2);
                    fib_lim4=size(y2_top,2);
                else
                    fib_lim4=round(actual_y_top(l,j)+box_size/2);
                end
                
                box_ratio(l,j)=(fib_lim2-fib_lim1)*(fib_lim4-fib_lim3)/6400;
                
                A_fib = roi_fn(fib_lim1:fib_lim2,fib_lim3:fib_lim4);
                A_post=I7(fib_lim1:fib_lim2,fib_lim3:fib_lim4);
 %               fib_corr = A_fib*double(double(A_post));
                
                [boundaries_fib2,L_fib2]=bwboundaries(A_fib,4,'noholes');
                s_fib2=regionprops(L_fib2,'Area','Orientation','MajorAxisLength','MinorAxisLength');
                fib2_a=cat(1,s_fib2.Area);
                fib2_o=cat(1,s_fib2.Orientation);
                fib2_m1=cat(1,s_fib2.MajorAxisLength);
                fib2_m2=cat(1,s_fib2.MinorAxisLength);
                
%                if fib2_a < (80*80*0.50)
                    fib2_area(l,j)=fib2_area(l,j)+sum(fib2_a);
%                end

                 if isnan(mean(fib2_o));
                 else
                     fib2_Orientation(l,j)=-mean(fib2_o);
                 end
%                 
%                 if isnan(mean(fib2_m1));
%                 else    
%                     fib2_MajorAxisLength(l,j)=fib2_MajorAxisLength(l,j)+mean(fib2_m1);
%                 end
%                 
%                 if isnan(mean(fib2_m2));                
%                 else
%                     fib2_MinorAxisLength(l,j)=fib2_MinorAxisLength(l,j)+mean(fib2_m2);
%                 end
                
            end
            
        end 
    end

sum_fib=0;
sum_ext_fib=0;
unocc_post_fib=0;

for i=1:row;
    for j=1:col;
        if and(cell_posts_x(i,j)==0,cell_posts_y(i,j)==0);
            unocc_post_fib=unocc_post_fib+fib2_area(i,j);
        end
    end
end

avg_unocc_fib=unocc_post_fib/border_num;

for i=1:row;
    for j=1:col;
        
        if fib2_area(i,j)-avg_unocc_fib > 0;
            fib2_area(i,j)=fib2_area(i,j)-avg_unocc_fib;
        else
            fib2_area(i,j)=0;
        end
        
         if and(cell_posts_mag(i,j) ~= 0,fib2_area(i,j) ~=0);                
             sum_fib=sum_fib+fib2_area(i,j)/box_ratio(i,j);          
         end
        
        if exterior_posts_mag(i,j) ~=0;
            sum_ext_fib=sum_ext_fib+fib2_area(i,j)/box_ratio(i,j);
        end
        
        
    end
end



average_fib=sum_fib/cell_num;
sum_int_fib=sum_fib-sum_ext_fib;


k1=1;
k2=1;
k3=1;

force_exterior_dir_list=[];
fn_exterior_dir_list=[];
force_interior_dir_list=[];
fn_interior_dir_list=[];
force_unocc_dir_list=[];
fn_unocc_dir_list=[];

force_exterior_list=[];
fn_exterior_list=[];
force_interior_list=[];
fn_interior_list=[];
force_unocc_list=[];
fn_unocc_list=[];

for i=1:row;
    for j=1:col;
%         fib2_x(i,j)=fib2_MajorAxisLength(i,j)*-sin(fib2_Orientation(i,j)/180*pi);
%         fib2_y(i,j)=fib2_MajorAxisLength(i,j)*cos(fib2_Orientation(i,j)/180*pi);
%         fn_corr(i,j)=fib2_area(i,j)*force_mag(i,j)/(average_fib*avg_post_force);
        
        if exterior_posts_mag(i,j) ~= 0;
            fn_exterior_list(k1)=fib2_area(i,j)*convert^2;
            force_exterior_list(k1)=force_mag(i,j);
            if and(fib2_area(i,j)*convert^2 > 2.5,or(force_mag(i,j) > 5, border_posts_mag(i,j) > 5));
                fn_exterior_dir_list(k1)=fib2_Orientation(i,j);
                if and(force_dir(i,j) <= 90, force_dir(i,j) >= -90);
                    force_exterior_dir_list(k1)=force_dir(i,j);
                elseif force_dir(i,j) > 90;
                    force_exterior_dir_list(k1)=force_dir(i,j)-180;
                elseif force_dir(i,j) < -90;
                    force_exterior_dir_list(k1)=force_dir(i,j)+180;
                end
            end
            k1=k1+1;
        elseif cell_posts_mag(i,j) ~=0;
            fn_interior_list(k2)=fib2_area(i,j)*convert^2;
            force_interior_list(k2)=force_mag(i,j);
            if and(fib2_area(i,j)*convert^2 > 2.5,or(force_mag(i,j) > 5, border_posts_mag(i,j) > 5));
                fn_interior_dir_list(k2)=fib2_Orientation(i,j);
                if and(force_dir(i,j) <= 90, force_dir(i,j) >= -90);
                    force_interior_dir_list(k2)=force_dir(i,j);
                elseif force_dir(i,j) > 90;
                    force_interior_dir_list(k2)=force_dir(i,j)-180;
                elseif force_dir(i,j) < -90;
                    force_interior_dir_list(k2)=force_dir(i,j)+180;
                end
            end
            k2=k2+1;
        else
            fn_unocc_list(k3)=fib2_area(i,j)*convert^2;
            force_unocc_list(k3)=force_unocc_mag(i,j);
            if and(fib2_area(i,j)*convert^2 > 2.5,or(force_mag(i,j) > 5, border_posts_mag(i,j) > 5));
                fn_unocc_dir_list(k3)=fib2_Orientation(i,j);
                if and(border_posts_dir(i,j) <= 90, border_posts_dir(i,j) >= -90);
                    force_unocc_dir_list(k3)=border_posts_dir(i,j);
                elseif border_posts_dir(i,j) > 90;
                    force_unocc_dir_list(k3)=border_posts_dir(i,j)-180;
                elseif border_posts_dir(i,j) < -90;
                    force_unocc_dir_list(k3)=border_posts_dir(i,j)+180;
                end        
            end
            k3=k3+1;
        end
    end
end



