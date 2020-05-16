% set(title1,'String','Composite Actin Image with Vectors');

post_top_flag=0;
post_bot_flag=0;
actin_flag=1;
fn_flag=0;
vinc_flag=0;
fn_plot_flag=0;
vinc_plot_flag=0;


if def_flag == 1;
    for i=1:row;
        for j=1:col;
            set(deflectionsx(i,j),'Visible','off');
            set(deflectionsy(i,j),'Visible','off');
        end
    end
end

if top_flag == 1;
    for i=1:row;
        for j=1:col;
            set(centers_topx_h(i,j),'Visible','off');
            set(centers_topy_h(i,j),'Visible','off');
        end
    end
end

if bot_flag == 1;
    for i=1:row;
        for j=1:col;
            set(centers_botx_h(i,j),'Visible','off');
            set(centers_boty_h(i,j),'Visible','off');
        end
    end
end

hold off;

subplot(1,1,1);
f4_top=imadjust(f_top,[double(min(min(f_top)))/2^16 double(max(max(f_top)))/2^16],[0 1]);
f4a_top=imrotate(f4_top,rot_angle);
f4b_top=imcrop(f4a_top,[p1 offset]);
comp4(:,:,1)=f4b_top(:,:);

f4_bot=imadjust(f_bot,[double(min(min(f_bot)))/2^16 double(max(max(f_bot)))/2^16],[0 1]);
f4a_bot=imrotate(f4_bot,rot_angle);
f4b_bot=imcrop(f4a_bot,[p1 offset]);
comp4(:,:,3)=f4b_bot(:,:);

post_actin_h1=imshow(comp4);

hold on;

for i=1:length(s_actin);
    b=boundaries_actin{i};
    if length(b)>(pi*3/convert);
        post_actin_h(i)=plot(b(:,2),b(:,1),'y');
        b_actin(k)=i;
        k=k+1;
    end
end

mag_factor = 10;

     h=quiver(actual_y_bot,actual_x_bot,mag_factor*exterior_posts_y/convert,mag_factor*exterior_posts_x/convert,0,'k');
        if size(h,1) == 2;
            set(h(1),'Linewidth',2);
            set(h(2),'Linewidth',2);
        else
            set(h(1),'Linewidth',2);
        end
     h2=quiver(actual_y_bot,actual_x_bot,mag_factor*border_posts_y/convert,mag_factor*border_posts_x/convert,0,'r');
        if size(h2,1) == 2;
            set(h2(1),'Linewidth',2);
            set(h2(2),'Linewidth',2);
        else
            set(h2(1),'Linewidth',2);
        end
     h3=quiver(actual_y_bot,actual_x_bot,mag_factor*interior_posts_y/convert,mag_factor*interior_posts_x/convert,0,'w');
        if size(h3,1) == 2;
            set(h3(1),'Linewidth',2);
            set(h3(2),'Linewidth',2);
        else
            set(h3(1),'Linewidth',2);
        end
