% CAL 1/24/03
%
% Displays posts image with centers labeled and numbered

% set(title1,'String','Top Image with Labeled Centers');

post_top_flag=1;
post_bot_flag=0;
actin_flag=0;
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

post_top_h=imshow(II);
maximize(gcf);
hold on;
post_top_h2=plot(centr2_top(:,1),centr2_top(:,2),'b*');

%for i=1:length(top_bound);
%   tb1=top_bound{i};
%   plot(tb1(:,2),tb1(:,1),'b','Linewidth',1);
%end



hold off;



%axis([-50 1400 -50 1200]);
%title(file);

axis tight;

