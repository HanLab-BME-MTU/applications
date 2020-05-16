% Display final results, set up buttons for different tabs

perim=0;
for i=1:length(s_actin);
perim=perim+length(boundaries_actin{i})*convert;
end

toggle_flag = 0;

end_h=uicontrol('Style','Text','String',{'Total Force:','Total Area:','Force Per Post:','Std Dev, Unoccupied Posts:', 'Std Dev, Bottom Centers:', 'Force Imbalance:', 'Force Imbalance/Post:'; sum_total_force, sum(cat(1,s_actin.Area))*convert^2, avg_post_force, stdev, std_bot, sum_net_force_mag, avg_post_err; 'nN','um^2','nN','um','um','nN','nN'; ' ',' ',' ',' ',' ',' ',' '},'Position',[20/1152 400/800 170/1152 370/800],'Units','normalized','FontSize',8,'FontName','Times New Roman');
end_h2=uicontrol('Style','Text','String',{'Mean Post Area, Top:','Mean Post Area, Bottom:','Mean Eccentricity, Top:','Mean Eccentricity, Bottom:','# of Posts Excluded, Top:','# of Posts Excluded, Bottom:','Elapsed Time:'; mean(mean(final_areaA))*convert^2,mean(mean(final_areaB))*convert^2,mean(mean(final_eccA)),mean(mean(final_eccB)),size(bad_posts_top,1),size(bad_posts_bot,1),time; 'um^2','um^2',' ',' ',' ',' ' 's';' ',' ',' ',' ',' ',' ',' '},'Position',[20/1152 10/800 170/1152 370/800],'Units','normalized','FontSize',8,'FontName','Times New Roman');
end_h3=uicontrol('Style','Text','String',{'Cell Perimeter','Fibronectin Area','Exterior Post Force','Interior Post Force','Exterior Fibronectin'; perim,sum_fib*convert^2,sum_ext_force,sum_int_force,sum_ext_fib*convert^2; 'um','um^2','nN','nN','um^2';' ',' ',' ',' ',' '},'Position',[20/1152 10/800 170/1152 370/800],'Units','normalized','FontSize',8,'FontName','Times New Roman','Visible','Off');


top_h=uicontrol('Style','Pushbutton','String','Top Image','Position',[1030/1152 700/800 120/1152 40/800],'Callback','show_posts_top_w_NS','Units','normalized','FontName','Times New Roman');

bot_h=uicontrol('Style','Pushbutton','String','Bottom Image','Position',[1030/1152 650/800 120/1152 40/800],'Callback','show_posts_bot_w_NS','Units','normalized','FontName','Times New Roman');

actin_h=uicontrol('Style','Pushbutton','String','Cell/Post/Vector Composite','Position',[1030/1152 600/800 120/1152 40/800],'Callback','show_actin_vectors','Units','normalized','FontName','Times New Roman','FontSize',7);

fn_h=uicontrol('Style','Pushbutton','String','FN Image','Position',[1030/1152 550/800 120/1152 40/800],'Enable','off','Units','normalized','FontName','Times New Roman','Callback','show_fn');

if fn_flag1==1;
    set(fn_h,'Enable','on');
end

vinc_h=uicontrol('Style','Pushbutton','String','Vinc Image','Position',[1030/1152 500/800 120/1152 40/800],'Enable','off','Units','normalized','FontName','Times New Roman','Callback','show_vinc');

if vinc_flag1==1;
     set(vinc_h,'Enable','on');
end

def_tab_h=uicontrol('Style','Pushbutton','String','Deflection Table','Position',[1030/1152 450/800 120/1152 40/800],'Callback','show_def','Units','normalized','FontName','Times New Roman');

cntr_bot_tab_h=uicontrol('Style','Pushbutton','String','Top Centroids Table','Position',[1030/1152 400/800 120/1152 40/800],'Callback','show_cntr_top','Units','normalized','FontName','Times New Roman');

cntr_top_tab_h=uicontrol('Style','Pushbutton','String','Bottom Centroids Table','Position',[1030/1152 350/800 120/1152 40/800],'Callback','show_cntr_bot','Units','normalized','FontName','Times New Roman','FontSize',7);

rerun_tab_h=uicontrol('Style','Pushbutton','String','Re-run Data','Position',[1030/1152 300/800 120/1152 40/800],'Callback','analyze_posts_repeat','Units','normalized','FontName','Times New Roman');

actin2_tab_h=uicontrol('Style','Pushbutton','String','Actin Image','Position',[1030/1152 250/800 120/1152 40/800],'Units','normalized','FontName','Times New Roman','Callback','show_actin');

new_stats_tab_h=uicontrol('Style','Pushbutton','String','Toggle Stats','Position',[1030/1152 200/800 120/1152 40/800],'Units','normalized','FontName','Times New Roman','Callback','new_stats');


fn_plot_tab_h=uicontrol('Style','Pushbutton','String','FN Area Plot','Position',[1030/1152 150/800 120/1152 40/800],'Enable','off','Units','normalized','FontName','Times New Roman','Callback','fn_plot');

if fn_flag1==1;
    set(fn_plot_tab_h,'Enable','on');
end

vinc_plot_tab_h=uicontrol('Style','Pushbutton','String','Vinc Area Plot','Position',[1030/1152 100/800 120/1152 40/800],'Enable','off','Units','normalized','FontName','Times New Roman','Callback','vinc_plot');

if vinc_flag1==1;
    set(vinc_plot_tab_h,'Enable','on');
end






% set(main_title,'Visible','off');
set(vinc_title,'Visible','off');
set(vinc1_name,'Visible','off');
set(vinc2_name,'Visible','off');
set(fn_title,'Visible','off');
set(fn1_name,'Visible','off');
set(fn2_name,'Visible','off');
set(prefix_title,'Visible','off');
set(prefix_name,'Visible','off');
set(accept1,'Visible','off');
set(lim_title,'Visible','off');
set(lim_title1,'Visible','off');
set(lim_title2,'Visible','off');
set(lim_title3,'Visible','off');
set(lim_title4,'Visible','off');
set(lim_title5,'Visible','off');
set(lim_title6,'Visible','off');
set(lim1,'Visible','off');
set(lim2,'Visible','off');
set(lim3,'Visible','off');
set(lim4,'Visible','off');
set(lim5,'Visible','off');
set(lim6,'Visible','off');
