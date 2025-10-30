% Set up initial screen for input values:

figure;
maximize(gcf);

set(gcf,'Units','normalized');

vinc_flag1 = 0;
fn_flag1 = 0;
back_flag=0;

prefix_names2={};

title1=uicontrol('Style','Text','String','Waiting for User Inputs','HorizontalAlignment','center','Units','normalized','FontName','Times New Roman');
title1.Position= [500/1152 750/800 200/1152 20/800];

main_title=uicontrol('Style','Text','String','This software needs 3 files with a common prefix identifier and the following suffixes: "Top.tif", "Bottom.tif", and "Actin.tif".   Optional Fibronectin analysis and Vinculin analysis need "FN.tif", and "Vinc.tif".   Files must be in the current directoty and must be 12-bit TIFFs.','Units','normalized','FontName','Times New Roman','FontSize',10);
main_title.Position=[500/1152 400/800 200/1152 160/800];

vinc_title=uicontrol('Style','Text','String','Analyze Vinculin Image?','Units','normalized','FontName','Times New Roman');
vinc_title.Position=[20/1152 740/800 170/1152 40/800];

vinc1_name=uicontrol('Style','radiobutton','String','Yes','Units','normalized','Max',1,'Min',0,'Callback','vinc1','FontName','Times New Roman');
vinc1_name.Position=[30/1152 740/800 50/1152 20/800];

vinc2_name=uicontrol('Style','radiobutton','String','No','Units','normalized','Max',1,'Min',0,'Callback','vinc2','FontName','Times New Roman','Value',1);
vinc2_name.Position=[110/1152 740/800 50/1152 20/800];

fn_title=uicontrol('Style','Text','String','Analyze Fibronectin Image?','Units','normalized','FontName','Times New Roman');
fn_title.Position=[20/1152 690/800 170/1152 40/800];

fn1_name=uicontrol('Style','radiobutton','String','Yes','Max',1,'Min',0,'Callback','fn1','Units','normalized','FontName','Times New Roman');
fn1_name.Position=[30/1152 690/800 50/1152 20/800];

fn2_name=uicontrol('Style','radiobutton','String','No','Max',1,'Min',0,'Callback','fn2','Units','normalized','FontName','Times New Roman','Value',1);
fn2_name.Position=[110/1152 690/800 50/1152 20/800];

if repeat_flag == 1;
    set(vinc1_name,'Enable','off');
    set(vinc2_name,'Enable','off');
    set(fn1_name,'Enable','off');
    set(fn2_name,'Enable','off');
end

if repeat_flag==0;
    area_lim=0.05;
    ecc_lim=0.02;
    cent_lim=0.5;
    area_max=1000;
    area_min=400;
    ecc_max=0.6;
end

lim_title=uicontrol('Style','Text','String','Edge Detection Limits','Units','normalized','FontName','Times New Roman');
lim_title.Position=[20/1152 280/800 170/1152 400/800];

lim_title1=uicontrol('Style','Text','String','% Change in Area:','Units','normalized','FontName','Times New Roman');
lim_title1.Position=[20/1152 440/800 170/1152 20/800];
lim1=uicontrol('Style','Edit','BackgroundColor','w','String',area_lim,'Callback','arealim','Units','normalized','FontName','Times New Roman');
lim1.Position=[50/1152 420/800 110/1152 20/800];
lim_title2=uicontrol('Style','Text','String','% Change in Eccentricity:','Units','normalized','FontName','Times New Roman');
lim_title2.Position=[20/1152 380/800 170/1152 20/800];
lim2=uicontrol('Style','Edit','BackgroundColor','w','String',ecc_lim,'Callback','ecclim','Units','normalized','FontName','Times New Roman');
lim2.Position=[50/1152 360/800 110/1152 20/800];
lim_title3=uicontrol('Style','Text','String','Pixel Change in Centroid:','Units','normalized','FontName','Times New Roman');
lim_title3.Position=[20/1152 320/800 170/1152 20/800];
lim3=uicontrol('Style','Edit','BackgroundColor','w','String',cent_lim,'Callback','centlim','Units','normalized','FontName','Times New Roman');
lim3.Position=[50/1152 300/800 110/1152 20/800];
lim_title4=uicontrol('Style','Text','String','Maximum Pixel Area Per Post:','Units','normalized','FontName','Times New Roman');
lim_title4.Position=[20/1152 620/800 170/1152 20/800];
lim4=uicontrol('Style','Edit','BackgroundColor','w','String',area_max,'Callback','arealim1','Units','normalized','FontName','Times New Roman');
lim4.Position=[50/1152 600/800 110/1152 20/800];
lim_title5=uicontrol('Style','Text','String','Minimum Pixel Area Per Post:','Units','normalized','FontName','Times New Roman');
lim_title5.Position=[20/1152 560/800 170/1152 20/800];
lim5=uicontrol('Style','Edit','BackgroundColor','w','String',area_min,'Callback','arealim2','Units','normalized','FontName','Times New Roman');
lim5.Position=[50/1152 540/800 110/1152 20/800];
lim_title6=uicontrol('Style','Text','String','Maximum Eccentricity Per Post:','Units','normalized','FontName','Times New Roman');
lim_title6.Position=[20/1152 500/800 170/1152 20/800];
lim6=uicontrol('Style','Edit','BackgroundColor','w','String',ecc_max,'Callback','ecclim1','Units','normalized','FontName','Times New Roman');
lim6.Position=[50/1152 480/800 110/1152 20/800];

prefix_title=uicontrol('Style','Text','String','Enter Prefix Identifier:','Units','normalized','FontName','Times New Roman');
prefix_title.Position=[20/1152 210/800 170/1152 60/800];

prefix_name=uicontrol('Style','Edit','BackgroundColor','w','Callback','prefix','Units','normalized','FontName','Times New Roman');
prefix_name.Position=[50/1152 230/800 110/1152 20/800];
if repeat_flag ==1
    set(prefix_name,'string',file);
end

back_title=uicontrol('Style','Text','String','Background Subtract Images?','Units','normalized','FontName','Times New Roman');
back_title.Position=[20/1152 160/800 170/1152 40/800];
back1_name=uicontrol('Style','radiobutton','String','Yes','Units','normalized','Max',1,'Min',0,'Callback','back1','FontName','Times New Roman');
back1_name.Position=[30/1152 160/800 50/1152 20/800];
back2_name=uicontrol('Style','radiobutton','String','No','Units','normalized','Max',1,'Min',0,'Callback','back2','FontName','Times New Roman','Value',1);
back2_name.Position=[110/1152 160/800 50/1152 20/800];

if repeat_flag == 0
    postd=3;
    circ_rad=postd/2;
end

post_dia_name=uicontrol('Style','Text','String','Post Diameter (in microns):','Units','normalized','FontName','Times New Roman');
post_dia_name.Position=[20/1152 110/800 170/1152 40/800];
post_dia=uicontrol('Style','Edit','BackgroundColor','w','String',postd,'Callback','postdia','Units','normalized','FontName','Times New Roman');
post_dia.Position=[50/1152 110/800 110/1152 20/800];

if repeat_flag == 0
    pin_spacing=9;
end
accept1=uicontrol('Style','pushbutton','String','OK','Callback','accept','Units','normalized','FontName','Times New Roman');
accept1.Position=[80/1152 30/800 50/1152 20/800];

post_spacing_name=uicontrol('Style','Text','String','Post to Post Spacing (in microns):','Units','normalized','FontName','Times New Roman');
post_spacing_name.Position=[20/1152 60/800 170/1152 40/800];
post_spacing=uicontrol('Style','Edit','BackgroundColor','w','String',pin_spacing,'Callback','postspace','Units','normalized','FontName','Times New Roman');
post_spacing.Position=[50/1152 60/800 110/1152 20/800];

if repeat_flag ==1
    set(post_dia,'Enable','off');
    set(post_spacing,'Enable','off');
end

dir_struct=dir;

j=1;

for i=1:length(dir_struct)
    name_len=length(dir_struct(i).name);
    if name_len > 4
        if dir_struct(i).name(name_len-3:name_len) == '.tif'
            dir_names(j)={dir_struct(i).name(1:name_len-4)};
            j=j+1;
        end
    end
end

j=1;

%for i=1:length(dir_names);
%    name_len2=length(dir_names{i});
%    if dir_names(i).(name_len2-4:name_len2)==string('Actin');
%        prefix_names(j)={dir_names{i}(1:name_len2-5)};
%        j=j+1;
%    end
%end

uiwait;

