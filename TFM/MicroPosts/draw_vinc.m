% Draw vinc: determine vinc assembly in cell


vinc=imread(vinc_file);

vinc1=imadjust(vinc);

vinc2a=imrotate(vinc1,rot_angle);
vinc2b=imcrop(vinc2a,[p1 offset]);

gray_factor_vinc=1;

% Automatically threshold Actin image

vinc3=im2bw(vinc2b,gray_factor_vinc*graythresh(vinc2b));

vinc4=vinc3;

vinc4b=vinc4;

subplot(2,1,1); imshow(vinc2b); subplot(2,1,2); imshow(vinc4);
maximize(gcf);

% Prompt for adjustments to threshold value

vinc_ok=questdlg('Is vinculin shown accurately?');

if vinc_ok(1:2) == 'No';
    %gf=inputdlg('Input new scaling factor (default=1.0):');
    %gray_factor=str2num(gf{1});
    hui_v=uicontrol('style','slider','Tag','ThreshSlide','Value',1,'min',0,'max',1/graythresh(vinc2b),'Position',[500/1152 30/800 200/1152 20/800],'Units','normalized','Callback','m0v2');
    hui_vc=uicontrol('Style', 'text','Tag','ThreshLabel','String', 'Threshold (0-1):','Position',[500/1152 50/800 200/1152 20/800],'Units','normalized');
    hui_va=uicontrol('style','edit','Tag','ThreshValue','String',sprintf('%3.2f',gray_factor_vinc*graythresh(vinc2b)),'Position',[560/1152 10/800 80/1152 20/800],'Units','normalized','callback','thresh_input_v');
    hui_vb=uicontrol('style','pushbutton','String','Reset','Position',[400/1152 30/800 40/1152 20/800],'Units','normalized','Enable','on','Callback','reset_gray_v');
    
%     if gray_factor*graythresh(actin2) < 1.0;
%         actin3=im2bw(actin2,gray_factor*graythresh(actin2));
%     else
%         actin3=im2bw(actin2,1);
%     end
%     actin4=imfill(actin3,'holes');
%     subplot(2,1,1); imshow(actin2); subplot(2,1,2); imshow(actin4);
%     maximize(gcf);

    hui_v2=uicontrol('style','pushbutton','String','OK','Position',[800/1152 30/800 40/1152 20/800],'Units','normalized','Enable','on','Callback','accept2v');
    
    uiwait;
end

subplot(1,1,1);