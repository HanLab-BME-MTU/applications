% Draw Actin: determine outline of cell

se=strel('disk',3);

actin=imread(actin_file);

% actin1=imadjust(actin);
actin1=imadjust(actin,[double(min(min(actin)))/2^16 double(max(max(actin)))/2^16],[0 1]);

actin2=actin1;

gray_factor=1;

% Automatically threshold Actin image

actin3=im2bw(actin2,gray_factor*graythresh(actin2));
actin4=imfill(actin3,'holes');

subplot(2,1,1); imshow(actin2); subplot(2,1,2); imshow(actin4);
maximize(gcf);

% Prompt for adjustments to threshold value

actin_ok=questdlg('Is cell outline shown accurately?');

if actin_ok(1:2) == 'No';
    %gf=inputdlg('Input new scaling factor (default=1.0):');
    %gray_factor=str2num(gf{1});
    hui=uicontrol('style','slider','Tag','ThreshSlide','Value',1,'min',0,'max',1/graythresh(actin2),'Position',[500/1152 30/800 200/1152 20/800],'Units','normalized','Callback','m02');
    hui_c=uicontrol('Style', 'text','Tag','ThreshLabel','String', strcat('Threshold (0-',num2str(1/graythresh(actin2)),'):'),'Position',[500/1152 50/800 200/1152 20/800],'Units','normalized');
    hui_a=uicontrol('style','edit','Tag','ThreshValue','String',sprintf('%3.2f',gray_factor*graythresh(actin2)),'Position',[560/1152 10/800 80/1152 20/800],'Units','normalized','callback','thresh_input');
    hui_b=uicontrol('style','pushbutton','String','Reset','Position',[400/1152 30/800 40/1152 20/800],'Units','normalized','Enable','on','Callback','reset_gray');
    
%     if gray_factor*graythresh(actin2) < 1.0;
%         actin3=im2bw(actin2,gray_factor*graythresh(actin2));
%     else
%         actin3=im2bw(actin2,1);
%     end
%     actin4=imfill(actin3,'holes');
%     subplot(2,1,1); imshow(actin2); subplot(2,1,2); imshow(actin4);
%     maximize(gcf);

    hui2=uicontrol('style','pushbutton','String','OK','Position',[800/1152 30/800 40/1152 20/800],'Units','normalized','Enable','on','Callback','accept2');
    
    uiwait;
end

