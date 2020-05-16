% Draw FN: determine FN assembly in cell


fib=imread(fn_file);

fib1=imadjust(fib);

fib2a=imrotate(fib1,rot_angle);
fib2b=imcrop(fib2a,[p1 offset]);

gray_factor_fib=1;

% Automatically threshold Actin image
% 
%  l=0;
%  for M = 1:80:(hw(1)*80+1-80);
%     l = l + 1;
%     j = 0;
%         for N = 1:80:(hw(2)*80+1-80);
%             j = j+1;
%             
%             fib2_window=fib2b(M:M+79,N:N+79);
%             fib3_window=im2bw(fib2_window,gray_factor_fib*graythresh(fib2_window));
%             fib3(M:M+79,N:N+79)=fib3_window;
%         end
%  end
 
fib3=im2bw(fib2b,gray_factor_fib*graythresh(fib2b));

fib4=fib3;

fib4b=fib4;

subplot(2,1,1); imshow(fib2b); subplot(2,1,2); imshow(fib4);
maximize(gcf);

% Prompt for adjustments to threshold value

fn_ok=questdlg('Is fibronectin shown accurately?');

if fn_ok(1:2) == 'No';
    %gf=inputdlg('Input new scaling factor (default=1.0):');
    %gray_factor=str2num(gf{1});
    hui_f=uicontrol('style','slider','Tag','ThreshSlide','Value',1,'min',0,'max',1/graythresh(fib2b),'Position',[500/1152 30/800 200/1152 20/800],'Units','normalized','Callback','m0f2');
    hui_fc=uicontrol('Style', 'text','Tag','ThreshLabel','String', 'Threshold (0-1):','Position',[500/1152 50/800 200/1152 20/800],'Units','normalized');
    hui_fa=uicontrol('style','edit','Tag','ThreshValue','String',sprintf('%3.2f',gray_factor_fib*graythresh(fib2b)),'Position',[560/1152 10/800 80/1152 20/800],'Units','normalized','callback','thresh_input_f');
    hui_fb=uicontrol('style','pushbutton','String','Reset','Position',[400/1152 30/800 40/1152 20/800],'Units','normalized','Enable','on','Callback','reset_gray_f');
    
%     if gray_factor*graythresh(actin2) < 1.0;
%         actin3=im2bw(actin2,gray_factor*graythresh(actin2));
%     else
%         actin3=im2bw(actin2,1);
%     end
%     actin4=imfill(actin3,'holes');
%     subplot(2,1,1); imshow(actin2); subplot(2,1,2); imshow(actin4);
%     maximize(gcf);

    hui_f2=uicontrol('style','pushbutton','String','OK','Position',[800/1152 30/800 40/1152 20/800],'Units','normalized','Enable','on','Callback','accept2f');
    
    uiwait;
end

subplot(1,1,1);