function r3dview3D(img)
% R3DVIEW3D displays *.r3d image series
%
% SYNOPSIS r3dview3D(img)
%
% INPUT img : image=[x,y,z,wvl,time], if the image is omitted,
%                  a file selection dialog is opened
%

% 12/3/01 dT

% if no image, get image from file
if(nargin==0)
   [img, fname]=r3dread;
   if isempty(img)
      return;
   end;
else
   fname='img';
end;

% XY border pixels should be ignored (invalid camera pixels)
img=img(2:size(img,1)-1,2:size(img,2)-1,:,:,:);

% intialize params
%iptsetpref('ImshowBorder', 'tight');
imSize=size(img);

colordef black;

%open control panel
contFigH=r3dview3Dpanel;

%set first image
scrSize=get(0,'ScreenSize');
imgFigH=figure('Name',fname,'Position',...
   [(scrSize(3)-scrSize(3)/3-10) (scrSize(4)-scrSize(4)/3-50) scrSize(3)/3 scrSize(4)/3],'NumberTitle','off');

%set props
set(imgFigH,'Tag','r3dview3D');
set(imgFigH,'UserData',img);
thresslideH=findobj(contFigH,'Tag','thresslider');
threseditH=findobj(contFigH,'Tag','thresedit');
upthresslideH=findobj(contFigH,'Tag','upthresslider');
upthreseditH=findobj(contFigH,'Tag','upthresedit');
timeslideH=findobj(contFigH,'Tag','timeslider');
timeeditH=findobj(contFigH,'Tag','timeedit');
%set threshold and edit values
tImg=squeeze(img(:,:,:,1,1));
minIn=min(tImg(:));
maxIn=max(tImg(:));
difIn=maxIn-minIn+1;

%update histo
hi=hist(tImg(:),difIn);
%set thres to max
thres=maxIn;

set(thresslideH,'Min',minIn);
set(thresslideH,'Max',maxIn);
set(thresslideH,'SliderStep',[1/difIn 10/difIn]);
set(thresslideH,'Value',thres);
set(threseditH,'String',num2str(thres));
set(upthresslideH,'Min',minIn);
set(upthresslideH,'Max',maxIn);
set(upthresslideH,'SliderStep',[1/difIn 10/difIn]);
set(upthresslideH,'Value',thres);
set(upthreseditH,'String',num2str(thres));

%set timeslider and edit values
set(timeslideH,'Min',0.9999999);
set(timeslideH,'Max',size(img,5));
set(timeslideH,'SliderStep',[1/size(img,5) 10/size(img,5)]);
set(timeslideH,'Value',1);
set(timeeditH,'String','1');

%show img
figure('Tag','r3dhisto','Position',[10 (scrSize(4)-scrSize(4)/3-50) scrSize(3)/3 scrSize(4)/3]);
plot([1:length(hi)]+minIn,hi);
figure(imgFigH);
xlabel('X');
ylabel('Y');
zlabel('Z');
plot3data(tImg,thres);
axesH=findobj(gcf,'Type','axes');
set(axesH,'DataAspectRatio',[1 1 1],'DataAspectRatioMode','manual');
