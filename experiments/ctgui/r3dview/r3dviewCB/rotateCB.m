function rotateCB
% callback for trackViewer

imgFigureH=findobj('Type','figure','Name','View Panel');
if isempty(imgFigureH)
   return;
end;
img=GetUserData(imgFigureH,'img');

%rotate image
img=permute(img,[2 3 1 4 5]);
SetUserData(imgFigureH,img,1);

%readjust stackslider and edit values
stackslideH=findobj(gcbf,'Tag','stackslider');
stackeditH=findobj(gcbf,'Tag','stackedit');

set(stackslideH,'Min',1);
set(stackslideH,'Max',size(img,3));
set(stackslideH,'SliderStep',[1/size(img,3) 10/size(img,3)]);
set(stackslideH,'Value',1);
set(stackeditH,'String','1');

%get timepoint
timesliderH=findobj(gcbf,'Style','slider','Tag','timeslider');
tp=get(timesliderH,'Value');

%display z-slice
imgH=findobj(imgFigureH,'Type','image');
set(imgH,'CData',img(:,:,1,1,tp));
figure(imgFigureH);
