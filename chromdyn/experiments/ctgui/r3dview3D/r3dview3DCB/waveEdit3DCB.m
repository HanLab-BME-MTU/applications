function timeSlide
% callback for trackViewer

figureH=findobj('Type','figure','Tag','r3dView');
if isempty(figureH)
   return;
end;
axesH=findobj('Type','axes');
if isempty(axesH)
   return;
end;
DataAspectRatioMode=get(axesH,'DataAspectRatioMode');
if(strcmp(DataAspectRatioMode,'auto'))
   set(axesH,'DataAspectRatio',[1 1 1],'DataAspectRatioMode','manual');
else
   set(axesH,'DataAspectRatioMode','auto');
end;

   