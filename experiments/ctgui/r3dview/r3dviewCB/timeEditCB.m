function timeEditCB
% callback for trackViewer

imgFigureH=findobj('Type','figure','Name','View Panel');
if isempty(imgFigureH)
   return;
end;
val=str2num(get(gcbo,'String'));
timesliderH=findobj(gcbf,'Style','slider','Tag','timeslider');
maxV=get(timesliderH,'Max');
minV=get(timesliderH,'Min');
val=min(maxV,val);
val=max(minV,val);
set(gcbo,'String',num2str(val));
set(timesliderH,'Value',val);
img=GetUserData(imgFigureH,'img');

stackliderH=findobj(gcbf,'Style','slider','Tag','stackslider');
sp=get(stackliderH,'Value');

%display z-slice
imgH=findobj(imgFigureH,'Type','image');
set(imgH,'CData',img(:,:,sp,1,val));
figure(imgFigureH);
   