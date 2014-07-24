function stackEditCB
% callback for trackViewer

imgFigureH=findobj('Type','figure','Name','View Panel');
if isempty(imgFigureH)
   return;
end;
val=str2num(get(gcbo,'String'));
stacksliderH=findobj(gcbf,'Style','slider','Tag','stackslider');

timesliderH=findobj(gcbf,'Style','slider','Tag','timeslider');
tp=get(timesliderH,'Value');

maxV=get(stacksliderH,'Max');
minV=get(stacksliderH,'Min');
%check for boundary
val=min(maxV,val);
val=max(minV,val);
set(gcbo,'String',num2str(val));
set(stacksliderH,'Value',val);
img=GetUserData(imgFigureH,'img');

%display z-slice
imgH=findobj(imgFigureH,'Type','image');
set(imgH,'CData',img(:,:,val,1,tp));
figure(imgFigureH);
   