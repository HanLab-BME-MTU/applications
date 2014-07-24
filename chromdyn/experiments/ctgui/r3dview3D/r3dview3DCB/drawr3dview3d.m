function drawr3dview3d
imgFigureH=findobj('Type','figure','Tag','r3dview3D');
if isempty(imgFigureH)
   return;
end;
%get image
img=get(imgFigureH,'UserData');
% get timepoint
timesliderH=findobj(gcbf,'Style','slider','Tag','timeslider');
tp=get(timesliderH,'Value');

% threshold values
upthresslideH=findobj(gcbf,'Style','slider','Tag','upthresslider');
thresslideH=findobj(gcbf,'Style','slider','Tag','thresslider');
lthres=get(thresslideH,'Value');
upthres=get(upthresslideH,'Value');

%display z-slice
figure(imgFigureH);
cla;
plot3data(img(:,:,:,1,tp),lthres,upthres);
