function timeEdit3DCB
% callback for trackViewer

histFigH=findobj('Type','figure','Tag','r3dhisto');
imgFigureH=findobj('Type','figure','Tag','r3dview3D');
if isempty(imgFigureH)
   return;
end;
val=str2num(get(gcbo,'String'));
timesliderH=findobj(gcbf,'Style','slider','Tag','timeslider');
maxV=round(get(timesliderH,'Max'));
minV=round(get(timesliderH,'Min'));
val=min(maxV,val);
val=max(minV,val);
set(gcbo,'String',num2str(val));
set(timesliderH,'Value',val);
img=get(imgFigureH,'UserData');

%update threshold and edit values
thresslideH=findobj(gcbf,'Style','slider','Tag','thresslider');
threseditH=findobj(gcbf,'Style','edit','Tag','thresedit');
upthresslideH=findobj(gcbf,'Style','slider','Tag','upthresslider');
upthreseditH=findobj(gcbf,'Style','edit','Tag','upthresedit');

%use current thres value
lthres=get(thresslideH,'Value');
upthres=get(upthresslideH,'Value');

tImg=squeeze(img(:,:,:,1,val));
minIn=min(tImg(:));
maxIn=max(tImg(:));
difIn=maxIn-minIn+1;
%check boundary
upthres=min(upthres,maxIn);
lthres=min(lthres,upthres);
lthres=max(lthres,minIn);
upthres=max(upthres,lthres);

%lower thres
set(thresslideH,'Min',minIn);
set(thresslideH,'Max',maxIn);
set(thresslideH,'SliderStep',[1/difIn 10/difIn]);
set(thresslideH,'Value',lthres);
set(threseditH,'String',num2str(lthres));
%upper thres
set(upthresslideH,'Min',minIn);
set(upthresslideH,'Max',maxIn);
set(upthresslideH,'SliderStep',[1/difIn 10/difIn]);
set(upthresslideH,'Value',upthres);
set(upthreseditH,'String',num2str(upthres));

%update histo
hi=hist(tImg(:),difIn);

%display z-slice
if isempty(histFigH)
   figure('Tag','r3dhisto');
else
   figure(histFigH);
end;
plot([1:length(hi)]+minIn,hi);

% refresh view panels
drawr3dview3d;