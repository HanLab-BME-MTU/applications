function thresSlide3DCB
% callback for trackViewer

%get handles
upthresslideH=findobj(gcbf,'Style','slider','Tag','upthresslider');
upthreseditH=findobj(gcbf,'Style','edit','Tag','upthresedit');
thresslideH=findobj(gcbf,'Style','slider','Tag','thresslider');
threseditH=findobj(gcbf,'Style','edit','Tag','thresedit');

tag=get(gcbo,'Tag');
switch tag
case 'thresslider'
   %lower thres
   upthres=get(upthresslideH,'Value');
   lthres=round(get(gcbo,'Value'));
   %make sure that lthers<=upthres
   upthres=max(upthres,lthres);
case 'upthresslider'
   %upper thres
   lthres=get(thresslideH,'Value');
   upthres=round(get(gcbo,'Value'));
   %make sure that upthres>=lthres
   lthres=min(lthres,upthres);
end;

set(thresslideH,'Value',lthres);
set(threseditH,'String',num2str(lthres));
set(upthresslideH,'Value',upthres);
set(upthreseditH,'String',num2str(upthres));

% refresh view panels
drawr3dview3d;

