function thresEdit3DCB
% callback for trackViewer

%get handles
upthresslideH=findobj(gcbf,'Style','slider','Tag','upthresslider');
upthreseditH=findobj(gcbf,'Style','edit','Tag','upthresedit');
thresslideH=findobj(gcbf,'Style','slider','Tag','thresslider');
threseditH=findobj(gcbf,'Style','edit','Tag','thresedit');

tag=get(gcbo,'Tag');
switch tag
case 'thresedit'
   %lower thres
   lthres=str2num(get(gcbo,'String'));   
   %make sure that lthers<=upthres and lthres>minV   
   maxV=get(upthresslideH,'Value');
   minV=get(thressliderH,'Min');
   %check for boundary
   lthres=min(maxV,lthres);
   lthres=max(minV,lthres);
   set(gcbo,'String',num2str(lthres));
   set(thressliderH,'Value',lthres);
case 'upthresedit'
   %upper thres
   upthres=str2num(get(gcbo,'String'));
   %make sure that upthres>=lthres and upthres<=maxV
   maxV=get(thresslideH,'Max');
   minV=get(thresslideH,'Value');
   %check for boundary
   upthres=min(maxV,upthres);
   upthres=max(minV,upthres);
   set(gcbo,'String',num2str(upthres));
   set(upthresslideH,'Value',upthres);
end;

% refresh view panels
drawr3dview3d;