function manuelpostpro(hObject)


picturecountH =findall(0,'Style','text','Tag','picturecount');
pictureslideH =findall(0,'Style','slider','Tag','pictureslide');

if ~isempty(picturecountH)
    delete(picturecountH);
end
if ~isempty(pictureslideH)
    delete(pictureslideH);
end


handles = guidata(hObject)

handles.listofcells='...';



imagedirectory=handles.jobvalues.imagedirectory;
imagename=handles.jobvalues.imagename;
firstimage=handles.jobvalues.firstimage;
lastimage=handles.jobvalues.lastimage;
increment=handles.jobvalues.increment;

ma=floor((lastimage-firstimage)/increment+0.001);
handles.ma=ma;
slider_step(1) = 1/(ma-1);
slider_step(2) = 3/(ma-1);


guidata(hObject, handles);

figure,

picturecountH=uicontrol('Style','text',...
                        'Units','normalized',...
                        'Tag','picturecount',...
                        'Position',[0.02,0.93,0.05,0.06])
                     
set(picturecountH,'String',num2str(firstimage));



pictureslideH=uicontrol('Style','slider',...
                        'Units','normalized',... 
                        'Value',0.00000001,...
                        'Min',0,...
                        'Max',1,...
                        'SliderStep',slider_step,...
                        'Callback','changeframe',...
                        'Tag','pictureslide',...
                        'Position',[0.02,0.02,0.05,0.9])
                     

