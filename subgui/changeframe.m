function changeframe


delete(gca)

pictureslideH=findall(0,'Tag','pictureslide');
hObject=findall(0,'Tag','GUI_pp_jobbrowse_pb');
handles=guidata(hObject);

piccount =findall(0,'Style','text','Tag','picturecount');
                        
                        



imagedirectory=handles.jobvalues.imagedirectory;
imagename=handles.jobvalues.imagename;
firstimage=handles.jobvalues.firstimage;
lastimage=handles.jobvalues.lastimage;
increment=handles.jobvalues.increment;
ma=handles.ma;
ImageNamesList = handles.jobvalues.imagenameslist;



slidervalue=get(pictureslideH,'Value');
slidervalue=round(slidervalue*(ma-1)+1);

whichpic=(slidervalue-1)*increment + firstimage;

set(piccount,'String',num2str(whichpic))


cd(imagedirectory)


name = char(ImageNamesList(whichpic));


%get the current picture
picnew=imreadnd2(name,0,handles.jobvalues.intensityMax);

hold on;
imshow(picnew,[]), title(num2str(whichpic))

  
         
indRealCell=find(handles.MPM(:,2*slidervalue-1) | handles.MPM(:,2*slidervalue));

cellsWithNums=zeros(size(handles.MPM,1),3);
cellsWithNums(:,3)=[1:1:size(handles.MPM,1)]';
cellsWithNums(:,1:2)=handles.MPM(:,2*slidervalue-1:2*slidervalue);


plot(cellsWithNums(indRealCell,1),cellsWithNums(indRealCell,2),'r.');
txt= text(cellsWithNums(indRealCell,1),cellsWithNums(indRealCell,2),num2str(cellsWithNums(indRealCell,3)),'Color','r');
set(txt,'ButtonDownFcn','manrelink'); 

hold off;



