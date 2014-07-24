function nothingatall
% Comment


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



slidervalue=get(pictureslideH,'Value');
slidervalue=round(slidervalue*(ma-1)+1);

whichpic=(slidervalue-1)*increment + firstimage;

set(piccount,'String',num2str(whichpic))


bodyfilename=imagename(1:end-7);
ext=imagename(end-3:end);

% Format
s=3; %s=length(num2str(no));
strg=sprintf('%%.%dd',s);

% Create numerical index
indxStr=sprintf(strg,whichpic);


cd(imagedirectory)
name=[bodyfilename indxStr ext]; 

%get the current picture
picnew=imreadnd2(name,0,handles.jobs(projNum).intensityMax);

hold on;
imshow(picnew,[]), title(num2str(whichpic))

  
         
frup=find(handles.MPM(:,2*slidervalue-1));

frallap=zeros(size(handles.MPM,1),3);
frallap(:,3)=[1:1:size(handles.MPM,1)]';
frallap(:,1:2)=handles.MPM(:,2*slidervalue-1:2*slidervalue);


plot(frallap(frup,1),frallap(frup,2),'.');
txt= text(frallap(frup,1),frallap(frup,2),num2str(frallap(frup,3)),'Color','r');
set(txt,'ButtonDownFcn','manrelink'); 

hold off;




