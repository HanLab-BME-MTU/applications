function leveldeterminer(hObject)

handles = guidata(hObject);

projNum = get(handles.GUI_st_job_lb,'Value');
path=handles.jobs(projNum).imagedirectory;
main=handles.jobs(projNum).imagename;
first=handles.jobs(projNum).firstimage;
last=handles.jobs(projNum).lastimage;



cd(path);

ext=main(end-3:end);
main=main(1:end-7);

s=3; %s=length(num2str(no));
strg=sprintf('%%.%dd',s);
% Create numerical index

indxStr=sprintf(strg,first);   
name=[main indxStr ext]; 
firstimg=imread(name);

indxStr=sprintf(strg,last);
name=[main indxStr ext]; 
    
lastimg=imread(name);


figure

intensiback = impixel(firstimg,[]), title('Click on the background (a lot of times). Make sure your clicks are spread out evenly. Then press enter') ;
handles.jobs(projNum).fi_background= sum(intensiback(:,1)/size(intensiback,1));
clear intensiback
close


figure 

intensinuc = impixel(firstimg,[]), title('Click on the nucleoi (a lot of times). Make sure your clicks on a lot of different ones. Then press enter');
handles.jobs(projNum).fi_nucleus= sum(intensinuc(:,1)/size(intensinuc,1));
clear intensinuc
close


figure

intensiback = impixel(lastimg,[]), title('Click on the background (a lot of times). Make sure your clicks are spread out evenly. Then press enter') ;
handles.jobs(projNum).la_background= sum(intensiback(:,1)/size(intensiback,1));
clear intensiback
close


figure,title('Click on the nucleoi (a lot of times).')

intensinuc = impixel(lastimg,[]),title('Click on the nucleoi (a lot of times).') ;
handles.jobs(projNum).la_nucleus= sum(intensinuc(:,1)/size(intensinuc,1));
clear intensinuc
close








                  
% Update handles structure
guidata(hObject, handles);
