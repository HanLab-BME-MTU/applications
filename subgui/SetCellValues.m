function SetCellValues(hObject,jobdeterminer)

change=0;
handles = guidata(hObject);

projNum = get(handles.GUI_st_job_lb,'Value');
imageDirectory=handles.jobs(projNum).imagedirectory;
imageName=handles.jobs(projNum).imagename;
firstImageNum=handles.jobs(projNum).firstimage;
last=handles.jobs(projNum).lastimage;



cd(imageDirectory);

ext=imageName(end-3:end);
main=imageName(1:end-7);

s=3; %s=length(num2str(no));
strg=sprintf('%%.%dd',s);
% Create numerical index

indxStr=sprintf(strg,firstImageNum);   
name=[main indxStr ext]; 
firstImg=imread(name);


[img_h,img_w]=size(firstImg);

indxStr=sprintf(strg,last);
name=[main indxStr ext]; 
    
lastImg=imread(name);



%get the intensity values of these pictures

%first picture
%background
if jobdeterminer==1
    figure, imshow(firstImg),title('Make a polygon around the smallest nucloi, by clicking on its rimm, going around clock- or anticlockwise,then press ENTER') ; 
    BW =roipoly;
    close
    if ~length(find(BW))==0
        handles.jobs(projNum).minsize= length(find(BW));
        change=1;
    end
end
if jobdeterminer==2
    figure, imshow(firstImg), title('Make a polygon around the biggest nucloi, by clicking on its rimm, going around clock- or anticlockwise,then press ENTER') ;
     BW =roipoly;
    close
    if ~length(find(BW))==0
        handles.jobs(projNum).maxsize= length(find(BW));
        change=1;
    end
   
end
  
if jobdeterminer==3
    figure, imshow(firstImg), title('click on the centers of the two cells closest to each other, then press ENTER') ;
     [x,y] =getpts;
     close
     
     if length(x)>1
         handles.jobs(projNum).minsdist= round(sqrt((x(1)-x(2))^2+(y(1)-y(2))^2));
         change=1;
    end
    
end

if change==1
% Update handles structure
guidata(hObject, handles);

%%%%%%%%save altered values to disk%%%%%%%%%%%%
cd(handles.jobs(projNum).savedirectory)
jobvalues=handles.jobs(projNum);
save ('jobvalues','jobvalues')
clear jobvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fillFields(handles,handles.jobs(projNum));
end