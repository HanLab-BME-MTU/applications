function SetCellValues(hObject,jobdeterminer)


%this routine enables the user to specify:
%-Minimal size of nuclei (1)
%-Maximal size of nuclei (2)
%-Minimal distance between two nuclei (3)
%interactively. Which one of the three it is going to be, depends on the
%value of jobdeterminer

change=0;
handles = guidata(hObject);

projNum = get(handles.GUI_st_job_lb,'Value');
imageDirectory=handles.jobs(projNum).imagedirectory;
imageName=handles.jobs(projNum).imagename;
FirstImageNum=handles.jobs(projNum).firstimage;
LastImaNum=handles.jobs(projNum).lastimage;
ImageNamesList = handles.jobs(projNum).imagenameslist



cd(ImageDirectory);

name = char(ImageNamesList(FirstImaNum));
firstImg=imreadnd2(name,0,handles.jobs(projNum).intensityMax);

name = char(ImageNamesList(LastImaNum)); 
lastImg=imreadnd2(name,0,handles.jobs(projNum).intensityMax);

[img_h,img_w]=size(firstImg);



%-Minimal size of nuclei (1)
if jobdeterminer == 1
    figure, imshow(firstImg),title('Make a polygon around the smallest nucleus, by clicking on its rimm, going around clock- or anticlockwise,then press ENTER') ; 
    BW =roipoly;
    close
    if ~length(find(BW))==0
        handles.jobs(projNum).minsize= length(find(BW));
        change=1;
    end
end



%-Maximal size of nuclei (2)
if jobdeterminer == 2
    figure, imshow(firstImg), title('Make a polygon around the biggest nucleus, by clicking on its rimm, going around clock- or anticlockwise,then press ENTER') ;
     BW =roipoly;
    close
    if ~length(find(BW))==0
        handles.jobs(projNum).maxsize= length(find(BW));
        change=1;
    end
   
end
  


%-Minimal distance between two nuclei (3)
if jobdeterminer == 3
    figure, imshow(firstImg), title('click on the centers of the two cells closest to each other, then press ENTER') ;
     [x,y] =getpts;
     close
     
     if length(x)>1
         handles.jobs(projNum).minsdist= round(sqrt((x(1)-x(2))^2+(y(1)-y(2))^2));
         change=1;
    end
    
end


%change equals one, if anything has been specified by the user
if change==1
	% Update handles structure
	guidata(hObject, handles);
	
	%%%%%%%%save altered values to disk%%%%%%%%%%%%
	cd(handles.jobs(projNum).savedirectory)
	jobvalues=handles.jobs(projNum);
	save ('jobvalues','jobvalues')
	clear jobvalues
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%now we update the value in PolyTrack
	fillFields(handles,handles.jobs(projNum));
end