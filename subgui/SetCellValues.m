function SetCellValues(hObject,jobdeterminer)
% SetCellValues let's the user specify values of images
%
% SYNOPSIS       SetCellValues(hObject,jobdeterminer)
%
% INPUT          hObject : handle to an object of GUI (PolyTrack)
%                jobdeterminer : what kind of value interestes the user 
%                                1 - minimal size nucloi
%                                2 - maximal size nucloi
%                                3 - minimal distance between two nuclei
%
% OUTPUT         none (get written directly into handles
%                minsize or maxsize or minsdist
%
% DEPENDENCIES   SetCellValues uses {nothing}
%                                  
%                SetCellValues is used by { PolyTrack }
%
% REMARK         SetCellValues fetches directly in GUI PolyTrack:
% 					projNum : currently treated job
%                   handles : structure with information used within GUI
% 			      from  handles.jobs(projNum):
% 					imagedirectory : where are the images 
% 					imagename : what are the images called
% 					imagenameslist : list of images within imagedirectory with imagename
% 					firstimage : which images shall we start with (refers to imagenameslist)
% 					lastimage : which images shall be the last one (refers to imagenameslist)
% 					intensityMax : highest value image can have (calc from bitdepth)
% 					
%
% Colin Glass, Feb 04         


%this routine enables the user to specify:
%-Minimal size of nuclei (1)
%-Maximal size of nuclei (2)
%-Minimal distance between two nuclei (3)
%interactively. Which one of the three it is going to be, depends on the
%value of jobdeterminer

change=0;
handles = guidata(hObject);

projNum = get(handles.GUI_st_job_lb,'Value');
ImageDirectory=handles.jobs(projNum).imagedirectory;
ImageName=handles.jobs(projNum).imagename;
FirstImaNum=handles.jobs(projNum).firstimage;
LastImaNum=handles.jobs(projNum).lastimage;
ImageNamesList = handles.jobs(projNum).imagenameslist;
intensityMax = handles.jobs(projNum).intensityMax;


cd(ImageDirectory);

name = char(ImageNamesList(FirstImaNum));
firstImg=imreadnd2(name,0,intensityMax);

name = char(ImageNamesList(LastImaNum)); 
lastImg=imreadnd2(name,0,intensityMax);

[img_h,img_w]=size(firstImg);



%-Minimal size of nuclei (1)
if jobdeterminer == 1
    figure, imshow(firstImg,[]),title('Make a polygon around the smallest nucleus, by clicking on its rimm, going around clock- or anticlockwise,then press ENTER') ; 
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