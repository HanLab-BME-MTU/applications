function leveldeterminer(hObject)
% leveldeterminer let's the user specify values of images
%
% SYNOPSIS       leveldeterminer(hObject)
%
% INPUT          hObject : handle to an object of GUI (PolyTrack)

% OUTPUT         none (get written directly into handles
%                a level for each of the following:
%                background, nucloi, halos
%                each of those for both the first and last image 
%
% DEPENDENCIES   leveldeterminer uses {nothing}
%                                  
%                leveldeterminer is used by { PolyTrack }
%
% REMARK         leveldeterminer fetches directly in GUI PolyTrack:
% 					projNum : currently treated job
%                   handles : structure with information used within GUI
% 				  from handles.jobs(projNum):
% 					imagedirectory : where are the images 
% 					imagenameslist : list of images within imagedirectory with imagename
% 					firstimage : which images shall we start with (refers to imagenameslist)
% 					lastimage : which images shall be the last one (refers to imagenameslist)
% 					intensityMax : highest value image can have (calc from bitdepth)
% 					
%
% Colin Glass, Feb 04         



%This programm allows the user to specify the intensity values of his
%images interactively. Values of interest are: background, nuclei, halos
%each value get's specified for the first and the last frame of the current
%analysis, so that for every frame the values can be interpolated

handles = guidata(hObject);

projNum = get(handles.GUI_st_job_lb,'Value');
ImageDirectory = handles.jobs(projNum).imagedirectory;
FirstImaNum = handles.jobs(projNum).firstimage;
LastImaNum = handles.jobs(projNum).lastimage;
ImageNamesList = handles.jobs(projNum).imagenameslist




cd(ImageDirectory);

name = char(ImageNamesList(FirstImaNum));
firstImage=imreadnd2(name,0,handles.jobs(projNum).intensityMax);

name = char(ImageNamesList(LastImaNum));
lastImage=imreadnd2(name,0,handles.jobs(projNum).intensityMax);

[img_h,img_w]=size(firstImage);


%get the intensity values of these pictures
%%%%%%%%%%%%%%first picture%%%%%%%%%%%%%%%%%%%%%%%%%%

%background
backfirst = figure, imshow(firstImage,[]), title('Click on the background (approx 8 times). Make sure your clicks are spread out evenly. Right-click or press enter to finish.') ;

[X,Y] = getpts(backfirst);

intense = [];
for dots = 1:size(X,1)
    
	if X(dots) > 3 & X(dots) <img_w-3 & Y(dots) > 3 & Y(dots) <img_h-3  
      
        intense(end+1)=sum(sum(firstImage(Y(dots)-3:Y(dots)+3,X(dots)-3:X(dots)+3)))/49;
    end
end
handles.jobs(projNum).fi_background= sum(intense)/length(intense);
close
clear X
clear Y



%nucloi
nucfirst=figure, imshow(firstImage,[]), title('Click on the nucleoi (approx 8 times). Make sure your clicks on a lot of different ones. Right-click or press enter to finish.');
[X,Y] = getpts(nucfirst);
intense=[];
for dots=1:size(X,1)
    
	if X(dots) > 3 & X(dots) <img_w-3 & Y(dots) > 3 & Y(dots) <img_h-3  
      
        intense(end+1)=sum(sum(firstImage(Y(dots)-3:Y(dots)+3,X(dots)-3:X(dots)+3)))/49;
    end
end
handles.jobs(projNum).fi_nucleus= sum(intense)/length(intense);
close
clear X
clear Y



%halos
halofirst=figure, imshow(firstImage,[]), title('Click on the halos(approx 8 times). Make sure your clicks on a lot of different ones. Right-click or press enter to finish.');

[X,Y] = getpts(halofirst);
intense=[];
for dots=1:size(X,1)
    
	if X(dots) > 3 & X(dots) <img_w-3 & Y(dots) > 3 & Y(dots) <img_h-3  
      
        intense(end+1)=sum(sum(firstImage(Y(dots)-3:Y(dots)+3,X(dots)-3:X(dots)+3)))/49;
    end
end
handles.jobs(projNum).fi_halolevel= sum(intense)/length(intense);
close
clear X
clear Y




%%%%%%%%%last picture%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%background
backlast=figure, imshow(lastImage,[]), title('Click on the background (approx 8 times). Make sure your clicks are spread out evenly. Right-click or press enter to finish.') ;



[X,Y] = getpts(backlast);
intense=[];
for dots=1:size(X,1)
    
	if X(dots) > 3 & X(dots) <img_w-3 & Y(dots) > 3 & Y(dots) <img_h-3  
      
        intense(end+1)=sum(sum(lastImage(Y(dots)-3:Y(dots)+3,X(dots)-3:X(dots)+3)))/49;
    end
end
handles.jobs(projNum).la_background= sum(intense)/length(intense);
close
clear X
clear Y



%nucloi
nuclast=figure, imshow(lastImage,[]), title('Click on the nucleoi (approx 8 times). Make sure your clicks on a lot of different ones. Right-click or press enter to finish.');

[X,Y] = getpts(nuclast);
intense=[];
for dots=1:size(X,1)
    
	if X(dots) > 3 & X(dots) <img_w-3 & Y(dots) > 3 & Y(dots) <img_h-3  
      
        intense(end+1)=sum(sum(lastImage(Y(dots)-3:Y(dots)+3,X(dots)-3:X(dots)+3)))/49;
    end
end
handles.jobs(projNum).la_nucleus= sum(intense)/length(intense);
close
clear X
clear Y



%halos
halolast=figure, imshow(lastImage,[]), title('Click on the halos (approx 8 times). Make sure your clicks on a lot of different ones. Right-click or press enter to finish.');

[X,Y] = getpts(halolast);
intense=[];
for dots=1:size(X,1)
    
	if X(dots) > 3 & X(dots) <img_w-3 & Y(dots) > 3 & Y(dots) <img_h-3  
      
        intense(end+1) = sum( sum( lastImage(Y(dots)-3:Y(dots)+3,X(dots)-3: X(dots)+3) ) ) / 49;
    end
end
handles.jobs(projNum).la_halolevel= sum(intense)/length(intense);
close
clear X
clear Y




                  
% Update handles structure
guidata(hObject, handles);


