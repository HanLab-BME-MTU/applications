function leveldeterminer(hObject)

handles = guidata(hObject);

projNum = get(handles.GUI_st_job_lb,'Value');
ImageDirectory = handles.jobs(projNum).imagedirectory;

FirstImaNum = handles.jobs(projNum).firstimage;
LastImaNum = handles.jobs(projNum).lastimage;
ImageNamesList = handles.jobs(projNum).imagenameslist

cd(ImageDirectory);

name = char(ImageNamesList(FirstImaNum));

firstImage=imreadnd2(name,0,handles.jobs(projNum).intensityMax);


[img_h,img_w]=size(firstImage);


name = char(ImageNamesList(LastImaNum));

lastImage=imreadnd2(name,0,handles.jobs(projNum).intensityMax);



%get the intensity values of these pictures

%first picture
%background
backfirst = figure, imshow(firstImage,[]), title('Click on the background (approx 8 times). Make sure your clicks are spread out evenly. Then press enter') ;

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
nucfirst=figure, imshow(firstImage,[]), title('Click on the nucleoi (approx 8 times). Make sure your clicks on a lot of different ones. Then press enter');
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
halofirst=figure, imshow(firstImage,[]), title('Click on the halos(approx 8 times). Make sure your clicks on a lot of different ones. Then press enter');

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



%%%%%%%%%last picture
%background
backlast=figure, imshow(lastImage,[]), title('Click on the background (approx 8 times). Make sure your clicks are spread out evenly. Then press enter') ;

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
nuclast=figure, imshow(lastImage,[]), title('Click on the nucleoi (approx 8 times). Make sure your clicks on a lot of different ones. Then press enter');

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
halolast=figure, imshow(lastImage,[]), title('Click on the halos (approx 8 times). Make sure your clicks on a lot of different ones. Then press enter');

[X,Y] = getpts(halolast);
intense=[];
for dots=1:size(X,1)
    
	if X(dots) > 3 & X(dots) <img_w-3 & Y(dots) > 3 & Y(dots) <img_h-3  
      
        intense(end+1)=sum(sum(lastImage(Y(dots)-3:Y(dots)+3,X(dots)-3:X(dots)+3)))/49;
    end
end
handles.jobs(projNum).la_halolevel= sum(intense)/length(intense);
close
clear X
clear Y




                  
% Update handles structure
guidata(hObject, handles);


