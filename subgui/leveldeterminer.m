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
firstimg=imreadnd2(name,0,handles.jobs(projNum).intensityMax);


[img_h,img_w]=size(firstimg);

indxStr=sprintf(strg,last);
name=[main indxStr ext]; 
    
lastimg=imreadnd2(name,0,handles.jobs(projNum).intensityMax);



%get the intensity values of these pictures

%first picture
%background
backfirst=figure, imshow(firstimg,[]), title('Click on the background (approx 8 times). Make sure your clicks are spread out evenly. Then press enter') ;

[X,Y] = getpts(backfirst);
intense=[];
for dots=1:size(X,1)
    
	if X(dots) > 3 & X(dots) <img_w-3 & Y(dots) > 3 & Y(dots) <img_h-3  
      
        intense(end+1)=sum(sum(firstimg(Y(dots)-3:Y(dots)+3,X(dots)-3:X(dots)+3)))/49;
    end
end
handles.jobs(projNum).fi_background= sum(intense)/length(intense);
close
clear X
clear Y

%nucloi
nucfirst=figure, imshow(firstimg,[]), title('Click on the nucleoi (approx 8 times). Make sure your clicks on a lot of different ones. Then press enter');
[X,Y] = getpts(nucfirst);
intense=[];
for dots=1:size(X,1)
    
	if X(dots) > 3 & X(dots) <img_w-3 & Y(dots) > 3 & Y(dots) <img_h-3  
      
        intense(end+1)=sum(sum(firstimg(Y(dots)-3:Y(dots)+3,X(dots)-3:X(dots)+3)))/49;
    end
end
handles.jobs(projNum).fi_nucleus= sum(intense)/length(intense);
close
clear X
clear Y

%halos
halofirst=figure, imshow(firstimg,[]), title('Click on the halos(approx 8 times). Make sure your clicks on a lot of different ones. Then press enter');

[X,Y] = getpts(halofirst);
intense=[];
for dots=1:size(X,1)
    
	if X(dots) > 3 & X(dots) <img_w-3 & Y(dots) > 3 & Y(dots) <img_h-3  
      
        intense(end+1)=sum(sum(firstimg(Y(dots)-3:Y(dots)+3,X(dots)-3:X(dots)+3)))/49;
    end
end
handles.jobs(projNum).fi_halolevel= sum(intense)/length(intense);
close
clear X
clear Y



%%%%%%%%%last picture
%background
backlast=figure, imshow(lastimg,[]), title('Click on the background (approx 8 times). Make sure your clicks are spread out evenly. Then press enter') ;

[X,Y] = getpts(backlast);
intense=[];
for dots=1:size(X,1)
    
	if X(dots) > 3 & X(dots) <img_w-3 & Y(dots) > 3 & Y(dots) <img_h-3  
      
        intense(end+1)=sum(sum(lastimg(Y(dots)-3:Y(dots)+3,X(dots)-3:X(dots)+3)))/49;
    end
end
handles.jobs(projNum).la_background= sum(intense)/length(intense);
close
clear X
clear Y

%nucloi
nuclast=figure, imshow(lastimg,[]), title('Click on the nucleoi (approx 8 times). Make sure your clicks on a lot of different ones. Then press enter');

[X,Y] = getpts(nuclast);
intense=[];
for dots=1:size(X,1)
    
	if X(dots) > 3 & X(dots) <img_w-3 & Y(dots) > 3 & Y(dots) <img_h-3  
      
        intense(end+1)=sum(sum(lastimg(Y(dots)-3:Y(dots)+3,X(dots)-3:X(dots)+3)))/49;
    end
end
handles.jobs(projNum).la_nucleus= sum(intense)/length(intense);
close
clear X
clear Y

%halos
halolast=figure, imshow(lastimg,[]), title('Click on the halos (approx 8 times). Make sure your clicks on a lot of different ones. Then press enter');

[X,Y] = getpts(halolast);
intense=[];
for dots=1:size(X,1)
    
	if X(dots) > 3 & X(dots) <img_w-3 & Y(dots) > 3 & Y(dots) <img_h-3  
      
        intense(end+1)=sum(sum(lastimg(Y(dots)-3:Y(dots)+3,X(dots)-3:X(dots)+3)))/49;
    end
end
handles.jobs(projNum).la_halolevel= sum(intense)/length(intense);
close
clear X
clear Y




                  
% Update handles structure
guidata(hObject, handles);


