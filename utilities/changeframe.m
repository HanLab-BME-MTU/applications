function changeframe
% changeframe finds the rigth image and coordinates. These it shows in the
%             figure opened by manuelpostpro
%
% SYNOPSIS       changeframe
%
% INPUT          none (it gets values from the slider created in
%                      manuelpostpro)
%
% OUTPUT         none
%
% DEPENDENCIES   changeframe uses {nothing}
%                                  
%                changeframe is used by { manuelpostpro (as a callback) }
%
% Colin Glass, Feb 04         


%this is the callback of the slider, created in manualpostpro
%What we do here is:
%- find out which frame the user currently wants to look at
%- show this frame in the figure (created in manualpostpro)
%- plot the coordinates of the cells into this picture
%- plot the cell numbers (near the respective cells) with the 
%  callback manrelink
%if the user clicks on the number of the cell, manrelink comes into play.
%Look there for more



%delete the figure currently shown in the figure
delete(gca)

%look for objects needed for information
pictureslideH=findall(0,'Tag','pictureslide');
hObject=findall(0,'Tag','GUI_pp_jobbrowse_pb');
piccount =findall(0,'Style','text','Tag','picturecount');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%fetch and extract values

handles=guidata(hObject);

%the good old jobvalues... we need them once more                
imagedirectory=handles.postpro.imagepath;
imagename=handles.jobvalues.imagename;
firstimage=handles.jobvalues.firstimage;
lastimage=handles.jobvalues.lastimage;
increment=handles.jobvalues.increment;
ma=handles.ma;
ImageNamesList = handles.jobvalues.imagenameslist;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%get the current value of the slider, so that we know which frame the user
%wants
slidervalue=get(pictureslideH,'Value');
slidervalue=round(slidervalue*(ma-1)+1);

%calculate the number of the frame
whichpic=(slidervalue-1)*increment + firstimage;

%write the current frame number in the little window above the slider
set(piccount,'String',num2str(whichpic));


%go get the frame
cd (imagedirectory);

name = char(ImageNamesList(whichpic));
picnew=imreadnd2(name,0,handles.jobvalues.intensityMax);


%show the frame
hold on;
imshow(picnew,[]), title(num2str(whichpic))

  
%get the cells corresponding to this frame. We create a third column
%(first two being [x,y]) with the row indices of MPM. We are NOT interested in
%the indices the cells will have in the vector cellsWithNums!!! Why?
%because in MPM every cell has it's own row, so the row indices of MPM is
%the actual number of the cell

%identify the real cells (at least one coord different from zero)
indRealCell=find(handles.MPM(:,2*slidervalue-1) | handles.MPM(:,2*slidervalue));

cellsWithNums=zeros(size(handles.MPM,1),3);

%row indices
cellsWithNums(:,3)=[1:1:size(handles.MPM,1)]';

%here we grab all rows in MPM, so that the row indices correspond to the cells
cellsWithNums(:,1:2)=handles.MPM(:,2*slidervalue-1:2*slidervalue);

%NOW we only take the cells identified as real cells (at least one coord
%different from zero)
plot(cellsWithNums(indRealCell,1),cellsWithNums(indRealCell,2),'r.');
txt= text(cellsWithNums(indRealCell,1),cellsWithNums(indRealCell,2),num2str(cellsWithNums(indRealCell,3)),'Color','r');


if  handles.whichcallback == 1
    set(txt,'ButtonDownFcn','manrelink'); 
   
elseif  handles.whichcallback == 2
    if ~isempty(handles.selectedcells)
        handles.selectedcells=[];
        guidata(hObject, handles);
    end
    set(txt,'ButtonDownFcn','cellselect');
    
end

hold off;



