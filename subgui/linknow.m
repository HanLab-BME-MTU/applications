function linknow
% linknow  finds out which cell is currently calling back and links it to
%          the cell selected in linklist (listbox)
%
% SYNOPSIS       linknow
%
% INPUT          none (gets data from the current object, which is the object with the callback manrelink (text made by changeframe))
%
% OUTPUT         none
%
% DEPENDENCIES   linknow uses {nothing}
%                                  
%                linknow is used by { manrelink (as a callback) }
%

%this is the callback of a click within the listbox, created by manrelink.
%Here we link the cell clicked on in the listbox to the current callback
%cell.

hObject=findall(0,'Tag','GUI_pp_jobbrowse_pb');
handles=guidata(hObject);

linklisthelpingH= findall(0,'Style','listbox','Tag','linklist');

%make sure there actually are cells listed in the listbox. If not we no
%cell to link to.
listofcells=handles.listofcells;
if ~iscell(listofcells)
    return
end

%Number of the current callback cell (the one the user has clicked on
%last). NOT on the cell selected within listbox!
%We have saved that number in userdata while creating listbox
celltoLink=str2num(get(linklisthelpingH,'UserData'));

%the linkage will take place in the current frame, not the one where the
%cell within the listbox was selected
piccount =findall(0,'Style','text','Tag','picturecount');
frametoLink=str2num(get(piccount,'String'));


%details on the cell we wish to have an other cell asigned to. This is the
%cell listed in the listbox, not the current callback cell.
%find out which cell of the list was selected
whichinlist=get(linklisthelpingH,'Value');

%get the string from thre listbox. That string contains the cell number and the
%frame in which it was added to the listbox, seprated by a '/'
FrameandCell=char(listofcells(whichinlist));
sepindex=findstr(FrameandCell,'/');

%isolate cell number and frame number
celltoKeep=str2num(FrameandCell(1:sepindex-1));
frametoKeep=str2num(FrameandCell(sepindex+1:end));


%the numbering for frame and for processed frame can differ, since not every
%existing frame had to be processed. We need to calculate which processed
%frame the current frame corresponds to.
mpmframe=((frametoLink-handles.jobvalues.firstimage)/handles.jobvalues.increment)+1

%link. (mpmframe*2 because for every processed frame there two columns of
%coordinates
handles.MPM(celltoKeep,mpmframe*2-1:end)=handles.MPM(celltoLink,mpmframe*2-1:end);

%erase redundant track
keep = questdlg('Shall the cell, which is about to be linked be erased completely, or shall its track previous to linking be kept? ','','ERASE!!!','keep','ERASE!!!');
if strcmp(keep,'keep')
     handles.MPM(celltoLink,mpmframe*2-1:end)=0;
    
else
     handles.MPM(celltoLink,:)=0;
end


%erase the linked cell from the linklist
if length(listofcells)==1
    listofcells = char('...');
    none=1
else
    listofcells(whichinlist) = [];
end

handles.listofcells=listofcells;


%save results
guidata(hObject,handles);


chuckhelpH =findall(0,'Style','pushbutton','Tag','chuckout');
linkbuttonH =findall(0,'Style','pushbutton','Tag','linkbutton');

delete(chuckhelpH)
delete(linklisthelpingH)
delete(linkbuttonH)
