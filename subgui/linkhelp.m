function linkhelp
% linkhelp finds out which cell is currently calling back and writes the
%          details of this cell into the linklist (listbox)
%
% SYNOPSIS       linkhelp
%
% INPUT          none (gets data from the current object, which is the object with the callback manrelink (text made by changeframe))
%
% OUTPUT         none
%
% DEPENDENCIES   linkhelp uses {nothing}
%                                  
%                linkhelp is used by { manrelink (as a callback) }
%
% Colin Glass, Feb 04   
%this is the callback of the linkbutton (manrelink). It adds the current
%callback cell to the linklist, so that it can later be selected for
%linkage

hObject=findall(0,'Tag','GUI_pp_jobbrowse_pb');
handles=guidata(hObject);

%get the current list of cells to be linked
linklisthelpingH =findall(0,'Style','listbox','Tag','linklist');
whichcell= get(linklisthelpingH,'UserData');
listofcells=handles.listofcells;


piccount =findall(0,'Style','text','Tag','picturecount');
whichpic =get(piccount,'String');
%make a string containing the cell number and the frame number. This string
%will be writen into the linklist. Later we can extract this data again
detailslink=[whichcell,'/',whichpic];


%add string to the existing list
if ~iscell(listofcells)
    listofcells = cellstr(detailslink);
else
    if isempty(find(ismember(listofcells,detailslink)))
         listofcells(end+1) = cellstr(detailslink);
    else
        
        silly = questdlg('Why are you trying to link the same cell twice at the same time?','','I am tired','I am silly','I am silly');
	
    end
end

%save changed data
handles.listofcells=listofcells;

guidata(hObject,handles);

chuckhelpH =findall(0,'Style','pushbutton','Tag','chuckout');
linkbuttonH =findall(0,'Style','pushbutton','Tag','linkbutton');

delete(chuckhelpH)
delete(linklisthelpingH)
delete(linkbuttonH)