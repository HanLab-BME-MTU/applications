function linkhelp


hObject=findall(0,'Tag','GUI_pp_jobbrowse_pb');
handles=guidata(hObject);


linklisthelpingH =findall(0,'Style','listbox','Tag','linklist');
whichcell= get(linklisthelpingH,'UserData');

listofcells=handles.listofcells;


piccount =findall(0,'Style','text','Tag','picturecount');
whichpic =get(piccount,'String');
detailslink=[whichcell,'/',whichpic];


if ~iscell(listofcells)
    listofcells = cellstr(detailslink);
else
    if isempty(find(ismember(listofcells,detailslink)))
         listofcells(end+1) = cellstr(detailslink);
    else
        
        silly = questdlg('Why are you trying to link the same cell twice at the same time?','','I am tired','I am silly','I am silly');
	
    end
end


handles.listofcells=listofcells;

guidata(hObject,handles);

chuckhelpH =findall(0,'Style','pushbutton','Tag','chuckout');
linkbuttonH =findall(0,'Style','pushbutton','Tag','linkbutton');

delete(chuckhelpH)
delete(linklisthelpingH)
delete(linkbuttonH)