function linknow

hObject=findall(0,'Tag','GUI_pp_jobbrowse_pb');
handles=guidata(hObject);

linklisthelpingH= findall(0,'Style','listbox','Tag','linklist');

listofcells=handles.listofcells;
if ~iscell(listofcells)
    return
end

%details on the cell we wish to asign to an other cell
celltoLink=str2num(get(linklisthelpingH,'UserData'));

piccount =findall(0,'Style','text','Tag','picturecount');
frametoLink=str2num(get(piccount,'String'));


%details on the cell we wish to have an other cell asigned to
whichinlist=get(linklisthelpingH,'Value');

FrameandCell=char(listofcells(whichinlist));
sepindex=findstr(FrameandCell,'/');

celltoKeep=str2num(FrameandCell(1:sepindex-1));
frametoKeep=str2num(FrameandCell(sepindex+1:end));



handles.MPM(celltoKeep,frametoLink*2-1:end)=handles.MPM(celltoLink,frametoLink*2-1:end);

keep = questdlg('Shall the cell, which is about to be linked be erased completely, or shall its track previous to linking be kept? ','','ERASE!!!','keep','ERASE!!!');
if strcmp(keep,'keep')
     handles.MPM(celltoLink,frametoLink*2-1:end)=0;
    
else
     handles.MPM(celltoLink,:)=0;
end



if length(listofcells)==1
    listofcells = char('...');
    none=1
else
    listofcells(whichinlist) = [];
end

handles.listofcells=listofcells;

guidata(hObject,handles);


chuckhelpH =findall(0,'Style','pushbutton','Tag','chuckout');
linkbuttonH =findall(0,'Style','pushbutton','Tag','linkbutton');

delete(chuckhelpH)
delete(linklisthelpingH)
delete(linkbuttonH)
