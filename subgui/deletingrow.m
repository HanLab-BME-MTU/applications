function deletingrow

%this is the callback of the chuckout button (manrelink).
%it erases the current callback cell

hObject=findall(0,'Tag','GUI_pp_jobbrowse_pb');
handles=guidata(hObject);


chuckhelpH =findall(0,'Style','pushbutton','Tag','chuckout');

%find out which cell the user wishes to delete
whichcell= str2num(get(chuckhelpH,'UserData'));

%delete that cell
handles.MPM(whichcell,:)=0;

%save changed data
guidata(hObject,handles);



linklisthelpingH= findall(0,'Style','listbox','Tag','linklist');
linkbuttonH =findall(0,'Style','pushbutton','Tag','linkbutton');

delete(chuckhelpH)
delete(linklisthelpingH)
delete(linkbuttonH)
