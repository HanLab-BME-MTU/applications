function deletingrow

hObject=findall(0,'Tag','GUI_pp_jobbrowse_pb');
handles=guidata(hObject);


chuckhelpH =findall(0,'Style','pushbutton','Tag','chuckout');

whichcell= str2num(get(chuckhelpH,'UserData'));


handles.MPM(whichcell,:)=0;

guidata(hObject,handles);



linklisthelpingH= findall(0,'Style','listbox','Tag','linklist');
linkbuttonH =findall(0,'Style','pushbutton','Tag','linkbutton');

delete(chuckhelpH)
delete(linklisthelpingH)
delete(linkbuttonH)
