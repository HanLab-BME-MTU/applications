function cellselect
% cellselect saves a cell within handles as a selected cell for future use
%
% SYNOPSIS       cellselect
%
% INPUT          none (gets data from the current object, which is the object with the callback cellselect (text made by changeframe))
%
% OUTPUT         none (writes which cell has been selected into handles)
%
% DEPENDENCIES   cellselect uses {nothing}
%                                  
%                cellselect is used by { changeframe (as a callback) }
%
% Colin Glass, Feb 04     


hObject=findall(0,'Tag','GUI_pp_jobbrowse_pb');
handles=guidata(hObject);

%first we have to find out which cell the user has clicked on
whichcell=str2num(get(gco,'String'));
set(gco,'Color','y')

handles.selectedcells(end+1,1)=whichcell;

guidata(hObject, handles);