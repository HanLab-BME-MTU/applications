function cellselect
% cellselect saves a cell within handles as a selected cell for future use
%
% SYNOPSIS       cellselect
%
% INPUT          none (gets data from the current object, which is the object 
%                      with the callback cellselect (text made by changeframe))
%
% OUTPUT         none (writes which cell has been selected into handles)
%
% DEPENDENCIES   cellselect uses { nothing }
%                                  
%                cellselect is used by { changeframe (as a callback) }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source

% Get the current handles structure by finding hObject
hObject = findall (0, 'Tag', 'GUI_pp_jobbrowse_pb');
handles = guidata (hObject);

% Find out which cell the user has clicked on
selectedCell = str2num (get (gco, 'String'));
set (gco, 'Color', 'y')

% And add it to the list of selected cells
handles.selectedcells (end + 1, 1) = selectedCell;

% Update the handles structure
guidata(hObject, handles);
