function deletingrow
% deletingrow  finds out which cell is currently calling this function and erase it from MPM
%
% SYNOPSIS       deletingrow
%
% INPUT          none (gets data from the current object, which is the object with the 
%                      callback manrelink (text made by changeframe))
%
% OUTPUT         none
%
% DEPENDENCIES   deletingrow uses { nothing }
%                                  
%                deletingrow is used by { manrelink (as a callback) }
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source

% Fetch the current handles structure by finding hObject
hObject = findall (0, 'Tag', 'GUI_pp_jobbrowse_pb');
handles = guidata (hObject);

% Get the information from the delete button
deleteButtonHandle = findall (0, 'Style', 'pushbutton', 'Tag', 'chuckout');

% Find out which cell the user wants to delete
selectedCell = str2num (get (deleteButtonHandle, 'UserData'));

% Delete the found cell from MPM
handles.MPM (selectedCell,:) = 0;

% Find and get rid of all the associated GUI objects as well
linkedListHandle = findall (0, 'Style', 'listbox', 'Tag', 'linklist');
linkButtonHandle = findall (0, 'Style', 'pushbutton', 'Tag', 'linkbutton');
delete (linkedListHandle)
delete (linkButtonHandle)
delete (deleteButtonHandle)

% Update the handles structure
guidata (hObject, handles);

