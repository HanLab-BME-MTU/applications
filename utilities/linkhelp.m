function linkhelp
% linkhelp finds out which cell is currently calling back and writes the
% details of this cell into the linklist (listbox)
%
% This is the callback of the linkbutton (manrelink). It adds the current
% callback cell to the linklist, so that it can later be selected for linkage
%
% SYNOPSIS       linkhelp
%
% INPUT          none (gets data from the current object, which is the object with 
%                      the callback manrelink (text made by changeframe))
%
% OUTPUT         none
%
% DEPENDENCIES   linkhelp uses { nothing }
%                                  
%                linkhelp is used by { manrelink (as a callback) }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source

% Fetch the current handles structure by finding hObject
hObject = findall (0, 'Tag', 'GUI_pp_jobbrowse_pb');
handles = guidata (hObject);

% Get the current list of cells to be linked
cellLinkerList = findall (0, 'Style', 'listbox', 'Tag', 'linklist');
selectedCell = get (cellLinkerList, 'UserData');
linkedList = handles.listofcells;

% Get the frame counter handle and fetch frame number
frameCounterHandle = findall (0, 'Style', 'text', 'Tag', 'picturecount');
frameCount = get (frameCounterHandle, 'String');

% Create a string containing the cell and frame numbers. This string
% will be writen into the linked list. Later we can extract this data for further processing
linkedCell = [selectedCell, '/', frameCount];


% add string to the existing list
if ~iscell (linkedList)
   linkedList = cellstr (linkedCell);
else
   if isempty (find (ismember (linkedList, linkedCell)))
      linkedList (end+1) = cellstr (linkedCell);
   else
      fprintf (1, 'You are trying to link the same cell twice. This link will be ignored.\n');
   end
end

% Store the list info in the handles structure
handles.listofcells = linkedList;

% Find and get rid of all the associated GUI objects as well
deleteButtonHandle = findall (0, 'Style', 'pushbutton', 'Tag', 'chuckout');
linkButtonHandle = findall (0, 'Style', 'pushbutton', 'Tag', 'linkbutton');
delete (deleteButtonHandle)
delete (cellLinkerList)
delete (linkButtonHandle)

% Update the handles structure
guidata (hObject, handles);

