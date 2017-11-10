function linknow
% linknow  finds out which cell is currently calling back and links it to
%          the cell selected in linklist (listbox)
%
% SYNOPSIS       linknow
%
% INPUT          none (gets data from the current object, which is the object with 
%                      the callback manrelink (text made by ptShowSlidingFrames))
%
% OUTPUT         none
%
% DEPENDENCIES   linknow uses { nothing }
%                                  
%                linknow is used by { manrelink (as a callback) }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source

% Fetch the current handles structure by finding hObject
hObject = findall (0, 'Tag', 'GUI_pp_jobbrowse_pb');
handles = guidata (hObject);

% Make sure there actually are cells listed in the list. If not there's nothing to link 
% to so just return
linkedList = handles.listofcells;
if ~iscell(linkedList)
   fprintf (1, 'No cells present in the linked list. Returning...\n');
   return;
end

% Get the number of the cell the user has clicked on.
% That number was saved in userdata while creating the list
linkedListHandle = findall (0, 'Style', 'listbox', 'Tag', 'linklist');
cellToLink = str2num (get (linkedListHandle, 'UserData'));

% The linking will take place in the current frame, not the one where the
% cell within the listbox was selected
frameCounterHandle = findall (0, 'Style', 'text', 'Tag', 'picturecount');
frameCount = str2num (get (frameCounterHandle, 'String'));

% Details on the cell we wish to have another cell assigned to. This is the
% cell listed in the listbox, not the current callback cell.
selectedListCell = get (linkedListHandle, 'Value');

% Get the string from the listbox. That string contains the cell number and the
% frame in which it was added to the listbox, separated by a '/'
linkString = char (linkedList (selectedListCell));
separatorIndex = findstr (linkString, '/');

% Isolate cell number and frame number
linkToCell = str2num (linkString (1:separatorIndex-1));
linkToFrame = str2num (linkString (separatorIndex+1:end));

% The numbering for frame and for processed frame can differ, since not every
% existing frame had to be processed. We need to calculate which processed
% frame the current frame corresponds to.
mpmframe = ((frameCount - handles.jobvalues.firstimage) / handles.jobvalues.increment) + 1;

% Link mpmframe*2 because for every processed frame there two columns of coordinates
handles.MPM (linkToCell, mpmframe * 2 - 1:end) = handles.MPM (cellToLink, mpmframe * 2 - 1:end);

% Erase redundant track or not? Ask the user. Erasing it will be the default.
answer = questdlg ('Should the previous cell track be kept or can it be erased? ', '', 'Erase', 'Keep', 'Erase');
if strcmp (answer, 'Keep')
   % Keep the previous track
   handles.MPM(cellToLink, mpmframe * 2 - 1:end) = 0;
else
   % Erase it from MPM
   handles.MPM (cellToLink,:) = 0;
end

% Erase the linked cell from the linked list
if length (linkedList) == 1
   % The cell was the only entry so initialize the whole list (to empty)
   linkedList = char ('...');
else
   % Erase the cell entry
   linkedList (selectedListCell) = [];
end

% Store linked list for later use
handles.listofcells=linkedList;

% Find and get rid of all the associated GUI objects as well
deleteButtonHandle = findall (0, 'Style', 'pushbutton', 'Tag', 'chuckout');
linkButtonHandle = findall (0, 'Style', 'pushbutton', 'Tag', 'linkbutton');
delete (linkButtonHandle)
delete (deleteButtonHandle)
delete (linkedListHandle)

% Update the handles structure
guidata(hObject,handles);

