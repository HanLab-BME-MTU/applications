function manrelink
% manrelink creates objects with callbacks to manipulate cell tracks
%
% This function creates three objects:
% - listbox (linklist) where the program keeps track of the cells you have selected for
%   linking together. When you click on one of the listed cells, it will be linked to
%   the current callback cell
% - pushbutton (linkbutton). If you click on this, the current callback cell
%   will be added to the list (in listbox)
% - pushbutton (chuckout). If you click on this, the current callback cell
%   will be erased
%
% SYNOPSIS       manrelink
%
% INPUT          none (gets data from the current object, which is the object with 
%                      the callback manrelink (text made by changeframe))
%
% OUTPUT         none
%
% DEPENDENCIES   manrelink uses { nothing }
%                                  
%                manrelink is used by { changeframe (as a callback) }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source

% Fetch the current handles structure by finding hObject
hObject = findall (0, 'Tag', 'GUI_pp_jobbrowse_pb');
handles = guidata (hObject);

% First we have to find out which cell the user has clicked on
selectedCell = get (gco,'String')

% First we look if these three objects already exist somewhere and, if that
% is the case, erase them.
linkedListHandle = findall (0, 'Style', 'listbox', 'Tag', 'linklist');
if ~isempty (linkedListHandle)
   delete (linkedListHandle);
end;

linkButtonHandle = findall (0, 'Style', 'pushbutton', 'Tag', 'linkbutton');
if ~isempty (linkButtonHandle)
   delete (linkButtonHandle);
end;

deleteButtonHandle = findall (0, 'Style', 'pushbutton', 'Tag', 'chuckout');
if ~isempty (deleteButtonHandle)
   delete (deleteButtonHandle);
end;

% Get the current list that should go into listbox
linkedList = handles.listofcells;

% Find the coordinates of the pixel the user clicked on, so that we can create the objects there
% If this is not on a cell, gcbf returns [] and nothing happens
figure (gcbf);
coordinates = get (gcf, 'currentPoint')
%imgPos = get (gca, 'Position');

% Calculate the coordinates the objects should have
objectCoord = [coordinates(1) coordinates(2) 60 60];
relinkButtonCoordinates = [coordinates(1) coordinates(2)+60 60 30];
deleteButtonCoordinates = [coordinates(1) coordinates(2)+90 60 30];

% Create the gui objects including their respective callbacks
linkedListHandle = uicontrol ('Style', 'listbox', ...
                              'Callback', 'linknow', ...
                              'String', linkedList, ...
                              'UserData', selectedCell, ...
                              'Tag', 'linklist', ...
                              'Position', objectCoord);

linkButtonHandle = uicontrol ('Style', 'pushbutton', 'String', 'Link', ...
                              'Callback', 'linkhelp', ...
                              'Tag', 'linkbutton', ...
                              'UserData', selectedCell, ...
                              'Position', relinkButtonCoordinates);

deleteButtonHandle = uicontrol ('Style', 'pushbutton', 'String', 'Delete', ...
                                'Callback', 'deletingrow', ...
                                'Tag', 'chuckout', ...
                                'UserData', selectedCell, ...
                                'Position', deleteButtonCoordinates);

