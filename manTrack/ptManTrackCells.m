function ptManTrackCells(src,eventdata)
% ptManTrackCells is the callback for the clickable cell nuclei
%
% SYNOPSIS       ptManTrackCells
%
% INPUT          none
%
% OUTPUT         none
%
% DEPENDENCIES   ptManTrackCells uses {nothing}
%                                  
%                ptManTrackCells is used by { ptShowImageSequence, ptShowNextInSequence }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Feb 05          First Release

% get a handle to the current axes
h = get(gca);

% Get the handles struct of the main gui
handles = guidata(gcbo);

% Fetch the values from the handlesMt
resultsDirectory = handles.resultsDirectory;
imageRange = handles.imageRange;
MPM = handles.MPM;
validFrames = handles.validFrames;
manualMPM = handles.manualMPM;
rowSize = handles.rowSize;
colSize = handles.colSize;

% Get the current frame number
imageCounterHandle = findall (0, 'Style', 'text', 'Tag', 'imageCounter');
imageNr = str2num(get(imageCounterHandle, 'String'));

% Get the coordinates clicked by the user
row = round(h.CurrentPoint(1,2));
col = round(h.CurrentPoint(1,1));

% Store the first frame where the user started clicking
if isempty(find(manualMPM))
    handles.firstTrackedFrame = imageNr;
end

% Store in an array for later
manualMPM(1,2*imageNr-1) = col;
manualMPM(1,2*imageNr) = row;

% Let's put a yellow dot on the coords clicked by the user
firstCoords = MPM(:,2*imageNr-1:2*imageNr);
realCellIndex = find(firstCoords(:,1) & firstCoords(:,2));
firstCoords = firstCoords(realCellIndex,:);

% Pick out the ones that are closed to the clicked point
[value,indx] = min(sqrt((firstCoords(:,1)-col).^2 + ...
                        (firstCoords(:,2)-row).^2));

% Plot the point yellow
hold on;
plot(firstCoords(indx,1),firstCoords(indx,2),'y.');
hold off;
                        
% Store updated manual MPM
handles.manualMPM = manualMPM;

% Update handles
guidata(gcbo,handles);