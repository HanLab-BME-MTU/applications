function ptShowNextInSequenceForDetect
% ptShowNextInSequenceForDetect finds the rigth image and coordinates. These it shows in the
% figure opened by ptShowImageSequenceForDetect
%
% SYNOPSIS       ptShowNextInSequenceForDetect
%
% INPUT          none (it gets values from the slider created in ptShowImageSequenceForDetect)
%
% OUTPUT         none (it updates the handles object directly)
%
% DEPENDENCIES   ptShowNextInSequenceForDetect uses {nothing}
%                                  
%                ptShowNextInSequenceForDetect is used by { ptShowImageSequenceForDetect }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Jun 29          First Release

% get the structure in the subfunction
handles = guidata(gcbo);

% Fetch the values from the handlesMt
imageDirectory = handles.imageDirectory;
imageName = handles.imageName;
firstImage = handles.firstImage;
lastImage = handles.lastImage;
imageNameList = handles.imageNameList;
intensityMax = handles.intensityMax;
imageRange = handles.imageRange;
saveDirectory = handles.saveDirectory;
bodyName = handles.bodyName;

% Delete the axes currently shown in the figure on the screen
delete (gca)

% Get a handle to the main window
hManDetect = findall(0,'Tag','ptManDetect','Name','ptManDetect');

if isempty(hManDetect)
    % cannot find: return
    msgStr = ['Error: ptManTrack GUI cannot be found. Exiting...'];
    h = errordlg(msgStr);
    uiwait(h);
    
    hSeqFigure = findall(0,'Tag','ptSequenceFigure');
    if ~isempty(hSeqFigure)
        close(hSeqFigure);
    end
    
    return;
else    
    % Get the mantrack handle struct
    handlesMt = guidata(hManDetect);

    % Look for objects needed for information
    sliderHandle = findall (0, 'Tag', 'imageSlider');
    %hObject = findall (0, 'Tag', 'GUI_add_pb');
    imageCounterHandle = findall (0, 'Style', 'text', 'Tag', 'imageCounter');

    % Get the current value of the slider, so that we know which frame the user wants to process
    sliderValue = get (sliderHandle, 'Value');
    sliderValue = round (sliderValue * imageRange);

    % Calculate the frame number to show
    %imageNumber = (sliderValue - 1)+firstImage;
    imageNumber = (sliderValue - 1)+firstImage;

    % Write the current frame number in the little window above the slider
    set (imageCounterHandle, 'String', num2str (imageNumber));

    % Read the image frame from disk
    cd (imageDirectory);
    fileName = char(imageNameList(imageNumber));
    image = imreadnd2(fileName, 0, intensityMax);

    % Show the frame on the screen in the current figure
    hold on;
    imshow (image, []), title (num2str (imageNumber));
    hold off;

    % Make sure the axis are correct
    axis([1 size(image,2) 1 size(image,1)]);
    
    % Get points from the image
    h = findall(0,'Tag','ptSequenceFigure');
    [x,y] = getpts(h);
    coords = [x y];
    roundCoords = round(coords);

    % Save if any points were clicked
    if ~isempty(coords)
        save([saveDirectory filesep 'manDetectResults_' bodyName(1:end-2) '_' ...
             num2str(imageNumber) '.mat'], 'roundCoords');
        save([saveDirectory filesep 'manDetectResults_' bodyName(1:end-2) '_' ...
             num2str(imageNumber) '_coord.mat'], 'coords'); 
    end

    % That's it: wait for the next user action
end
