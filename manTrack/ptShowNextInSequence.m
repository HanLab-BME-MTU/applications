function ptShowNextInSequence
% ptShowNextInSequence finds the rigth image and coordinates. These it shows in the
% figure opened by ptShowImageSequence
%
% SYNOPSIS       ptShowNextInSequence
%
% INPUT          none (it gets values from the slider created in ptShowImageSequence)
%
% OUTPUT         none (it updates the handles object directly)
%
% DEPENDENCIES   ptShowNextInSequence uses {nothing}
%                                  
%                ptShowNextInSequence is used by { ptShowImageSequence }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Feb 05          First Release

% Delete the axes currently shown in the figure on the screen
delete (gca)

% Get a handle to the main window
hManTrack = findall(0,'Tag','ptManTrack','Name','ptManTrack');

if isempty(hManTrack)
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
    handlesMt = guidata(hManTrack);

    % Look for objects needed for information
    sliderHandle = findall (0, 'Tag', 'imageSlider');
    %hObject = findall (0, 'Tag', 'GUI_add_pb');
    imageCounterHandle = findall (0, 'Style', 'text', 'Tag', 'imageCounter');

    % Fetch the values from the handlesMt
    imageDirectory = handlesMt.imageDirectory;
    imageName = handlesMt.imageName;
    firstImage = handlesMt.firstImage;
    lastImage = handlesMt.lastImage;
    imageNameList = handlesMt.imageNameList;
    intensityMax = handlesMt.intensityMax;
    imageRange = handlesMt.imageRange;

    % Get the current value of the slider, so that we know which frame the user wants to process
    sliderValue = get (sliderHandle, 'Value');
    sliderValue = round (sliderValue * imageRange);

    % Calculate the frame number to show
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

    % Identify the real cells (at least one coord different from zero)
    realCellIndex = find (handlesMt.MPM(:, 2 * sliderValue - 1) | handlesMt.MPM(:, 2 * sliderValue));

    % Find the row indices from a transposed MPM matrix
    cellsWithNums = zeros (size (handlesMt.MPM, 1), 3);
    cellsWithNums(:,3) = [1:1:size(handlesMt.MPM,1)]';

    % Grab all rows in MPM, so that the row indices correspond to the cells
    cellsWithNums(:,1:2) = handlesMt.MPM(:, 2 * sliderValue - 1 : 2 * sliderValue);

    % Now take the cells identified as real cells (at least one coord different from zero)
    % and plot those as red dots. The cell number is written as colored text on the current axes.
    hold on;
    dots = plot (cellsWithNums (realCellIndex, 1), cellsWithNums (realCellIndex, 2), 'r.');
    set(dots,'Tag','dots','ButtonDownFcn','ptManTrackCells');
    hold off;
    
    % That's it: wait for the next user action
end
