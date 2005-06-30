function ptShowImageSequenceForDetect()
% ptShowImageSequenceForDetect opens a figure and puts a slider into it; the callback belonging
% to that slider (ptShowNextInSequenceForDetect) takes care of showing the correct frame
%
% SYNOPSIS       ptShowImageSequenceForDetect ()
%
% INPUT          none
%
% OUTPUT         none
%
% DEPENDENCIES   ptShowImageSequenceForDetect uses {nothing}
%                                  
%                ptShowImageSequenceForDetect is used by { ptManDetect }
%
% Revision History
% Name                  Date            Comment
% --------------------------------------------------------
% Andre Kerstens        Jun 29          First version

% Get a handle to the main window
hManDetect = findall(0,'Tag','ptManDetect','Name','ptManDetect');

if isempty(hManDetect)
    % cannot find: return
    msgStr = ['Error: ptManTrack GUI cannot be found. Exiting...'];
    h = errordlg(msgStr);
    uiwait(h);
    return;
else    
    % Get the mandetect handle struct
    handlesMt = guidata(hManDetect);

    % Fetch the jobvalues and image directory
    imageDirectory = handlesMt.movieDirectory;
    imageName = handlesMt.imageName;
    firstImage = handlesMt.firstImage;
    lastImage = handlesMt.lastImage;
    imageNameList = handlesMt.imageNameList;
    intensityMax = handlesMt.intensityMax;
    saveDirectory = get(handlesMt.edit_savedir_ptmt,'String');
    bodyName = handlesMt.bodyName;
           
    % if the slider (and less important, the little windos with the number of
    % the current frame) already exist, delete them now.
    imageCounterHandle = findall (0, 'Style', 'text', 'Tag', 'imageCounter');
    sliderHandle = findall (0, 'Style', 'slider', 'Tag', 'imageSlider');

    % Calculate the image range taking into account the increment between frames
    %imageRange = size (validFrames,2);
    imageRange = lastImage - firstImage + 1;

    % Generate the slider step values for the uicontrol
    % First the arrow slide step (1):
    slider_step(1) = 1 / imageRange;
    % Then the stepsize (5):
    slider_step(2) = 5 / imageRange;

    % Draw a new figure on the screen
    h = figure;
    %set(h,'Tag','ptSequenceFigure','CloseRequestFcn',@ptManDetectClose);
    set(h,'Tag','ptSequenceFigure');

    % Draw the frame counter in the figure; it is identified by the tag
    % imageCount
    imageCounterHandle = uicontrol('Style', 'text',...
                                   'Units', 'normalized',...
                                   'Tag', 'imageCounter',...
                                   'Position', [0.02,0.93,0.05,0.06]);

    % Set the frame counter to the first image number
    set (imageCounterHandle, 'String', num2str(firstImage));

    % Draw the slider in the figure; it is identified by the tag imageSlider and calls
    % the function ptShowSlidingFrames when moved
    sliderHandle = uicontrol('Style', 'slider', ...
                             'Units', 'normalized', ... 
                             'Value', 1/(imageRange), ...
                             'Min', 1/(imageRange), ...
                             'Max', 1, ...
                             'SliderStep', slider_step, ...
                             'Callback', 'ptShowNextInSequenceForDetect', ...
                             'Tag', 'imageSlider', ...
                             'Position', [0.02,0.02,0.05,0.9]);

    % Now we show the first image including red dots for the cells overlaid on the image
    % Read the image frame from disk
    cd (imageDirectory);
    fileName = char(imageNameList(firstImage));
    image = imreadnd2(fileName, 0, intensityMax);

    % Get image size
    imInfo = imfinfo([imageDirectory filesep fileName]);
    rowSize = imInfo.Height
    colSize = imInfo.Width
    
    % Show the frame on the screen in the current figure
    hold on;
    imshow (image, []), title (num2str (firstImage));
    hold off;

    % Make sure the axis are correct
    axis([1 size(image,2) 1 size(image,1)]);

    % Get points from the image
    [x,y] = getpts(h);
    coords = [x y];
    roundCoords = round(coords);

    % Save if any points were clicked
    if ~isempty(coords)
        save([saveDirectory filesep 'manDetectResults_' bodyName(1:end-2) '_' ...
             num2str(firstImage) '.mat'], 'roundCoords');
    end
    
    % create structure of handles
    handlesNew = guihandles(sliderHandle); 
    
    % add relevant data
    handlesNew.imageDirectory   = imageDirectory;
    handlesNew.imageName        = imageName;
    handlesNew.firstImage       = firstImage;
    handlesNew.lastImage        = lastImage;
    handlesNew.imageNameList    = imageNameList;
    handlesNew.intensityMax     = intensityMax;
    handlesNew.imageRange       = imageRange;
    handlesNew.rowSize          = rowSize;
    handlesNew.colSize          = colSize;
    handlesNew.saveDirectory    = saveDirectory;
    handlesNew.bodyName         = bodyName;
    
    % save the structure
    guidata(sliderHandle,handlesNew)     
end