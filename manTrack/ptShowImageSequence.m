function ptShowImageSequence ()
% ptShowImageSequence opens a figure and puts a slider into it; the callback belonging
% to that slider (ptShowNextInSequence) takes care of showing the correct frame
%
% SYNOPSIS       ptShowImageSequence ()
%
% INPUT          none
%
% OUTPUT         none
%
% DEPENDENCIES   ptShowImageSequence uses {nothing}
%                                  
%                ptShowImageSequence is used by { ptManTrack }
%
% Revision History
% Name                  Date            Comment
% --------------------------------------------------------
% Andre Kerstens        Mar 04          Cleaned up source

% Get a handle to the main window
hManTrack = findall(0,'Tag','ptManTrack','Name','ptManTrack');

if isempty(hManTrack)
    % cannot find: return
    msgStr = ['Error: ptManTrack GUI cannot be found. Exiting...'];
    h = errordlg(msgStr);
    uiwait(h);
    return;
else    
    % Get the mantrack handle struct
    handlesMt = guidata(hManTrack);

    % Fetch the jobvalues and image directory
    imageDirectory = handlesMt.movieDirectory;
    imageName = handlesMt.imageName;
    firstImage = handlesMt.firstImage;
    lastImage = handlesMt.lastImage;
    imageNameList = handlesMt.imageNameList;
    intensityMax = handlesMt.intensityMax;
    resultsDirectory = handlesMt.resultsDirectory;
    MPM = handlesMt.MPM;
    validFrames = handlesMt.validFrames;
    saveDirectory = get(handlesMt.edit_savedir_ptmt,'String');
    bodyName = handlesMt.bodyName;
           
    % if the slider (and less important, the little windos with the number of
    % the current frame) already exist, delete them now.
    imageCounterHandle = findall (0, 'Style', 'text', 'Tag', 'imageCounter');
    sliderHandle = findall (0, 'Style', 'slider', 'Tag', 'imageSlider');

    % Calculate the image range taking into account the increment between frames
    %imageRange = lastImage - firstImage + 1;
    imageRange = size (validFrames,2);
    %handlesMt.imageRange = imageRange;

    % Generate the slider step values for the uicontrol
    % First the arrow slide step (1):
    slider_step(1) = 1 / imageRange;
    % Then the stepsize (5):
    slider_step(2) = 5 / imageRange;

    % Update the handlesMt structure
    %guidata(hManTrack, handlesMt);

    % Draw a new figure on the screen
    h = figure;
    set(h,'Tag','ptSequenceFigure','CloseRequestFcn',@ptManTrackClose);

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
                             'Callback', 'ptShowNextInSequence', ...
                             'Tag', 'imageSlider', ...
                             'Position', [0.02,0.02,0.05,0.9]);

    % Now we show the first image including red dots for the cells overlaid on the image
    % Read the image frame from disk
    cd (imageDirectory);
    fileName = char(imageNameList(firstImage));
    image = imreadnd2(fileName, 0, intensityMax);

    % Get image size
    imInfo = imfinfo([imageDirectory filesep fileName]);
    rowSize = imInfo.Height;
    colSize = imInfo.Width;
    
    % Show the frame on the screen in the current figure
    hold on;
    imshow (image, []), title (num2str (firstImage));
    hold off;

    % Make sure the axis are correct
    axis([1 size(image,2) 1 size(image,1)]);

    % Read the MPM file
%     try
%         load ([resultsDirectory filesep 'MPM.mat']);
        
        % Identify the real cells (at least one coord different from zero)
        realCellIndex = find(MPM(:,1) | MPM(:,2));

        % Find the row indices from a transposed MPM matrix
        cellsWithNums = zeros(size(MPM,1),3);
        cellsWithNums(:,3) = [1:1:size(MPM,1)]';

        % Grab all rows in MPM, so that the row indices correspond to the cells
        cellsWithNums(:,1:2) = MPM(:,1:2);

        % Now take the cells identified as real cells (at least one coord different from zero)
        % and plot those as red dots. The cell number is written as colored text on the current axes.
        hold on;
        dots = plot (cellsWithNums (realCellIndex, 1), cellsWithNums (realCellIndex, 2), 'r.');
        set(dots,'Tag','dots','ButtonDownFcn',@ptManTrackCells);
        txt = text (cellsWithNums (realCellIndex, 1), cellsWithNums (realCellIndex, 2), ...
                    num2str (cellsWithNums (realCellIndex, 3)), 'Color', 'r');

        % That's it: wait for the next user action
        hold off;
%     catch
%         msgStr = ['No MPM file present in ' resultsDirectory '. Exiting...'];
%         h = errordlg(msgStr);
%         uiwait(h);
%     end

    % create structure of handles
    handlesNew = guihandles(sliderHandle); 
    
    % add relevant data
    handlesNew.imageDirectory   = imageDirectory;
    handlesNew.resultsDirectory = resultsDirectory;
    handlesNew.imageName        = imageName;
    handlesNew.firstImage       = firstImage;
    handlesNew.lastImage        = lastImage;
    handlesNew.imageNameList    = imageNameList;
    handlesNew.intensityMax     = intensityMax;
    handlesNew.imageRange       = imageRange;
    handlesNew.MPM              = MPM;
    handlesNew.manualMPM        = zeros(1,size(MPM,2));
    handlesNew.validFrames      = validFrames;
    handlesNew.rowSize          = rowSize;
    handlesNew.colSize          = colSize;
    handlesNew.saveDirectory    = saveDirectory;
    handlesNew.bodyName         = bodyName;
    
    % save the structure
    guidata(sliderHandle,handlesNew)     
end