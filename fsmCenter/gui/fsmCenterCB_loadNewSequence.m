function fsmCenterCB_loadNewSequence
% fsmCenterCB_loadNewSequence
%
% This function is a callback for the pushImShow_Callback function of fsmCenter
%
% INPUT   None
%
% OUTPUT  None
%
% REMARK  A variable img is returned into Matlab base workspace
%
% Andre Kerstens, 12/10/2004

% Current directory
oldDir=cd;

% Read current image path from fsmCenter
hFsm=findall(0,'Tag','fsmCenter','Name','fsmCenter');
if ~isempty(hFsm)
    handles=guidata(hFsm);
    imagePath=get(handles.textCurrentImage,'String');
    if ~isempty(imagePath)
        if isdir(imagePath)
            cd(imagePath);
        end
    end
end

% Select image
[fileName, imageDirectory] = uigetfile(...
    {'*.tif;*.tiff;*.jpg;*.jpeg','Image Files (*.tif,*.tiff,*.jpg,*.jpeg)';
    '*.tif','TIF files (*.tif)'
    '*.tiff','TIFF files (*.tiff)'
    '*.jpg;','JPG files (*.jpg)'
    '*.jpeg;','JPEG files (*.jpeg)'
    '*.*','All Files (*.*)'},...
    'Select image to load');

if isa(fileName,'char') && isa(imageDirectory,'char')
    % Clear current window's children
    clf;
    
    [path,bodyName] = getFilenameBody(fileName);
    dirList = dir([imageDirectory, filesep, bodyName, '*']);    
    
    % If images do not exist then exit
    if isempty(dirList)
        h = errordlg(['No images starting with ' bodyName ' can be found in ' imageDirectory '.']);
        uiwait(h);
        close;
        return
    end

    % Rearrange images according to number
    imageNum = zeros(1, length(dirList));
    for i = 1:length(dirList)
        [path, bodyName, no] = getFilenameBody(dirList(i).name);
        imageNum(i) = str2double(no);
    end    
    imageNameList = {dirList(imageNum).name};
    
    % Find out what the first image nr is (not necessarily 1)
    [path,bodyName,no] = getFilenameBody(char(imageNameList(1)));

    firstImage = str2double(no);
            
    % Calculate the image range
    imageRange = length(dirList);

    % Calculate the last image nr
    lastImage = firstImage + imageRange - 1;

    % Generate the slider step values for the uicontrol
    % First the arrow slide step (1):
    sliderStep(1) = 1 / imageRange;

    % Then the through step size (5):
    sliderStep(2) = 5 / imageRange;

    % Draw the frame counter in the figure; it is identified by the tag picturecount
    uicontrol ('Style', 'text',...
        'Units', 'normalized',...
        'Tag', 'pictureCount',...
        'String', num2str(firstImage),...
        'Position', [0.02,0.93,0.06,0.06]);

    % Draw the slider in the figure; it is identified by the tag pictureslide and calls
    % the function ptShowSlidingFrames when moved
    uicontrol ('Style', 'slider', ...
        'Units', 'normalized',...
        'Value', sliderStep(1), ...
        'Min', sliderStep(1), ...
        'Max', 1, ...
        'SliderStep', sliderStep, ...
        'Callback', 'fsmShowImageSequence', ...
        'Tag', 'pictureSlide', ...
        'Position', [0.02,0.02,0.06,0.9]);
    
    % Read the image frame from disk
    image = double(imread([imageDirectory, filesep, char(imageNameList(1))]));
    
    % Get a handle to the figure
    imgHandle = gcf;

    % Show the frame on the screen in the current figure
    hold on;
    imshow(image, []);
    hold off;

    % Make sure the axis are correct
    axis([1 size(image,2) 1 size(image,1)]);
    
    % Add title
    titleStr = [imageDirectory, fileName, ' (', num2str(size(image, 1)), 'x', num2str(size(image, 2)), ')'];
    
    set(imgHandle, 'Name', titleStr, 'NumberTitle', 'off');
    
    % Add univocal tag
    set(imgHandle, 'Tag', 'ViewPanel');
    
    % Return image to MATLAB base workspace
    assignin('base','image',image);
    
    % Store the relevant values in the handles struct
    handles.imageSeq.firstImage = firstImage;
    handles.imageSeq.lastImage = lastImage;
    handles.imageSeq.imageNameList = imageNameList;
    handles.imageSeq.imageDirectory = imageDirectory;
    handles.imageSeq.bodyName = bodyName;
    handles.imageSeq.imageRange = imageRange;

    % Write back the handles struct
    guidata(hFsm, handles);
else
    close(gcf);
    return
end

% Add common menus
fsmCenterCB_addSeqToolsMenu(gcf);

% Create menu with additional tools
figure(gcf);
hMenu=uimenu('Label','More Tools');
uimenu(hMenu,'Label','Find edges','Callback','fsmCenterCB_findSeqEdges;','Accelerator','E','Separator','On');

% Go back to old directory
cd(oldDir);


