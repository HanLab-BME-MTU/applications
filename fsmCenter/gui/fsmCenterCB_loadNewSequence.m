function fsmCentercb_loadNewSequence
% fsmCentercb_loadNewSequence
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
[fileName,dirName] = uigetfile(...
    {'*.tif;*.tiff;*.jpg;*.jpeg','Image Files (*.tif,*.tiff,*.jpg,*.jpeg)';
    '*.tif','TIF files (*.tif)'
    '*.tiff','TIFF files (*.tiff)'
    '*.jpg;','JPG files (*.jpg)'
    '*.jpeg;','JPEG files (*.jpeg)'
    '*.*','All Files (*.*)'},...
    'Select image to load');
if(isa(fileName,'char') & isa(dirName,'char'))
    
    % Clear current window's children
    clf;

    % Fetch the jobvalues
    imageDirectory = dirName;
    imageName      = fileName;

    % Find out what part of the filename describes the images and which part
    % is just counting them through.    
    number = 0;
    countNum = 0;
    while ~isnan(number) & (countNum < 3)
       countNum = countNum + 1;
       number = str2num(fileName(end-(4+countNum):end-4));
    end

    % Extract the body of the filename and store in handles struct
    bodyName = fileName(1:(end-(4+countNum)));
    
    % Create a list of files present in the image directory selected by the user
    dirList = dir(imageDirectory);
    dirList = struct2cell(dirList);
    dirList = dirList(1,:);
    
    % Find all files within this directory with the same name as the selected filename
    ind = strmatch(bodyName, dirList);
    dirList = dirList(ind)';
    firstImage = 1;

    % If images do not exist then exit
    if isempty(dirList)
        h=errordlg(['No images starting with ' bodyName ' can be found in ' imageDirectory '.']);
        uiwait(h);
        close;
        return
    end
    
    % Rearrange images according to number
    for jRearange = 1:length(dirList)
       tmpName = char(dirList(jRearange));
       try
          imageNum(jRearange) = str2num(tmpName(length(bodyName)+1:end-4));
       catch
          dirList(jRearange) = [];
       end
    end  % for jRearange
    [dummy,imageNumInd] = sort(imageNum);
    imageNameList = dirList(imageNumInd);
    lastImage = length(dirList);
    
    % Calculate the image range
    imageRange = lastImage;

    % Generate the slider step values for the uicontrol
    % First the arrow slide step (1):
    slider_step(1) = 1 / imageRange;

    % Then the through step size (5):
    slider_step(2) = 5 / imageRange;

    % Draw the frame counter in the figure; it is identified by the tag picturecount
    imageCounterHandle = uicontrol ('Style', 'text',...
                                    'Units', 'normalized',...
                                    'Tag', 'pictureCount',...
                                    'Position', [0.02,0.93,0.05,0.06]);

    % Set the frame counter to the first image number
    set (imageCounterHandle, 'String', num2str(firstImage));

    % Draw the slider in the figure; it is identified by the tag pictureslide and calls
    % the function ptShowSlidingFrames when moved
    sliderHandle = uicontrol ('Style', 'slider', ...
                              'Units', 'normalized', ... 
                              'Value', 1/(imageRange), ...
                              'Min', 1/(imageRange), ...
                              'Max', 1, ...
                              'SliderStep', slider_step, ...
                              'Callback', 'fsmShowImageSequence', ...
                              'Tag', 'pictureSlide', ...
                              'Position', [0.02,0.02,0.05,0.9]);

    % Read the image frame from disk
    %cd (imageDirectory);
    fileName = char(imageNameList(firstImage));
    image = double(imread([imageDirectory,fileName]));
    
    % Get a handle to the figure
    imgHandle = gcf;

    % Show the frame on the screen in the current figure
    hold on;
    imshow(image, []);
    hold off;

    % Make sure the axis are correct
    axis([1 size(image,2) 1 size(image,1)]);
    
    % Add title
    titleStr=[dirName,fileName,' (',num2str(size(image,1)),'x',num2str(size(image,2)),')'];
    set(imgHandle,'Name',titleStr,'NumberTitle','off');
    
    % Add univocal tag
    set(imgHandle,'Tag','ViewPanel');
    
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
uimenu(hMenu,'Label','Find edges','Callback','fsmCenterCB_findEdges;','Accelerator','E','Separator','On');

% Go back to old directory
cd(oldDir);


