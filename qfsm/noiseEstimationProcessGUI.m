function varargout = noiseEstimationProcessGUI(varargin)
% noiseEstimationProcessGUI M-file for noiseEstimationProcessGUI.fig
%      noiseEstimationProcessGUI, by itself, creates a new noiseEstimationProcessGUI or raises the existing
%      singleton*.
%
%      H = noiseEstimationProcessGUI returns the handle to a new noiseEstimationProcessGUI or the handle to
%      the existing singleton*.
%
%      noiseEstimationProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in noiseEstimationProcessGUI.M with the given input arguments.
%
%      noiseEstimationProcessGUI('Property','Value',...) creates a new noiseEstimationProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before noiseEstimationProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to noiseEstimationProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help noiseEstimationProcessGUI

% Last Modified by GUIDE v2.5 17-Jun-2011 11:12:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @noiseEstimationProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @noiseEstimationProcessGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before noiseEstimationProcessGUI is made visible.
function noiseEstimationProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)


processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},...
    'initChannel',1);

% Choose default command line output for noiseEstimationProcessGUI
handles.output = hObject;

% ---------------------- Channel Setup -------------------------
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

% Save the image directories and names (for cropping preview)
userData.nFrames = userData.MD.nFrames_;
userData.imRectHandle.isvalid=0;
userData.firstImage = funParams.firstImage;
userData.lastImage = funParams.lastImage;
userData.cropROI = funParams.cropROI;
userData.filterSigma = funParams.filterSigma;
userData.previewFig=-1;

% Read the first image and update the sliders max value and steps
props = get(handles.listbox_selectedChannels, {'UserData','Value'});
userData.chanIndx = props{1}(props{2});
firstImage = userData.firstImage(userData.chanIndx);
lastImage = userData.lastImage(userData.chanIndx);
set(handles.edit_firstImage,'String',firstImage);
set(handles.edit_lastImage,'String',lastImage);
set(handles.edit_frameNumber,'String',firstImage);
set(handles.slider_frameNumber,'Min',firstImage,'Value',firstImage,'Max',lastImage,...
    'SliderStep',[1/double(lastImage-firstImage)  10/double(lastImage-firstImage)]);
userData.imIndx=firstImage;
userData.imData=userData.MD.channels_(userData.chanIndx).loadImage(userData.imIndx);
set(handles.edit_filterSigma,'String',userData.filterSigma(userData.chanIndx));

% ----------------------------------------------------------------

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = noiseEstimationProcessGUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(~, ~, handles)
% Delete figure
delete(handles.figure1);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, ~, handles)
% Notify the package GUI that the setting panel is closed
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'helpFig') && ishandle(userData.helpFig)
   delete(userData.helpFig) 
end

if ishandle(userData.previewFig), delete(userData.previewFig); end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);

% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

function close_previewFig(hObject, eventdata)
handles = guidata(get(hObject,'UserData'));
set(handles.checkbox_crop,'Value',0);
update_data(handles.checkbox_crop, eventdata, handles);

 % --- Executes on button press in checkbox_crop.
function update_data(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');

% Retrieve the channel index
props=get(handles.listbox_selectedChannels,{'UserData','Value'});
chanIndx = props{1}(props{2});
imIndx = get(handles.slider_frameNumber,'Value');

% If channel index has been modified, load new frame values
if (chanIndx~=userData.chanIndx)
    set(handles.edit_filterSigma,'String',userData.filterSigma(userData.chanIndx));
    firstImage = userData.firstImage(chanIndx);
    lastImage  = userData.lastImage(chanIndx);
    set(handles.edit_firstImage,'String',firstImage);
    set(handles.edit_lastImage,'String',lastImage);
    set(handles.edit_frameNumber,'String',firstImage);
    set(handles.slider_frameNumber,'Min',firstImage,'Value',firstImage,'Max',lastImage,...
    'SliderStep',[1/double(lastImage-firstImage)  10/double(lastImage-firstImage)]);
end

% Load a new image if either the image number or the channel has been changed
if (chanIndx~=userData.chanIndx) ||  (imIndx~=userData.imIndx)
    userData.imData=userData.MD.channels_(chanIndx).loadImage(imIndx);
    userData.updateImage=1;
    userData.chanIndx=chanIndx;
    userData.imIndx=imIndx;
    
    % Update ROI
    if userData.imRectHandle.isvalid
        userData.cropROI=getPosition(userData.imRectHandle);
    end
    
else
    userData.updateImage=0;
end

% In case of crop previewing mode
if get(handles.checkbox_crop,'Value')
    % Create figure if non-existing or closed
    if ~isfield(userData, 'previewFig') || ~ishandle(userData.previewFig)
        userData.previewFig = figure('Name','Select the background region to crop',...
            'DeleteFcn',@close_previewFig,'UserData',handles.figure1);
        axes('Position',[.05 .05 .9 .9]);
        userData.newFigure = 1;
    else
        figure(userData.previewFig);
        userData.newFigure = 0;
    end
    
    % Retrieve the image object handle
    imHandle = findobj(userData.previewFig,'Type','image');
    if userData.newFigure || userData.updateImage 
        if isempty(imHandle)
            imHandle=imshow(mat2gray(userData.imData));
            axis off;
        else
            set(imHandle,'CData',mat2gray(userData.imData));
        end
    end
        

   if userData.imRectHandle.isvalid
        % Update the imrect position
        setPosition(userData.imRectHandle,userData.cropROI)
   else 
        % Create a new imrect object and store the handle
        userData.imRectHandle = imrect(get(imHandle,'Parent'),userData.cropROI);
        fcn = makeConstrainToRectFcn('imrect',get(imHandle,'XData'),get(imHandle,'YData'));
        setPositionConstraintFcn(userData.imRectHandle,fcn);
    end
else
    if userData.imRectHandle.isvalid, 
        userData.cropROI=getPosition(userData.imRectHandle);
    end
    
    % Close the figure if applicable
    if ishandle(userData.previewFig), delete(userData.previewFig); end
    
end
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);

% --- Executes on slider movement.
function frameNumberEdition_Callback(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');

% Retrieve the value of the selected image
if strcmp(get(hObject,'Tag'),'edit_frameNumber')
    frameNumber = str2double(get(handles.edit_frameNumber, 'String'));
else
    frameNumber = get(handles.slider_frameNumber, 'Value');
end
frameNumber=round(frameNumber);

firstImage = str2double(get(handles.edit_firstImage,'String'));
lastImage = str2double(get(handles.edit_lastImage,'String'));

% Check the validity of the frame values
if isnan(firstImage) || isnan(lastImage) || isnan(frameNumber)
    warndlg('Please provide a valid frame value.','Setting Error','modal');
end
if firstImage>lastImage, firstImage=lastImage; end

firstImage = max(firstImage,1);
lastImage = min(lastImage,userData.nFrames);
frameNumber = min(max(frameNumber,firstImage),lastImage);

% Store value
userData.firstImage(userData.chanIndx) = firstImage;
userData.lastImage(userData.chanIndx) = lastImage;

set(handles.slider_frameNumber,'Min',firstImage,'Value',frameNumber,'Max',lastImage,...
    'SliderStep',[1/double(lastImage-firstImage)  10/double(lastImage-firstImage)]);
set(handles.edit_frameNumber,'String',frameNumber);
set(handles.edit_firstImage,'String',firstImage);
set(handles.edit_lastImage,'String',lastImage);

% Save data and update graphics
set(handles.figure1, 'UserData', userData);
guidata(hObject, handles);
update_data(hObject,eventdata,handles);


function edit_filterSigma_Callback(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');
% Retrieve the channel index
props=get(handles.listbox_selectedChannels,{'UserData','Value'});
chanIndx = props{1}(props{2});

value = str2double(get(hObject, 'String'));
if ~(value>0), 
    % Reset old value
    set(hObject,'String',userData.filterSigma(chanIndx))
else
    % Update the sigma value in the stored array-
    userData.filterSigma(chanIndx) = value;
    set(handles.edit_filterSigma, 'String', num2str(value));
    guidata(hObject, handles);
end
% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

% Input check
if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end

channelIndex = get (handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;

% Process Sanity check ( only check underlying data )
userData = get(handles.figure1, 'UserData');

try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

% Retrieve GUI-defined parameters
if userData.imRectHandle.isvalid
    userData.cropROI=getPosition(userData.imRectHandle);
end
funParams.cropROI=userData.cropROI;
funParams.firstImage=userData.firstImage;
funParams.lastImage=userData.lastImage;

% Save the filterSigma if different from psfSigma
% In order not to override filterSigma in batch movie set up
if ~isequal(userData.filterSigma,[userData.MD.channels_.psfSigma_])
    % funParams.filterSigma = userData.filterSigma;
    funParams.filterSigma = str2double(get(handles.edit_filterSigma,'String'));
end

% Set parameters and update main window
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);
