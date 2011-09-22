function varargout = stageDriftCorrectionProcessGUI(varargin)
% stageDriftCorrectionProcessGUI M-file for stageDriftCorrectionProcessGUI.fig
%      stageDriftCorrectionProcessGUI, by itself, creates a new stageDriftCorrectionProcessGUI or raises the existing
%      singleton*.
%
%      H = stageDriftCorrectionProcessGUI returns the handle to a new stageDriftCorrectionProcessGUI or the handle to
%      the existing singleton*.
%
%      stageDriftCorrectionProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in stageDriftCorrectionProcessGUI.M with the given input arguments.
%
%      stageDriftCorrectionProcessGUI('Property','Value',...) creates a new stageDriftCorrectionProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stageDriftCorrectionProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stageDriftCorrectionProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stageDriftCorrectionProcessGUI

% Last Modified by GUIDE v2.5 08-Sep-2011 13:16:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stageDriftCorrectionProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @stageDriftCorrectionProcessGUI_OutputFcn, ...
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


% --- Executes just before stageDriftCorrectionProcessGUI is made visible.
function stageDriftCorrectionProcessGUI_OpeningFcn(hObject,eventdata,handles,varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},...
    'initChannel',1);

% Set process parameters
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

set(handles.edit_referenceFramePath,'String',funParams.referenceFramePath);

userData.numParams = {'I0','sDN','GaussRatio','alpha','minCorLength','maxFlowSpeed'};
cellfun(@(x) set(handles.(['edit_' x]),'String',funParams.(x)),userData.numParams);
set(handles.checkbox_doPreReg,'Value',funParams.doPreReg);

% Save the image directories and names (for cropping preview)
userData_main = get(userData.mainFig, 'UserData');
userData.imageFileNames = userData_main.MD(userData_main.id).getImageFileNames();
userData.imDirs  = userData_main.MD(userData_main.id).getChannelPaths();
userData.nFrames = userData_main.MD(userData_main.id).nFrames_;
userData.imRectHandle.isvalid=0;
userData.cropROI = funParams.cropROI;
userData.previewFig=-1;

% Read the first image and update the sliders max value and steps
props = get(handles.listbox_selectedChannels, {'String','Value'});
userData.chanIndx = find(strcmp(props{1}{props{2}},userData.imDirs));
set(handles.edit_frameNumber,'String',1);
set(handles.slider_frameNumber,'Min',1,'Value',1,'Max',userData.nFrames,...
    'SliderStep',[1/double(userData.nFrames-1)  10/double(userData.nFrames-1)]);
userData.imIndx=1;
userData.imData=mat2gray(imread([userData.imDirs{userData.chanIndx} filesep...
        userData.imageFileNames{userData.chanIndx}{userData.imIndx}]));
    
set(handles.listbox_selectedChannels,'Callback',@(h,event) update_data(h,event,guidata(h)));
    
% Choose default command line output for stageDriftCorrectionProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = stageDriftCorrectionProcessGUI_OutputFcn(~, ~, handles) 
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

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

% Check user input
userData = get(handles.figure1, 'UserData');
if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
else
    channelIndex = get(handles.listbox_selectedChannels, 'Userdata');
    funParams.ChannelIndex = channelIndex;
end

% Retrieve reference frame path
funParams.referenceFramePath=get(handles.edit_referenceFramePath,'String');
if isempty(funParams.referenceFramePath)
    errordlg('Please select a reference frame.','Setting Error','modal')
    return;
end

% Read numeric information
for i = 1:numel(userData.numParams),
    value = get(handles.(['edit_' userData.numParams{i}]),'String');
    if isempty(value)
        errordlg('Please enter a valid value.','Setting Error','modal')
        return;
    end
    funParams.(userData.numParams{i})=str2double(value); 
end

% Read cropRoi if window if
if userData.imRectHandle.isvalid
    userData.cropROI=getPosition(userData.imRectHandle);
end
funParams.cropROI = userData.cropROI;
funParams.doPreReg = get(handles.checkbox_doPreReg,'Value');

% Process Sanity check ( only check underlying data )
try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

% Set parameters
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);

% --- Executes on button press in pushbutton_selectReferenceFrame.
function pushbutton_selectReferenceFrame_Callback(hObject, eventdata, handles)

[file path]=uigetfile({'*.tif;*.TIF;*.stk;*.STK;*.bmp;*.BMP;*.jpg;*.JPG',...
    'Image files (*.tif,*.stk,*.bmp,*.jpg)'},...
    'Select the reference frame');
if ~isequal(file,0) && ~isequal(path,0)
    set(handles.edit_referenceFramePath,'String',[path file]);
end

 % --- Executes on button press in checkbox_crop.
function update_data(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');

% Retrieve the channel index
props=get(handles.listbox_selectedChannels,{'String','Value'});
chanIndx = find(strcmp(props{1}{props{2}},userData.imDirs));
imIndx = get(handles.slider_frameNumber,'Value');

% Load a new image if either the image number or the channel has been changed
if (chanIndx~=userData.chanIndx) ||  (imIndx~=userData.imIndx)
    % Update image flag and dat
    userData.imData=mat2gray(imread([userData.imDirs{chanIndx} filesep...
        userData.imageFileNames{chanIndx}{imIndx}]));
    userData.updateImage=1;
    userData.chanIndx=chanIndx;
    userData.imIndx=imIndx;
        
    % Update roi
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
        userData.previewFig = figure('Name','Select the region to crop',...
            'DeleteFcn',@close_previewFig,'UserData',handles.figure1);
        userData.newFigure = 1;
    else
        figure(userData.previewFig);
        userData.newFigure = 0;
    end
    
    % Retrieve the image object handle
    imHandle =findobj(userData.previewFig,'Type','image');
    if userData.newFigure || userData.updateImage
        if isempty(imHandle)
            imHandle=imshow(userData.imData);
            axis off;
        else
            set(imHandle,'CData',userData.imData);
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
    % Save the roi if applicable
    if userData.imRectHandle.isvalid, 
        userData.cropROI=getPosition(userData.imRectHandle); 
    end
    % Close the figure if applicable
    if ishandle(userData.previewFig), delete(userData.previewFig); end
end
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);

function close_previewFig(hObject, eventdata)
handles = guidata(get(hObject,'UserData'));
set(handles.checkbox_crop,'Value',0);
update_data(handles.checkbox_crop, eventdata, handles);


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

% Check the validity of the frame values
if isnan(frameNumber)
    warndlg('Please provide a valid frame value.','Setting Error','modal');
end
frameNumber = min(max(frameNumber,1),userData.nFrames);

% Store value
set(handles.slider_frameNumber,'Value',frameNumber);
set(handles.edit_frameNumber,'String',frameNumber);

% Save data and update graphics
set(handles.figure1, 'UserData', userData);
guidata(hObject, handles);
update_data(hObject,eventdata,handles);


% --- Executes on button press in pushbutton_loadNoiseParameters.
function pushbutton_loadNoiseParameters_Callback(hObject, eventdata, handles)
[file path]=uigetfile({'*.mat;*.MAT',...
    'Mat files (*.mat)'},...
    'Select the file containing the noise model parameters');
if ~isequal(file,0) && ~isequal(path,0)
    noiseParams={'I0','sDN','GaussRatio'};
    vars=whos(noiseParams{:},'-file',[path file]);
    if numel(vars)~=numel(noiseParams),
        errordlg('Please select a file containing valid noise model parameters');
        return 
    end
    
    s=load([path file],noiseParams{:});
    for i=1:numel(noiseParams)
        set(handles.(['edit_' noiseParams{i}]),'String',s.(noiseParams{i}),...
            'Enable','on');
        
    end
end
