function varargout = steerableFilteringProcessGUI(varargin)
% steerableFilteringProcessGUI M-file for steerableFilteringProcessGUI.fig
%      steerableFilteringProcessGUI, by itself, creates a new steerableFilteringProcessGUI or raises the existing
%      singleton*.
%
%      H = steerableFilteringProcessGUI returns the handle to a new steerableFilteringProcessGUI or the handle to
%      the existing singleton*.
%
%      steerableFilteringProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in steerableFilteringProcessGUI.M with the given input arguments.
%
%      steerableFilteringProcessGUI('Property','Value',...) creates a new steerableFilteringProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before steerableFilteringProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to steerableFilteringProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help steerableFilteringProcessGUI

% Last Modified by GUIDE v2.5 07-Sep-2012 12:49:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @steerableFilteringProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @steerableFilteringProcessGUI_OutputFcn, ...
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


% --- Executes just before steerableFilteringProcessGUI is made visible.
function steerableFilteringProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',0);

% ---------------------- Channel Setup -------------------------
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

% Set up available input channels
set(handles.listbox_availableChannels,'String',userData.MD.getChannelPaths(), ...
    'UserData',1:numel(userData.MD.channels_));

channelIndex = funParams.ChannelIndex;

% Find any parent process
userData.parentProc = userData.crtPackage.getParent(userData.procID);
if isempty(userData.crtPackage.processes_{userData.procID}) && ~isempty(userData.parentProc)
    % Check existence of all parent processes
    emptyParentProc = any(cellfun(@isempty,userData.crtPackage.processes_(userData.parentProc)));
    if ~emptyParentProc
        % Intersect channel index with channel index of parent processes
        parentChannelIndex = @(x) userData.crtPackage.processes_{x}.funParams_.ChannelIndex;
        for i = userData.parentProc
            channelIndex = intersect(channelIndex,parentChannelIndex(i));
        end
    end
   
end

if ~isempty(channelIndex)
    channelString = userData.MD.getChannelPaths(channelIndex);
else
    channelString = {};
end

set(handles.listbox_selectedChannels,'String',channelString,...
    'UserData',channelIndex);

set(handles.popupmenu_imageflattenedflag,'String', {'Original Images','Flattened Images'});
set(handles.popupmenu_imageflattenedflag,'Value', funParams.ImageFlattenFlag);


set(handles.edit_subsample_number,'String',funParams.Sub_Sample_Num);

set(handles.edit_BaseSteerableFilterSigma,'String',funParams.BaseSteerableFilterSigma);
set(handles.edit_levelsofsteerablefilters,'String',funParams.Levelsofsteerablefilters);

% Update user data and GUI data
handles.output = hObject;
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = steerableFilteringProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% Delete figure
delete(handles.figure1);


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% Call back function of 'Apply' button
userData = get(handles.figure1, 'UserData');

% -------- Check user input --------
if isempty(get(handles.listbox_selectedChannels, 'String'))
   errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal') 
    return;
end
channelIndex = get (handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;

BaseSteerableFilterSigma = str2double(get(handles.edit_BaseSteerableFilterSigma, 'String'));
if isnan(BaseSteerableFilterSigma) || BaseSteerableFilterSigma < 0
    errordlg(['Please provide a valid input for '''...

    get(handles.text_BaseSteerableFilterSigma,'String') '''.'],'Setting Error','modal');
    return;
end
funParams.BaseSteerableFilterSigma=BaseSteerableFilterSigma;

Levelsofsteerablefilters = str2double(get(handles.edit_levelsofsteerablefilters, 'String'));
if isnan(Levelsofsteerablefilters) || Levelsofsteerablefilters < 0
    errordlg(['Please provide a valid input for '''...

    get(handles.text_Levelsofsteerablefilters,'String') '''.'],'Setting Error','modal');
    return;
end

funParams.Levelsofsteerablefilters=Levelsofsteerablefilters;

Sub_Sample_Num  = str2double(get(handles.edit_subsample_number, 'String'));
if isnan(Sub_Sample_Num) || Sub_Sample_Num < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_subsample,'String') '''.'],'Setting Error','modal');
    return;
end
funParams.Sub_Sample_Num  = Sub_Sample_Num;

funParams.ImageFlattenFlag = get(handles.popupmenu_imageflattenedflag,'Value');

% -------- Process Sanity check --------
% ( only check underlying data )

try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

processGUI_ApplyFcn(hObject, eventdata, handles,funParams);

% --- Executes on button press in checkbox_all.
function checkbox_all_Callback(hObject, eventdata, handles)

% Hint: get(hObject,'Value') returns toggle state of checkbox_all
contents1 = get(handles.listbox_availableChannels, 'String');

chanIndex1 = get(handles.listbox_availableChannels, 'Userdata');
chanIndex2 = get(handles.listbox_selectedChannels, 'Userdata');

% Return if listbox1 is empty
if isempty(contents1)
    return;
end

switch get(hObject,'Value')
    case 1
        set(handles.listbox_selectedChannels, 'String', contents1);
        chanIndex2 = chanIndex1;
    case 0
        set(handles.listbox_selectedChannels, 'String', {}, 'Value',1);
        chanIndex2 = [ ];
end
set(handles.listbox_selectedChannels, 'UserData', chanIndex2);


% --- Executes on button press in pushbutton_select.
function pushbutton_select_Callback(hObject, eventdata, handles)
% call back function of 'select' button

contents1 = get(handles.listbox_availableChannels, 'String');
contents2 = get(handles.listbox_selectedChannels, 'String');
id = get(handles.listbox_availableChannels, 'Value');

% If channel has already been added, return;
chanIndex1 = get(handles.listbox_availableChannels, 'Userdata');
chanIndex2 = get(handles.listbox_selectedChannels, 'Userdata');

for i = id
    if any(strcmp(contents1{i}, contents2) )
        continue;
    else
        contents2{end+1} = contents1{i};
        chanIndex2 = cat(2, chanIndex2, chanIndex1(i));

    end
end

set(handles.listbox_selectedChannels, 'String', contents2, 'Userdata', chanIndex2);


% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)
% Call back function of 'delete' button
contents = get(handles.listbox_selectedChannels,'String');
id = get(handles.listbox_selectedChannels,'Value');

% Return if list is empty
if isempty(contents) || isempty(id)
    return;
end

% Delete selected item
contents(id) = [ ];

% Delete userdata
chanIndex2 = get(handles.listbox_selectedChannels, 'Userdata');
chanIndex2(id) = [ ];
set(handles.listbox_selectedChannels, 'Userdata', chanIndex2);

% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if (id >length(contents) && id>1)
    set(handles.listbox_selectedChannels,'Value',length(contents));
end
% Refresh listbox
set(handles.listbox_selectedChannels,'String',contents);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% Notify the package GUI that the setting panel is closed
userData = get(handles.figure1, 'UserData');

if ishandle(userData.helpFig), delete(userData.helpFig); end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end


function edit_Levelsofsteerablefilters_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Levelsofsteerablefilters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Levelsofsteerablefilters as text
%        str2double(get(hObject,'String')) returns contents of edit_Levelsofsteerablefilters as a double


% --- Executes during object creation, after setting all properties.
function edit_Levelsofsteerablefilters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Levelsofsteerablefilters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_imageflattenedflag.
function popupmenu_imageflattenedflag_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_imageflattenedflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_imageflattenedflag contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_imageflattenedflag


% --- Executes during object creation, after setting all properties.
function popupmenu_imageflattenedflag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_imageflattenedflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',{'Original Images','Flattened Images'});
set(hObject,'Value',2);


function edit_subsample_number_Callback(hObject, eventdata, handles)
% hObject    handle to edit_subsample_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_subsample_number as text
%        str2double(get(hObject,'String')) returns contents of edit_subsample_number as a double


% --- Executes during object creation, after setting all properties.
function edit_subsample_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_subsample_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
