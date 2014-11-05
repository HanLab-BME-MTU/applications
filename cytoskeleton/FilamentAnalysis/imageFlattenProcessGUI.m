function varargout = imageFlattenProcessGUI(varargin)
% imageFlattenProcessGUI M-file for imageFlattenProcessGUI.fig
%      imageFlattenProcessGUI, by itself, creates a new imageFlattenProcessGUI or raises the existing
%      singleton*.
%
%      H = imageFlattenProcessGUI returns the handle to a new imageFlattenProcessGUI or the handle to
%      the existing singleton*.
%
%      imageFlattenProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in imageFlattenProcessGUI.M with the given input arguments.
%
%      imageFlattenProcessGUI('Property','Value',...) creates a new imageFlattenProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imageFlattenProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imageFlattenProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imageFlattenProcessGUI

% Last Modified by GUIDE v2.5 07-Sep-2012 11:20:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imageFlattenProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @imageFlattenProcessGUI_OutputFcn, ...
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


% --- Executes just before imageFlattenProcessGUI is made visible.
function imageFlattenProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

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

set(handles.edit_GaussFilterSigma,'String',funParams.GaussFilterSigma);
set(handles.edit_subsample_number,'String',funParams.Sub_Sample_Num);

set(handles.popupmenu_flatten_method,'String',{'Log','Sqrt','Power 2/3','Nothing'});
set(handles.popupmenu_flatten_method,'Value',funParams.method_ind);

set(handles.edit_timefilter_sigma,'String',funParams.TimeFilterSigma);


% flattenMethods = userData.crtProc.getMethods();
% set(handles.popupmenu_flatteningMethods,'String',{flattenMethods(:).name},...
%     'Value',funParams.MethodIndx);

% Update user data and GUI data
handles.output = hObject;
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = imageFlattenProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% Notify the package GUI that the setting panel is closed
userData = get(handles.figure1, 'UserData');

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


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
% update_data(hObject,eventdata,handles);

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
% update_data(hObject,eventdata,handles);


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call back function of 'Apply' button
userData = get(handles.figure1, 'UserData');

% -------- Check user input --------
if isempty(get(handles.listbox_selectedChannels, 'String'))
   errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal') 
    return;
end
channelIndex = get (handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;
 
gaussFilterSigma = str2double(get(handles.edit_GaussFilterSigma, 'String'));
if isnan(gaussFilterSigma) || gaussFilterSigma < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_GaussFilterSigma,'String') '''.'],'Setting Error','modal');
    return;
end
funParams.GaussFilterSigma=gaussFilterSigma;

funParams.method_ind = get(handles.popupmenu_flatten_method,'Value');

TimeFilterSigma = str2double(get(handles.edit_timefilter_sigma, 'String'));
if isnan(TimeFilterSigma) || TimeFilterSigma < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_timefilter,'String') '''.'],'Setting Error','modal');
    return;
end
funParams.TimeFilterSigma=TimeFilterSigma;

Sub_Sample_Num  = str2double(get(handles.edit_subsample_number, 'String'));
if isnan(Sub_Sample_Num) || Sub_Sample_Num < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_subsample,'String') '''.'],'Setting Error','modal');
    return;
end

funParams.Sub_Sample_Num  = Sub_Sample_Num;

          



% -------- Process Sanity check --------
% ( only check underlying data )

try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

% Set parameters and update main window
processGUI_ApplyFcn(hObject, eventdata, handles, funParams);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% Delete figure
delete(handles.figure1);


% --- Executes on selection change in popupmenu_flatten_method.
function popupmenu_flatten_method_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_flatten_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_flatten_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_flatten_method


% --- Executes during object creation, after setting all properties.
function popupmenu_flatten_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_flatten_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

 set(hObject,'String',{'Log','Sqrt','Power 2/3','Nothing'});
 set(hObject,'Value',3);
 


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

% Refresh listbox
set(handles.listbox_selectedChannels,'String',contents);


% --- Executes on button press in checkbox_timefilter.
function checkbox_timefilter_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_timefilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_timefilter



function edit_timefilter_sigma_Callback(hObject, eventdata, handles)
% hObject    handle to edit_timefilter_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_timefilter_sigma as text
%        str2double(get(hObject,'String')) returns contents of edit_timefilter_sigma as a double


% --- Executes during object creation, after setting all properties.
function edit_timefilter_sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_timefilter_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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
