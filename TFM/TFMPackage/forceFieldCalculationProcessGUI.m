function varargout = forceFieldCalculationProcessGUI(varargin)
% forceFieldCalculationProcessGUI M-file for forceFieldCalculationProcessGUI.fig
%      forceFieldCalculationProcessGUI, by itself, creates a new forceFieldCalculationProcessGUI or raises the existing
%      singleton*.
%
%      H = forceFieldCalculationProcessGUI returns the handle to a new forceFieldCalculationProcessGUI or the handle to
%      the existing singleton*.
%
%      forceFieldCalculationProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in forceFieldCalculationProcessGUI.M with the given input arguments.
%
%      forceFieldCalculationProcessGUI('Property','Value',...) creates a new forceFieldCalculationProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before forceFieldCalculationProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to forceFieldCalculationProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help forceFieldCalculationProcessGUI

% Last Modified by GUIDE v2.5 05-Apr-2021 00:28:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @forceFieldCalculationProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @forceFieldCalculationProcessGUI_OutputFcn, ...
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


% --- Executes just before forceFieldCalculationProcessGUI is made visible.
function forceFieldCalculationProcessGUI_OpeningFcn(hObject,eventdata,handles,varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:});

% Set process parameters
userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end
funParams = userData.crtProc.funParams_;
userData.numParams ={'PoissonRatio','meshPtsFwdSol','regParam','LcurveFactor'};
cellfun(@(x) set(handles.(['edit_' x]),'String',funParams.(x)),...
    userData.numParams)
set(handles.edit_YoungModulus,'String',funParams.YoungModulus/1000);
set(handles.edit_thickness,'String',funParams.thickness/1000);

% Create pop-up menu for force reconstruction method
solMethodBEMString ={'QR';'svd';'gsvd';'backslash';'1NormReg';'LaplacianReg';'1NormRegLaplacian'};
solMethodBEMData ={'QR';'svd';'gsvd';'backslash';'1NormReg';'LaplacianReg';'1NormRegLaplacian'};
solMethodBEMValue = find(strcmp(funParams.solMethodBEM,solMethodBEMData));
set(handles.popupmenu_solMethodBEM,'String',solMethodBEMString,...
    'UserData',solMethodBEMData,'Value',solMethodBEMValue);

% Create pop-up menu for force reconstruction method
methodString ={'FastBEM';'FTTC';'FEM'};
methodData ={'FastBEM';'FTTC';'FEM'};
methodValue = find(strcmp(funParams.method,methodData));
set(handles.popupmenu_method,'String',methodString,...
    'UserData',methodData,'Value',methodValue);
set(handles.edit_regParam, 'Value', funParams.regParam);

set(handles.useLcurve, 'Value', funParams.useLcurve);
set(handles.checkbox_everyframe, 'Value', funParams.useLcurveEveryFrame);
set(handles.lcorner, 'Value', strcmp(funParams.lcornerOptimal,'lcorner'));
set(handles.optimal, 'Value', strcmp(funParams.lcornerOptimal,'optimal'));
if funParams.useLcurve
    set(handles.checkbox_everyframe, 'Enable', 'on');
    set(handles.groupCornerOptimal.Children,  'Enable', 'on');
else
    set(handles.checkbox_everyframe, 'Enable', 'off');
    set(handles.groupCornerOptimal.Children,  'Enable', 'off');
end

set(handles.checkbox_lastToFirst, 'Value', funParams.lastToFirst);

try
    set(handles.edit_tolx,'String',funParams.tolx);
catch
    disp('tolx has not been set up.')
end

% Update BEM parameter panel
popupmenu_method_Callback(hObject,eventdata,handles);

% Set basis class lookup table path
set(handles.edit_basisClassTblPath,'String',funParams.basisClassTblPath);

% Choose default command line output for forceFieldCalculationProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = forceFieldCalculationProcessGUI_OutputFcn(~, ~, handles) 
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
if isempty(userData), userData = struct(); end

% Read numerical parameters
if isfield(userData,'numParams')
    for i=1:numel(userData.numParams)  
        value = get(handles.(['edit_' userData.numParams{i}]),'String');
        if isempty(value)
            errordlg(['Please enter a valid value for '...
                get(handles.(['text_' userData.numParams{i}]),'String') '.'],...
                'Setting Error','modal')
            return;
        end
        funParams.(userData.numParams{i})=str2double(value); 
    end
end

if isempty(get(handles.edit_YoungModulus,'String'))
    errordlg(['Please enter a valid value for ' ...
        get(handles.text_YoungModulus,'String') '.'],'Setting Error','modal');
    return;
end
funParams.YoungModulus = str2double(get(handles.edit_YoungModulus,'String'))*1000;

thickness = str2double(get(handles.edit_thickness,'String'));
if isnan(thickness) || thickness <=0
    errordlg(['Please enter a valid value for ' ...
        get(handles.text_thickness,'String') '.'],'Setting Error','modal');
    return;
end
funParams.thickness = thickness*1000;

% Read reconstruction method
props=get(handles.popupmenu_method,{'UserData','Value'});
funParams.method=props{1}{props{2}};


% Read BEM solution method
props=get(handles.popupmenu_solMethodBEM,{'UserData','Value'});
try
    funParams.solMethodBEM=props{1}{props{2}};
catch
    funParams.solMethodBEM=props{1}{1};
end
funParams.useLcurve = get(handles.useLcurve, 'Value');
funParams.useLcurveEveryFrame = get(handles.checkbox_everyframe, 'Value');
funParams.lastToFirst = get(handles.checkbox_lastToFirst, 'Value');
funParams.tolx = str2double(get(handles.edit_tolx,'String'));

% Read basis class lookup table path
funParams.basisClassTblPath=get(handles.edit_basisClassTblPath,'String');

% Read radio button about regularization parameter selection
funParams.lcornerOptimal = get(get(handles.groupCornerOptimal,'SelectedObject'),'Tag');

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


% --- Executes on button press in checkbox_filter.
function checkbox_filter_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    enableState='on';
else
    enableState='off';    
end
set(get(handles.uipanel_filterParameters,'Children'),'Enable',enableState);


% --- Executes on selection change in popupmenu_method.
function popupmenu_method_Callback(hObject, eventdata, handles)

props=get(handles.popupmenu_method,{'UserData','Value'});
if strcmpi(props{1}{props{2}},'fastbem') || strcmpi(props{1}{props{2}},'fem')
    set(get(handles.uipanel_BEM,'Children'),'Enable','on');
    set(handles.useLcurve,'Enable','on');
    set(get(handles.groupCornerOptimal,'Children'),'Enable','on');
    set(handles.edit_tolx,'Enable','on');
    set(handles.text_tolx,'Enable','on');
    popupmenu_solMethodBEM_Callback(hObject,eventdata,handles);
else %when the method is fttc
    set(get(handles.uipanel_BEM,'Children'),'Enable','off');
    set(handles.edit_tolx,'Enable','off');
    set(handles.text_tolx,'Enable','off');
%     set(handles.useLcurve,'Value',true);
    set(handles.useLcurve,'Enable','on');
    set(get(handles.groupCornerOptimal,'Children'),'Enable','on');
end


% --- Executes on selection change in popupmenu_solMethodBEM.
function popupmenu_solMethodBEM_Callback(hObject, eventdata, handles)

props=get(handles.popupmenu_solMethodBEM,{'String','Value'});
if ~isempty(props{2}) && any(strcmpi(props{1}{props{2}},{'svd','gsvd'}))
    set(handles.edit_regParam,'Enable','off');
    set(handles.edit_LcurveFactor,'Enable','on');
else
    set(handles.edit_regParam,'Enable','on');
    set(handles.edit_LcurveFactor,'Enable','off');
end

if ~isempty(props{2}) && any(strcmpi(props{1}{props{2}},{'1NormReg','1NormRegLaplacian'}))
    set(handles.edit_tolx,'Enable','on');
    set(handles.text_tolx,'Enable','on');
else
    set(handles.edit_tolx,'Enable','off');
    set(handles.text_tolx,'Enable','off');
end



% --- Executes on button press in pushbutton_basisClassTblPath.
function pushbutton_basisClassTblPath_Callback(hObject, eventdata, handles)

[file, path] = uigetfile({'*.mat;*.MAT',...
    'Mat files (*.mat)'},...
    'Select the file containing the basis class lookup table');
if ~isequal(file, 0) && ~isequal(path, 0)
    try
        vars = whos('basisClassTbl','-file',[path file]);
        if numel(vars) ~= 1
            errordlg('Please select a file containing a valid basis class lookup table');
            return 
        end
    catch
        disp(['Using this path (' [path file] ') for initial input table path.'])
    end
    set(handles.edit_basisClassTblPath,'String',[path file]);
    
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over popupmenu_method.
function popupmenu_method_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function popupmenu_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in useLcurve.
function useLcurve_Callback(hObject, eventdata, handles)
% hObject    handle to useLcurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
useLcurve=get(handles.useLcurve,{'UserData','Value'});
if useLcurve{2}
%     set(handles.edit_regParam,'Enable','off');
    set(handles.checkbox_everyframe, 'Enable', 'on');
    set(handles.edit_LcurveFactor,'Enable','on');
    set(handles.lcorner,'Enable','on');
    set(handles.optimal,'Enable','on');
else
%     set(handles.edit_regParam,'Enable','on');
    set(handles.checkbox_everyframe, 'Enable', 'off');
    set(handles.edit_LcurveFactor,'Enable','off');
%     set(handles.groupCornerOptimal,'Enable','off');
    set(handles.lcorner,'Enable','off');
    set(handles.optimal,'Enable','off');
end

% Hint: get(hObject,'Value') returns toggle state of useLcurve


% --- Executes on button press in checkbox_lastToFirst.
function checkbox_lastToFirst_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_lastToFirst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_lastToFirst


% --- Executes on button press in setROIfromForcemap.
function setROIfromForcemap_Callback(hObject, eventdata, handles)
% hObject    handle to setROIfromForcemap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MD=handles.figure1.UserData.MD;
setROIfromForcemap(MD);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton_basisClassTblPath.
function pushbutton_basisClassTblPath_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_basisClassTblPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_tolx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tolx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tolx as text
%        str2double(get(hObject,'String')) returns contents of edit_tolx as a double


% --- Executes during object creation, after setting all properties.
function edit_tolx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tolx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function setROIfromForcemap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setROIfromForcemap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkbox_everyframe.
function checkbox_everyframe_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_everyframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_everyframe
