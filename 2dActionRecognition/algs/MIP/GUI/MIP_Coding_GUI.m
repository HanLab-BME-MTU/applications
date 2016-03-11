function varargout = MIP_Coding_GUI(varargin)
% MIP_CODING_GUI MATLAB code for MIP_Coding_GUI.fig
%      MIP_CODING_GUI, by itself, creates a new MIP_CODING_GUI or raises the existing
%      singleton*.
%
%      H = MIP_CODING_GUI returns the handle to a new MIP_CODING_GUI or the handle to
%      the existing singleton*.
%
%      MIP_CODING_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MIP_CODING_GUI.M with the given input arguments.
%
%      MIP_CODING_GUI('Property','Value',...) creates a new MIP_CODING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MIP_Coding_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MIP_Coding_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MIP_Coding_GUI

% Last Modified by GUIDE v2.5 24-Oct-2012 20:18:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MIP_Coding_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MIP_Coding_GUI_OutputFcn, ...
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


% --- Executes just before MIP_Coding_GUI is made visible.
function MIP_Coding_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MIP_Coding_GUI (see VARARGIN)

% Choose default command line output for MIP_Coding_GUI
handles.output = hObject;

%Store GUI data on same handle as output.
handles.GUIData = handles.output;
%configure params files names:
default_params_file = fullfile('.','get_feature_params_default.m');
params_file = regexprep(default_params_file,'_default','_runtime');
copyfile(default_params_file,params_file);
%set to appdata:
setappdata(handles.GUIData,'default_params_file',default_params_file);
setappdata(handles.GUIData,'params_file',params_file);

%addpath:
run(params_file);
if(~isdeployed)
   addpath(genpath(params.GUI.add_path)); 
end

% Update handles structure
guidata(hObject, handles);

%update GUI fields:
set(handles.start_dir,'String',params.GUI.def_start_dir);
set(handles.out_dir,'String',params.GUI.def_out_dir);
%find input movies:
set(handles.Status,'String','Looking for Movies in selected dir');
drawnow;
find_files(handles);

% UIWAIT makes MIP_Coding_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MIP_Coding_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Browse_input.
function Browse_input_Callback(hObject, eventdata, handles)
% hObject    handle to Browse_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    cur_start_dir = get(handles.start_dir,'String');
    if(isempty(cur_start_dir))
        cur_start_dir = fullfile('..','..','Data');
    end
    dialog_title = 'Please Select parent directory for the data set';
    folder_name = uigetdir(cur_start_dir,dialog_title);
    set(handles.start_dir,'String',folder_name);
    set(handles.Status,'String','Looking for Movies in selected dir');
    drawnow;
    find_files(handles);

% --- Executes on selection change in files_list.
function files_list_Callback(hObject, eventdata, handles)
% hObject    handle to files_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns files_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from files_list


% --- Executes during object creation, after setting all properties.
function files_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to files_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in change_params.
function change_params_Callback(hObject, eventdata, handles)
% hObject    handle to change_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    %Open params file for edit:
    params_file = getappdata(handles.GUIData,'params_file');
    [status, result] = system(['notepad ',params_file,'']);
    set(handles.Status,'String','Params File updated');
    

% --- Executes on button press in run_coding.
function run_coding_Callback(hObject, eventdata, handles)
% hObject    handle to run_coding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    %get Files to run on:
    database_info = getappdata(handles.GUIData,'database_info');
    
    selection_str = get(get(handles.run_selection,'SelectedObject'),'String');
    switch selection_str
        case 'All files'
            database_info = database_info;%run on all files
        case 'Selected files'
            values = get(handles.files_list,'Value');
            database_info.db_data = database_info.db_data(values,:);
        otherwise
            return;%cannot get here...
    end
    set(handles.Status,'String','Running coding...');
    drawnow;
    %run the coding function:
    features_dir = get(handles.out_dir,'String');
    params_file = getappdata(handles.GUIData,'params_file');
    run(params_file);
    run_MIP_coding_on_files(features_dir,params,database_info,params.GUI.add_path);
    set(handles.Status,'String','Running coding... Done.');


% --- Executes when selected object is changed in run_selection.
function run_selection_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in run_selection 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Browse_output.
function Browse_output_Callback(hObject, eventdata, handles)
% hObject    handle to Browse_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    cur_out_dir = get(handles.out_dir,'String');
    if(isempty(cur_out_dir))
        cur_out_dir = fullfile('..','..','Data');
    end
    dialog_title = 'Please select output directory';
    folder_name = uigetdir(cur_out_dir,dialog_title);
    set(handles.out_dir,'String',folder_name);


% --- Executes on button press in readme.
function readme_Callback(hObject, eventdata, handles)
% hObject    handle to readme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    readmepath = fullfile('..','..','ReadMe');
    [status, result] = system(['notepad ',readmepath,'']);
    set(handles.Status,'String','Thank you for using MIP :)');
