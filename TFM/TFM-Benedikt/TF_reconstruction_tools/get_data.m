% This program (regularized fourier transform traction force
% reconstruction) was produced at the University of Heidelberg, BIOMS
% group of Ulich Schwarz. It calculates traction from a gel displacement
% field.
%
% Benedikt Sabass 20-5-2007

function varargout = get_data(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @get_data_OpeningFcn, ...
                   'gui_OutputFcn',  @get_data_OutputFcn, ...
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


% --- Executes just before get_data is made visible.
function get_data_OpeningFcn(hObject, eventdata, handles, varargin)

handles.data.strainfile_name = 0;
handles.data.imagedir_name = 0;
handles.data.targetdir_name = 0;
handles.data.maskdir_name = [];

handles.data.young = 10000;
handles.data.poisson = 0.5;
handles.data.pix_durch_my = 0.067;
handles.data.strain = []; 

% Choose default command line output for get_data
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes get_data wait for user response (see UIRESUME)
% uiwait(handles.get_data_figure);


% --- Outputs from this function are returned to the command line.
function varargout = get_data_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on button press in strainfile_browse.
function strainfile_browse_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile('*.mat', 'Select a .mat file');
if ~isequal(filename,0)
    handles.data.strainfile_name = fullfile(pathname, filename);
    if size(handles.data.strainfile_name,2) > 55
        disp_string = ['... ',handles.data.strainfile_name(end-55:end)];
    else
        disp_string = handles.data.strainfile_name;
    end
    set(handles.strainfile_win,'string',disp_string);
end
guidata(hObject, handles);

% --- Executes on button press in imagedir_browse.
function imagedir_browse_Callback(hObject, eventdata, handles)
handles.data.imagedir_name = uigetdir('','Location of cell images');
if ~isequal(handles.data.imagedir_name,0)
    if size(handles.data.imagedir_name,2) > 55
        disp_string = ['... ',handles.data.imagedir_name(end-55:end)];
    else
        disp_string = handles.data.imagedir_name;
    end
    set(handles.imagedir_win,'String',disp_string);
end
guidata(hObject, handles);

% --- Executes on button press in targetdir_browse.
function targetdir_browse_Callback(hObject, eventdata, handles)
handles.data.targetdir_name = uigetdir('','Location to save results in');
if ~isequal(handles.data.targetdir_name,0)
    if size(handles.data.targetdir_name,2) > 55
        disp_string = ['... ',handles.data.targetdir_name(end-55:end)];
    else
        disp_string = handles.data.targetdir_name;
    end
    set(handles.targetdir_win,'String',disp_string);
end
guidata(hObject, handles);

% --- Executes on button press in maskdir_browse.
function maskdir_browse_Callback(hObject, eventdata, handles)
handles.data.maskdir_name = uigetdir('','Location of mask images');
if ~isequal(handles.data.maskdir_name,0)
    if size(handles.data.maskdir_name,2) > 55
        disp_string = ['... ',handles.data.maskdir_name(end-55:end)];
    else
        disp_string = handles.data.maskdir_name;
    end
    set(handles.maskdir_win,'String',disp_string);
end
guidata(hObject, handles);



function young_win_Callback(hObject, eventdata, handles)
written = str2double(get(hObject,'String'));
if ~isnan(written) && written > 0
    handles.data.young = written*1000;
else
    errordlg('The young modulus must be given as a positive number.','Error');
    set(handles.young_win,'String', num2str(handles.data.young/1000));
end
guidata(hObject, handles);

 
% --- Executes during object creation, after setting all properties.
function young_win_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function poisson_win_Callback(hObject, eventdata, handles)
written = str2double(get(hObject,'String'));
if ~isnan(written) && written >= 0 && written <= 0.5
    handles.data.poisson = written;
else
    errordlg('The poisson ratio must be between 0 and 0.5.','Error');
    set(handles.poisson_win,'String', num2str(handles.data.poisson));
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function poisson_win_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function pix_durch_my_win_Callback(hObject, eventdata, handles)
written = str2double(get(hObject,'String'));
if ~isnan(written) && written > 0
    handles.data.pix_durch_my = written;
else
    errordlg('Please enter a positive number here.','Error');
    set(handles.pix_durch_my_win,'String', num2str(handles.data.pix_durch_my));
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pix_durch_my_win_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)

if isempty(handles.data.targetdir_name) || isequal(handles.data.targetdir_name,0) || ...
        isempty(handles.data.imagedir_name) || isequal(handles.data.imagedir_name,0) || ...
        isempty(handles.data.strainfile_name) || isequal(handles.data.strainfile_name,0)
        errordlg('Some of the file/directory names have not been specified properly.','Error');
        return;
end

hilf  = load('-mat', handles.data.strainfile_name);
hilf_name = fieldnames(hilf);
if isempty(hilf) ||  all(strcmp(hilf_name, 'strain') == 0) 
        errordlg('Specified file with displacement data does not contain data or is empty.','Error');
        return;
else
    hilf_name = fieldnames(hilf.strain);
    if all(strcmp(hilf_name, 'pos') == 0) || all(strcmp(hilf_name, 'vec') == 0)
        errordlg('Specified file with displacement data does not contain the fields .pos/.vec.','Error');
        return;
    else
        handles.data.strain = hilf.strain;
    end
end

button = questdlg('Choose the traction reconstruction routine to use:','What now?','TRPF','Reg-FTTC','Reg-FTTC');

switch button
    case 'TRPF'
        TRPF_preview_settings(handles.data);
    case 'Reg-FTTC'
        preview_settings(handles.data);
end

delete(handles.get_data_figure);
return;


% --- Executes on button press in abort.
function abort_Callback(hObject, eventdata, handles)
delete(handles.get_data_figure);
return;



% --- Executes during object creation, after setting all properties.
function pushbutton5_CreateFcn(hObject, eventdata, handles)


% --- Executes during object deletion, before destroying properties.
function pushbutton5_DeleteFcn(hObject, eventdata, handles)




