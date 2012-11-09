function varargout = panda(varargin)
% PANDA M-file for panda.fig
%      PANDA is a software to calculate cell protrusion/retraction
%      based on constrainted optimization in 2D time-series images. 
%
%      PANDA, by itself, creates a new PANDA or raises the existing
%      singleton*.
%
%      H = PANDA returns the handle to a new PANDA or the handle to
%      the existing singleton*.
%
%      PANDA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PANDA.M with the given input arguments.
%
%      PANDA('Property','Value',...) creates a new PANDA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pandaPanelCarton_OpeningFunction gets called.
%      An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to panda_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES, protrusionAnalysis, generatePRMap

% Last Modified by GUIDE v2.5 02-Jun-2009 12:43:33
% Author: Shann-Ching Chen, LCCB

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @panda_OpeningFcn, ...
                   'gui_OutputFcn',  @panda_OutputFcn, ...
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

% --------------------------------------------------------------------


% --- Executes just before panda is made visible.
function panda_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to panda (see VARARGIN)

handles.directory_name = pwd;
maskinfo = update_mask_num(hObject, handles);
handles.maskinfo = maskinfo;
handles.timevalue = 1;
handles.resolutionvalue = 1;
handles.segvalue = 30;
handles.dl_rate = 10;
% Choose default command line output for panda
handles.output = hObject;
handles.batch_processing = 0;
handles.CC = 0;

imagePath = ['bannerv1p1.jpg'];
handles.banneraxes = imread(imagePath);
axes(handles.axes_logo);
image(handles.banneraxes);
set(handles.axes_logo, 'Visible', 'off');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes panda wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = panda_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% ------------------------------------------------------------
% Update number of masks under the handles.directory_name
% ------------------------------------------------------------
function maskinfo = update_mask_num(hObject, handles)

Identifier = get(handles.idedit,'String');
FileType_list = get(handles.popupmenu1,'String');
FileType = FileType_list{get(handles.popupmenu1,'Value')};

fullmaskpath = [handles.directory_name filesep Identifier FileType];
maskdir = dir(fullmaskpath);
mask_num = length(maskdir);


maskinfo.fullmaskpath = fullmaskpath;
maskinfo.maskdir = maskdir;
maskinfo.mask_num = mask_num;
maskinfo.FileType = [Identifier FileType];

set(handles.imgnumtext, 'String', mask_num);

if(mask_num > 0)
    set(handles.loadbutton1, 'Enable', 'on'); 
    set(handles.outputdirectoryedit, 'String', handles.directory_name);
    set(handles.slider1, 'Enable', 'on');        
else
    set(handles.loadbutton1, 'Enable', 'off');
    set(handles.playbutton, 'Enable', 'off');
    set(handles.analysisbutton, 'Enable', 'off');
    set(handles.mapbutton, 'Enable', 'off');
    set(handles.browsebutton2, 'Enable', 'off');
    set(handles.outputdirectoryedit, 'Enable', 'off');
    set(handles.timeedit, 'Enable', 'off');
    set(handles.resolutionedit, 'Enable', 'off');
    set(handles.segedit, 'Enable', 'off');
    set(handles.slider1, 'Enable', 'off'); 
    set(handles.menu_dl, 'Enable', 'off');          
end

% --- Executes on button press in browsebutton1.
function browsebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to browsebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% imagePath = ['banner.jpg'];
% handles.banneraxes = imread(imagePath);
% 
% [info.Height info.Width info.Channel] = size(handles.banneraxes);
% % position = get(handles.banneraxes, 'Position');
% % axes(handles.banneraxes);
% image(handles.banneraxes)
% set(handles.banneraxes, ...
%     'Visible', 'off', ...
%     'Units', 'pixels', ...
%     'Position', [50 530 info.Width info.Height]);
% 



directory_name = uigetdir;
if directory_name == 0
    return;
end  
handles.directory_name = directory_name;
handles.result_directory_name = directory_name;
set(handles.inputdirectoryedit,'String',handles.directory_name);
maskinfo = update_mask_num(hObject, handles);
handles.maskinfo = maskinfo;

guidata(hObject, handles);



% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
maskinfo = update_mask_num(hObject, handles);
handles.maskinfo = maskinfo;

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function inputdirectoryedit_Callback(hObject, eventdata, handles)
% hObject    handle to inputdirectoryedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputdirectoryedit as text
%        str2double(get(hObject,'String')) returns contents of inputdirectoryedit as a double
new_directory_name = get(hObject,'String');

if isdir(new_directory_name)
    handles.directory_name = new_directory_name;
    set(hObject,'String', handles.directory_name);      
    maskinfo = update_mask_num(hObject, handles);
    handles.maskinfo = maskinfo;
    guidata(hObject, handles);
else
    errordlg([new_directory_name ' is not a valid directory name'],'Bad Input','modal')
    set(hObject,'String', handles.directory_name);  
end


% --- Executes during object creation, after setting all properties.
function inputdirectoryedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputdirectoryedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',pwd);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
mask_idx = get(handles.slider1,'Value');
mask_idx = round(mask_idx);
set(handles.slider1,'Value',mask_idx);
thismask = imread([handles.directory_name filesep handles.maskinfo.maskdir(mask_idx).name]);
axes(handles.axes1);
imshow(thismask, []);
title(['Mask index:' num2str(mask_idx) '    Mask Name:' handles.maskinfo.maskdir(mask_idx).name]);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function idedit_Callback(hObject, eventdata, handles)
% hObject    handle to idedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of idedit as text
%        str2double(get(hObject,'String')) returns contents of idedit as a double
maskinfo = update_mask_num(hObject, handles);
handles.maskinfo = maskinfo;

% Update handles structure
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function idedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to idedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','*');


% --- Executes on button press in loadbutton1.
function loadbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to loadbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mask_idx = 1;
slider_step(1) = 1/handles.maskinfo.mask_num;
slider_step(2) = 1/handles.maskinfo.mask_num;
thismask = imread([handles.directory_name filesep handles.maskinfo.maskdir(mask_idx).name]);
set(handles.slider1, 'Visible', 'on','Max',handles.maskinfo.mask_num,...
    'Value',1,'Min',1,'sliderstep',slider_step);

axes(handles.axes1);    imshow(thismask, []);
title(['Mask index:' num2str(mask_idx) '    Mask Name:' handles.maskinfo.maskdir(mask_idx).name]);
set(handles.playbutton, 'Visible', 'on', 'Enable', 'on');

set(handles.analysisbutton, 'Enable', 'on');
set(handles.mapbutton, 'Enable', 'on');
set(handles.browsebutton2, 'Enable', 'on');
set(handles.outputdirectoryedit, 'Enable', 'on');
set(handles.timeedit, 'Enable', 'on');
set(handles.resolutionedit, 'Enable', 'on');
set(handles.segedit, 'Enable', 'on');
set(handles.menu_dl, 'Enable', 'on');          
set(handles.start_position_setting, 'Enable', 'on');          

handles.result_directory_name = get(handles.outputdirectoryedit,'String');
handles.rowSize = size(thismask,1);
handles.colSize = size(thismask,2);
% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in browsebutton2.
function browsebutton2_Callback(hObject, eventdata, handles)
% hObject    handle to browsebutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result_directory_name = uigetdir;
if result_directory_name == 0
    return;
end  
handles.result_directory_name = result_directory_name;

set(handles.outputdirectoryedit,'String',handles.result_directory_name);
guidata(hObject, handles); 



function outputdirectoryedit_Callback(hObject, eventdata, handles)
% hObject    handle to outputdirectoryedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outputdirectoryedit as text
%        str2double(get(hObject,'String')) returns contents of outputdirectoryedit as a double

result_directory_name = get(hObject,'String');
if isdir(new_directory_name)
    handles.result_directory_name = result_directory_name;
    guidata(hObject, handles);
else
    errordlg([new_directory_name ' is not a valid directory name'],'Bad Input','modal')
    set(hObject,'String', handles.result_directory_name);  
end



% --- Executes during object creation, after setting all properties.
function outputdirectoryedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputdirectoryedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',pwd);




function timeedit_Callback(hObject, eventdata, handles)
% hObject    handle to timeedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timeedit as text
%        str2double(get(hObject,'String')) returns contents of timeedit as a double

val = str2double(get(hObject,'String'));
if isnumeric(val) & val > 0
    handles.timevalue = val;    
else    
    errordlg('You must enter a positive numeric value','Bad Input','modal')
    set(hObject,'String', handles.timevalue);      
end
guidata(hObject, handles); 


% --- Executes during object creation, after setting all properties.
function timeedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function resolutionedit_Callback(hObject, eventdata, handles)
% hObject    handle to resolutionedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resolutionedit as text
%        str2double(get(hObject,'String')) returns contents of resolutionedit as a double
val = str2double(get(hObject,'String'));
if isnumeric(val) & val > 0
    handles.resolutionvalue = val;    
else    
    errordlg('You must enter a positive numeric value','Bad Input','modal')
    set(hObject,'String', handles.resolutionvalue);      
end
guidata(hObject, handles); 


% --- Executes during object creation, after setting all properties.
function resolutionedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resolutionedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function segedit_Callback(hObject, eventdata, handles)
% hObject    handle to segedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of segedit as text
%        str2double(get(hObject,'String')) returns contents of segedit as a double
val = round(str2num(get(hObject,'String')));
if isnumeric(val) & val >= 10 & val <= 100
    handles.segvalue = val;
    set(hObject,'String', handles.segvalue);          
else    
    errordlg('You must enter a value between 10 and 100','Bad Input','modal')
    set(hObject,'String', handles.segvalue);      
end
guidata(hObject, handles); 


% --- Executes during object creation, after setting all properties.
function segedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to segedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in analysisbutton.
function analysisbutton_Callback(hObject, eventdata, handles)
% hObject    handle to analysisbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.statusbar, 'String', 'Status: Start Analysis');

handles.outdir = [handles.result_directory_name filesep 'analysis_dl' num2str(handles.dl_rate)]; 

OK = protrusionAnalysis(handles);
if OK == 0
    set(handles.statusbar, 'String', '');
else
    set(handles.statusbar, 'String', 'Status: Analysis is done.');
end


% --- Executes on button press in mapbutton.
function mapbutton_Callback(hObject, eventdata, handles)
% hObject    handle to mapbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.statusbar, 'String', 'Status: Examine Activity Map');
if handles.CC > 1
    handles.outdir = [handles.result_directory_name filesep 'analysis_CC2_dl' num2str(handles.dl_rate)]; 
else
    handles.outdir = [handles.result_directory_name filesep 'analysis_dl' num2str(handles.dl_rate)]; 
end
handles.outdir = [handles.result_directory_name filesep 'analysis_dl' num2str(handles.dl_rate)]; 
generatePRMap(handles);

% --------------------------------------------------------------------
function menu_close_Callback(hObject, eventdata, handles)
% hObject    handle to menu_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
user_response=questdlg('Do you want to exit?', ...
                    'Confirm Dialog', ...
                    'Yes','No','Yes');
                
switch user_response
case {'No'}
    % take no action
case 'Yes'
    % Prepare to close GUI application window
    %                  .
    %                  .
    %                  .
    delete(handles.figure1)
end


% --- Executes on button press in playbutton.
function playbutton_Callback(hObject, eventdata, handles)
% hObject    handle to playbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for mask_idx = 1:handles.maskinfo.mask_num
    thismask = imread([handles.directory_name filesep handles.maskinfo.maskdir(mask_idx).name]);
    axes(handles.axes1);
    imshow(thismask, []);
    set(handles.slider1,'Value',mask_idx);    
    title(['Mask index:' num2str(mask_idx) '    Mask Name:' handles.maskinfo.maskdir(mask_idx).name]);    
    pause(0.1);
end



% --- Executes on selection change in menu_dl.
function menu_dl_Callback(hObject, eventdata, handles)
% hObject    handle to menu_dl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns menu_dl contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menu_dl

%handles.dl_rate = get(hObject,'Value');
idx = get(hObject,'Value');
handles.dl_rate = idx*5+5;
set(handles.text8,'String',[filesep 'analysis_dl' num2str(handles.dl_rate) filesep]); 
if handles.CC > 1
    set(handles.text8,'String',[filesep 'analysis_CC2_dl' num2str(handles.dl_rate) filesep]);     
else
    set(handles.text8,'String',[filesep 'analysis_dl' num2str(handles.dl_rate) filesep]);         
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function menu_dl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_dl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in selectROI.
function selectROI_Callback(hObject, eventdata, handles)
% hObject    handle to selectROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CropPanel;



% --- Executes on button press in pushbutton_seg.
function pushbutton_seg_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SegPanel;


% --- Executes during object creation, after setting all properties.
function axes_logo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_logo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate banneraxes



% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1




% --------------------------------------------------------------------
function menu_Callback(hObject, eventdata, handles)
% hObject    handle to menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function start_position_setting_Callback(hObject, eventdata, handles)
% hObject    handle to start_position_setting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[parameters] = pandaParameter('myhandles',handles);
handles.ORIENT_CELL = parameters.ORIENT_CELL;
handles.BESTSTART   = parameters.BESTSTART;
handles.TOLERANCE   = parameters.TOLERANCE;

%[handles.ORIENT_CELL handles.BESTSTART handles.TOLERANCE]
guidata(hObject, handles);


% --- Executes on button press in centroid_correction.
function centroid_correction_Callback(hObject, eventdata, handles)
% hObject    handle to centroid_correction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of centroid_correction

if get(hObject,'Value');
    % with centroid correction
    CC = 2;
    set(handles.text8,'String',[filesep 'analysis_CC2_dl' num2str(handles.dl_rate) filesep]);     
else
    CC = 1;
    set(handles.text8,'String',[filesep 'analysis_dl' num2str(handles.dl_rate) filesep]);         
end
handles.CC = CC;
guidata(hObject, handles);
