function varargout = PolyTrack_PP(varargin)
% PolyTrack_PP M-file for PolyTrack_PP.fig
%
% PolyTrack_PP  contains the callback functions for the Polytrack Post Processing 
%               (PP) GUI. Next to Matlab standard initialization code and functions, 
%               a number of callbacks have been implemented that control the 
%               post processing of the data files created by the Polytrack program.
%
% SYNOPSIS      varargout = PolyTrack_PP(varargin)
%
% INPUT         varargin (optional)
%
% OUTPUT        varargout
%
% Last Modified by A. Kerstens 09-Mar-2004

% All kinds of matlab initialization stuff; leave as is...
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_postprocess_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_postprocess_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

%----------------------------------------------------------------------------

% --- Executes just before GUI_postprocess is made visible.
function GUI_postprocess_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_postprocess (see VARARGIN)

% Choose default command line output for GUI_postprocess
handles.output = hObject;

% Create a default postprocessing structure that defines all the fields used in the program
defaultPostPro = struct('minusframes', 5, 'plusframes', 2, 'minimaltrack', 5,...
                        'dragtail', 6, 'figureSize', []);

% Assign the default postprocessing values to the GUI handle so it can be passed around
handles.defaultPostPro = defaultPostPro;

% Set the color of the gui
set(hObject,'Color',[0,0,0.627]);

% Update handles structure
guidata(hObject, handles);

%----------------------------------------------------------------------------

% --- Outputs from this function are returned to the command line.
function varargout = GUI_postprocess_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%----------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_pp_jobpath_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_pp_jobpath_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%----------------------------------------------------------------------------

function GUI_pp_jobpath_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_pp_jobpath_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_pp_jobpath_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_pp_jobpath_ed as a double

handles = guidata(hObject);

% Get the job path and name
jobPath = get(hObject,'String');

if ~exist(jobPath, 'file')
   h=errordlg('The selected job path does not exist. Please select another job path... ');
   uiwait(h);
   return
end

% Update handles structure
guidata(hObject, handles);

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_pp_jobbrowse_pb.
function GUI_pp_jobbrowse_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_pp_jobbrowse_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Update handles structure
handles = guidata(hObject);

% Start with the default post processing structure
handles.postpro = handles.defaultPostPro;

% Use the gui to let the user select a jobvalues.mat filename
[filename, jobValPath] = uigetfile({'*.mat','mat-files'},'Select jobvalues.mat or MPM.mat');

% Check that the user actually selected a valid file
if ~strcmp(filename, 'jobvalues.mat') & ~strcmp(filename, 'MPM.mat')
   h = errordlg('Please select a file named jobvalues.mat or MPM.mat...');
   uiwait(h);          % Wait until the user presses the OK button
   return
end

% Change directory to the selected path
cd (jobValPath)

% Load MPM.mat file if selected from the gui and if it exists
if strcmp (filename, 'MPM.mat')
   if exist ('MPM.mat', 'file')
      load (filename);
      handles.MPM = MPM;
   else
      h = errordlg('The file MPM.mat does not exist...');
      uiwait(h);          % Wait until the user presses the OK button
      return;
   end
   
   % Step back one directory
   cd ..
   jobValPath = [pwd, filesep];

elseif strcmp (filename, 'jobvalues.mat')
   % Load the M.mat file which should be present
   if exist ('M.mat', 'file')
      load ('M.mat');
      % Call up the track link function to generate the MPM matrix
      handles.MPM = trackLinker(M);
   else
      h = errordlg('The file M.mat does not exist. Please make sure it is present as well...');
      uiwait(h);          % Wait until the user presses the OK button
      return;
   end
end
    
% Load the jobvalues file
if exist('jobvalues.mat', 'file')
   load('jobvalues.mat');
   handles.jobvalues = jobvalues;
else
   h = errordlg('The file jobvalues.mat does not exist...');
   uiwait(h);          % Wait until the user presses the OK button
   return;
end

% Load the cell properties file if it exists (it should!)
if exist ('cellprops.mat','file')
   load('cellprops.mat');
   if exist('cellProps','var')
      handles.cellProps = cellProps;
   else
      handles.cellProps = cellprops;
   end
else
   h = errordlg('The file cellprops.mat does not exist. Please make sure it is present as well...');
   uiwait(h);          % Wait until the user presses the OK button
   return;
end

% Now that all the loading is done, we'll start the real stuff
cd (jobValPath)

% Counters to keep track of where we are
done = 0;
counter = 1;

% Here is where a new data subdirectory has to be created. Since 
% this has to be a unique name, a counter is used to find the
% next available unique directory name
while done == 0
   newdirname=[];
   newdirname=['data', num2str(counter)];
 
   % If it doesn't exist yet, create it in the results directory
   if ~exist(newdirname, 'dir')
      mkdir(jobValPath, newdirname);
        
      % Save it in the handles struct and tell the loop we're done
      handles.postpro.saveallpath = [jobValPath, newdirname];
      done = 1;
   end
   counter = counter + 1;
end

% Now we have to fill up the rest of the postpro structure with
% our previously found data and parameters
handles.selectedcells = [];
handles.postpro.imagepath = handles.jobvalues.imagedirectory;
handles.postpro.maxdistpostpro = handles.jobvalues.maxsearch;
handles.postpro.analfirstimg = handles.jobvalues.firstimage;
handles.postpro.anallastimg = handles.jobvalues.lastimage;
handles.postpro.selectedcells = [];
handles.postpro.moviefirstimg = handles.jobvalues.firstimage;
handles.postpro.movielastimg = handles.jobvalues.lastimage;

% Update fields on the GUI with the latest values
set(handles.GUI_pp_jobpath_ed, 'String', jobValPath);
set(handles.GUI_pp_imagepath_ed, 'String', handles.postpro.imagepath);
set(handles.GUI_fm_saveallpath_ed, 'String', handles.postpro.saveallpath);
set(handles.GUI_ad_firstimage_ed, 'String', handles.postpro.analfirstimg);
set(handles.GUI_ad_lastimage_ed, 'String', handles.postpro.anallastimg);
set(handles.GUI_fm_movieimgone_ed, 'String', handles.postpro.analfirstimg);
set(handles.GUI_fm_movieimgend_ed, 'String', handles.postpro.anallastimg);
set(handles.GUI_app_relinkdist_ed, 'String', handles.postpro.maxdistpostpro);
set(handles.GUI_app_minusframes_ed, 'String', handles.postpro.minusframes);
set(handles.GUI_app_plusframes_ed, 'String', handles.postpro.plusframes);
set(handles.GUI_app_minimaltrack_ed, 'String', handles.postpro.minimaltrack);
set(handles.GUI_fm_tracksince_ed, 'String', handles.postpro.dragtail);

% And update the gui handles struct
guidata(hObject, handles);

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_pp_imagebrowse_pb.
function GUI_pp_imagebrowse_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_pp_imagebrowse_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Let the user browse for an image file and path
[filename,jobValPath] = uigetfile ({'*.tif','TIFF-files'},'Please select an image');

% And store this in the postpro struct
if exist('jobValPath', 'file')
   cd (jobValPath);

   if ~exist('filename', 'file')
      h = errordlg('The image file does not exist. Please select another file...');
      uiwait(h);          % Wait until the user presses the OK button
      return;
   else
      handles.postpro.imagepath = jobValPath;
   end
else
   h = errordlg('The image directory does not exist. Please select another directory...');
   uiwait(h);          % Wait until the user presses the OK button
   return;
end

% Update handles structure
guidata(hObject, handles);

%----------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_pp_imagepath_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_pp_imagepath_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%----------------------------------------------------------------------------

function GUI_pp_imagepath_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_pp_imagepath_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_pp_imagepath_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_pp_imagepath_ed as a double
handles = guidata(hObject);

% Assign the entered path to the postpro struct
path = get(hObject,'String');
handles.postpro.imagepath = path;

% Update handles structure
guidata(hObject, handles);

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_pp_manuelpostpro_pb.
function GUI_pp_manuelpostpro_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_pp_manuelpostpro_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);


handles.whichcallback = 1;

% Update handles structure
guidata(hObject, handles);



manuelpostpro(hObject);

%----------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_app_minusframes_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_app_minusframes_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%----------------------------------------------------------------------------

function GUI_app_minusframes_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_app_minusframes_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_app_minusframes_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_app_minusframes_ed as a double

handles = guidata(hObject);


numb=get(hObject,'String');


handles.postpro.minusframes= str2num(numb);


% Update handles structure
guidata(hObject, handles);

set(handles.GUI_app_minusframes_ed,'String',handles.postpro.minusframes);

%----------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_app_plusframes_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_app_plusframes_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%----------------------------------------------------------------------------

function GUI_app_plusframes_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_app_plusframes_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_app_plusframes_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_app_plusframes_ed as a double

handles = guidata(hObject);


numb=get(hObject,'String');


handles.postpro.plusframes= str2num(numb);


% Update handles structure
guidata(hObject, handles);

set(handles.GUI_app_plusframes_ed,'String',handles.postpro.plusframes);

%----------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_app_relinkdist_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_app_relinkdist_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%----------------------------------------------------------------------------

function GUI_app_relinkdist_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_app_relinkdist_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_app_relinkdist_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_app_relinkdist_ed as a double


handles = guidata(hObject);


numb=get(hObject,'String');


handles.postpro.maxdistpostpro= str2num(numb);


% Update handles structure
guidata(hObject, handles);

set(handles.GUI_app_relinkdist_ed,'String',handles.postpro.maxdistpostpro);

%----------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_app_minimaltrack_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_app_minimaltrack_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%----------------------------------------------------------------------------

function GUI_app_minimaltrack_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_app_minimaltrack_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_app_minimaltrack_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_app_minimaltrack_ed as a double


handles = guidata(hObject);


numb=get(hObject,'String');


handles.postpro.minimaltrack= str2num(numb);


% Update handles structure
guidata(hObject, handles);

set(handles.GUI_app_minimaltrack_ed,'String',handles.postpro.minimaltrack);

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_app_autopostpro_pb.
function GUI_app_autopostpro_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_app_autopostpro_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

saveallpath = handles.postpro.saveallpath;
minimaltrack = handles.postpro.minimaltrack;
maxdistpostpro = handles.postpro.maxdistpostpro ;
minusframes = handles.postpro.minusframes;
plusframes = handles.postpro.plusframes;

MPM = handles.MPM;
MPM = weedout (MPM, plusframes, minusframes, maxdistpostpro, minimaltrack, saveallpath);
handles.MPM = MPM;

guidata(hObject, handles);

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_ad_speed_rb.
function GUI_ad_speed_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_speed_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_ad_speed_rb

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_ad_trianghist_rb.
function GUI_ad_trianghist_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_trianghist_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_ad_trianghist_rb

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_ad_numberofthings_rb.
function GUI_ad_numberofthings_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_numberofthings_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_ad_numberofthings_rb

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_ad_areas_rb.
function GUI_ad_areas_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_areas_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_ad_areas_rb

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_ad_undefined3_rb.
function GUI_ad_undefined3_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_undefined3_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_ad_undefined3_rb

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_ad_perimeter_rb.
function GUI_ad_perimeter_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_perimeter_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_ad_perimeter_rb

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_ad_undefined1_rb.
function GUI_ad_undefined1_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_undefined1_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_ad_undefined1_rb

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_ad_undefined2_rb.
function GUI_ad_undefined2_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_undefined2_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_ad_undefined2_rb

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_ad_selectcells_pb.
function GUI_ad_selectcells_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_selectcells_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);


handles.whichcallback = 2;

% Update handles structure
guidata(hObject, handles);

manuelpostpro(hObject)

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_ad_selected_cb.
function GUI_ad_selected_cb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_selected_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_ad_selected_cb

%----------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_ad_firstimage_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_ad_firstimage_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%----------------------------------------------------------------------------

function GUI_ad_firstimage_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_firstimage_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_ad_firstimage_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_ad_firstimage_ed as a double
handles = guidata(hObject);

numb=get(hObject,'String');



handles.postpro.firstlastimg= str2num(numb);

% Update handles structure
guidata(hObject, handles);

%----------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_ad_lastimage_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_ad_lastimage_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%----------------------------------------------------------------------------

function GUI_ad_lastimage_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_lastimage_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_ad_lastimage_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_ad_lastimage_ed as a double
handles = guidata(hObject);

numb=get(hObject,'String');



handles.postpro.anallastimg= str2num(numb);

% Update handles structure
guidata(hObject, handles);

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_ad_analyze_pb.
function GUI_ad_analyze_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_analyze_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


stuffplotter(hObject);

handles = guidata(hObject);
if get(handles.GUI_ad_speed_rb,'Value')
    speed (hObject);
end

ptPoster (hObject);

%----------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_fm_movieimgone_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_fm_movieimgone_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%----------------------------------------------------------------------------

function GUI_fm_movieimgone_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_fm_movieimgone_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_fm_movieimgone_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_fm_movieimgone_ed as a double

%----------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_fm_movieimgend_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_fm_movieimgend_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%----------------------------------------------------------------------------

function GUI_fm_movieimgend_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_fm_movieimgend_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_fm_movieimgend_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_fm_movieimgend_ed as a double

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_fm_incloriginal_rb.
function GUI_fm_incloriginal_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_fm_incloriginal_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_fm_incloriginal_rb

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_fm_inclcentromers_rb.
function GUI_fm_inclcentromers_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_fm_inclcentromers_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_fm_inclcentromers_rb

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_fm_inclundefined_rb.
function GUI_fm_inclundefined_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_fm_inclundefined_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_fm_inclundefined_rb

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_fm_incltracks_rb.
function GUI_fm_incltracks_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_fm_incltracks_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_fm_incltracks_rb

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_fm_inclbody_rb.
function GUI_fm_inclbody_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_fm_inclbody_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_fm_inclbody_rb

%----------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_fm_tracksince_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_fm_tracksince_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%----------------------------------------------------------------------------

function GUI_fm_tracksince_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_fm_tracksince_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_fm_tracksince_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_fm_tracksince_ed as a double

handles = guidata(hObject);

numb = get(hObject,'String');

handles.postpro.dragtail= str2num(numb);

% Update handles structure
guidata(hObject, handles);

%----------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_fm_saveallpath_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_fm_saveallpath_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%----------------------------------------------------------------------------

function GUI_fm_saveallpath_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_fm_saveallpath_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_fm_saveallpath_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_fm_saveallpath_ed as a double

handles = guidata(hObject);

path = get(hObject,'String');

handles.postpro.saveallpath= path;

% Update handles structure
guidata(hObject, handles);

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_fm_universalstudios_pb.
function GUI_fm_universalstudios_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_fm_universalstudios_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

movieMaker(hObject) 

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_fm_browsesaveallpath_pb.
function GUI_fm_browsesaveallpath_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_fm_browsesaveallpath_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_fm_moviesize_pb.
function GUI_fm_moviesize_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_fm_moviesize_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  viewPrepFigH = figure%('NumberTitle','off','Name','select view');
%     plot3(allSpotCoords(:,1),allSpotCoords(:,2),allSpotCoords(:,3),'.','MarkerSize',16);
%     line(allSpotCoords(:,1),allSpotCoords(:,2),allSpotCoords(:,3),'Color','k');
%     set(gca,'XLim',axesXLim,'YLim',axesYLim,'ZLim',axesZLim,'Box','on','PlotBoxAspectRatio', boxAspectRatio);
%     view(azimuth,elevation);
%     
    h = helpdlg('Please choose a figure size, THEN close this window. This size will be used for the movie');
    %place dialogbox below figure
    %...somehow I don't get it with the positioning of the figure
    %     figPos = get(viewPrepFigH,'Position');
    %     dlgPos = get(h,'Position');
    %     dlgPos(2) = figPos(2) - dlgPos(4);
    set(h,'Position',[320.2500  272.2500  297.7500   79.5000]);
    uiwait(h);
%     
%     [azimuth,elevation] = view;
%     if azimuth < 0
%         azimuth = azimuth +360;
%     end
%     disp(['az and el in case you have to retry: ',num2str(azimuth),' ',num2str(elevation)])
%     
    handles.postpro.figureSize = get(viewPrepFigH,'Position');
    
    close(viewPrepFigH);

%----------------------------------------------------------------------------

% % % % % % % 
% % % % % % % % --- Executes on button press in GUI_app_autopostpro_cb.
% % % % % % % function GUI_app_autopostpro_cb_Callback(hObject, eventdata, handles)
% % % % % % % % hObject    handle to GUI_app_autopostpro_cb (see GCBO)
% % % % % % % % eventdata  reserved - to be defined in a future version of MATLAB
% % % % % % % % handles    structure with handles and user data (see GUIDATA)
% % % % % % % 
% % % % % % % % Hint: get(hObject,'Value') returns toggle state of GUI_app_autopostpro_cb
% % % % % % % 
% % % % % % % 

%----------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_ad_mintimeclust_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_ad_mintimeclust_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%----------------------------------------------------------------------------

function GUI_ad_mintimeclust_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_mintimeclust_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_ad_mintimeclust_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_ad_mintimeclust_ed as a double



handles = guidata(hObject);

numb = get(hObject,'String');

handles.postpro.mintimeclust= str2num(numb);

% Update handles structure
guidata(hObject, handles);

%----------------------------------------------------------------------------

function file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%----------------------------------------------------------------------------

function exit_item_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Destroy the GUI
delete(handles.polyTrack_PP_mainwindow);

