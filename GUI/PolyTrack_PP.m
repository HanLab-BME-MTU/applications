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
defaultPostPro = struct('minusframes', 5, 'plusframes', 2, 'minimaltrack', 5, ...
                        'dragtail', 6, 'dragtailfile', 'trackmovie.mov', 'figureSize', []);

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

      % Call up the track linker function to generate the MPM matrix
      MPM = ptTrackLinker (M);
      
      % If there are totally empty rows in MPM, erase them
      % First find all unique rows unequal to zero
      [notZeroEntryRows, notZeroEntryCols] = find (MPM);
      notZeroEntryRows = unique (notZeroEntryRows);
      
      % Get these entries from MPM and keep them
      cleanedMPM = MPM (notZeroEntryRows,:);

      % Assign the loaded and cleaned up MPM to the handles struct
      handles.MPM = cleanedMPM;
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

% Load the cell properties file if it exists
if exist ('cellProps.mat','file')
   load('cellProps.mat');
   if exist('cellProps','var')
      handles.postpro.cellProps = cellProps;
   end
else
   h = errordlg('The file cells.mat does not exist. Please make sure it is present...');
   uiwait(h);          % Wait until the user presses the OK button
   return;
end

% Load the cluster properties file if it exists
if exist ('clusterProps.mat','file')
   load('clusterProps.mat');
   if exist('clusterProps','var')
      handles.postpro.clusterProps = clusterProps;
   end
else
   h = errordlg('The file clusters.mat does not exist. Please make sure it is present...');
   uiwait(h);          % Wait until the user presses the OK button
   return;
end

% Now that all the loading is done, we'll start the processing
cd (jobValPath);

% Counters to keep track of where we are
counter = 1;

% Here is where a new data subdirectory has to be created. Since 
% this has to be a unique name, a counter is used to find the
% next available unique directory name
while 1
   % Initialize the new data dir name
   newDirectoryName = ['data', num2str(counter)];
 
   % If it doesn't exist yet, create it in the results directory
   if ~exist (newDirectoryName, 'dir')
      mkdir (jobValPath, newDirectoryName);
        
      % Save it in the handles struct and tell the loop we're done
      handles.postpro.saveallpath = [jobValPath, newDirectoryName];
      break;
   end
   % Else the directory existed already so we increase the counter
   counter = counter + 1;
end

% Now we have to fill up the rest of the postpro structure with
% our previously found data and parameters
handles.selectedcells = [];
handles.postpro.imagepath = handles.jobvalues.imagedirectory;
handles.postpro.increment = handles.jobvalues.increment;
handles.postpro.maxdistpostpro = handles.jobvalues.maxsearch;
handles.postpro.plotfirstimg = handles.jobvalues.firstimage;
handles.postpro.plotlastimg = handles.jobvalues.lastimage;
handles.postpro.selectedcells = [];
handles.postpro.moviefirstimg = handles.jobvalues.firstimage;
handles.postpro.movielastimg = handles.jobvalues.lastimage;
handles.postpro.jobpath = jobValPath;
handles.postpro.imagename = handles.jobvalues.imagename;

% Update fields on the GUI with the latest values
set (handles.GUI_pp_jobpath_ed, 'String', jobValPath);
set (handles.GUI_pp_imagepath_ed, 'String', handles.postpro.imagepath);
set (handles.GUI_fm_saveallpath_ed, 'String', handles.postpro.saveallpath);
set (handles.GUI_ad_firstimage_ed, 'String', handles.postpro.plotfirstimg);
set (handles.GUI_ad_lastimage_ed, 'String', handles.postpro.plotlastimg);
set (handles.GUI_fm_movieimgone_ed, 'String', handles.postpro.moviefirstimg);
set (handles.GUI_fm_movieimgend_ed, 'String', handles.postpro.movielastimg);
set (handles.GUI_app_relinkdist_ed, 'String', handles.postpro.maxdistpostpro);
set (handles.GUI_app_minusframes_ed, 'String', handles.postpro.minusframes);
set (handles.GUI_app_plusframes_ed, 'String', handles.postpro.plusframes);
set (handles.GUI_app_minimaltrack_ed, 'String', handles.postpro.minimaltrack);
set (handles.GUI_fm_tracksince_ed, 'String', handles.postpro.dragtail);
set (handles.GUI_fm_filename_ed, 'String', handles.postpro.dragtailfile);

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
[filename, jobValPath] = uigetfile ({'*.tif','TIFF-files'},'Please select an image');

% And store this in the postpro struct
if exist ('jobValPath', 'file')
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

% This is to signal that the ptManualPostProcessJob function is called by the
% GUI manual processing button
handles.whichcallback = 1;

% Update handles structure
guidata(hObject, handles);

% Do the manual postprocessing 
ptManualPostProcessJob (hObject);

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

% Get the string value of 'frames before gap' from the GUI
num = get (hObject, 'String');

% Convert this value to a number
handles.postpro.minusframes = str2num (num);

% Update handles structure
guidata (hObject, handles);

% Update the field on the GUI
set (handles.GUI_app_minusframes_ed, 'String', handles.postpro.minusframes);

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
handles = guidata (hObject);

% Get the string value of 'frames after gap' from the GUI
num = get (hObject, 'String');

% Convert this value to a number
handles.postpro.plusframes = str2num (num);

% Update handles structure
guidata (hObject, handles);

% Update the field on the GUI
set (handles.GUI_app_plusframes_ed, 'String', handles.postpro.plusframes);

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
handles = guidata (hObject);

% Get the string value of 'maximum track distance' from the GUI
num = get (hObject, 'String');

% Convert this value to a number
handles.postpro.maxdistpostpro = str2num (num);

% Update handles structure
guidata (hObject, handles);

% Update the field on the GUI
set (handles.GUI_app_relinkdist_ed, 'String', handles.postpro.maxdistpostpro);

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
handles = guidata (hObject);

% Get the string value of 'minimum track distance' from the GUI
num = get(hObject,'String');

% Convert this value to a number
handles.postpro.minimaltrack = str2num (num);

% Update handles structure
guidata (hObject, handles);

% Update the field on the GUI
set (handles.GUI_app_minimaltrack_ed, 'String', handles.postpro.minimaltrack);

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_app_autopostpro_pb.
function GUI_app_autopostpro_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_app_autopostpro_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Get all the latest values from the GUI
saveAllPath = handles.postpro.saveallpath;
minTrackDistance = handles.postpro.minimaltrack;
maxDistance = handles.postpro.maxdistpostpro;
minusFrames = handles.postpro.minusframes;
plusFrames = handles.postpro.plusframes;
MPM = handles.MPM;

% Get the latest MPM matrix and filter it using the values on the GUI (eg eliminating
% tracks that are too short)
updatedMPM = ptTrackFilter (MPM, plusFrames, minusFrames, maxDistance, minTrackDistance, saveAllPath);

% Update the handles structure with the filtered MPM matrix
handles.MPM = updatedMPM;

% Update the handles structure
guidata (hObject, handles);

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

% --- Executes on button press in GUI_ad_perimeter_rb.
function GUI_ad_perimeter_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_perimeter_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_ad_perimeter_rb

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_ad_selectcells_pb.
function GUI_ad_selectcells_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_selectcells_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

% This is to signal that the function ptManualPostProcessJob is called by the
% GUI select cells button
handles.whichcallback = 2;

% Update handles structure
guidata(hObject, handles);

% Call the function to select cells
% AK: have to check whether this still works!
ptManualPostProcessJob (hObject)

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

% Get the number of the first image
num = get (hObject,'String');

% Assign it to the postpro structure
handles.postpro.plotfirstimg = str2num (num);

% Update handles structure
guidata (hObject, handles);

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
handles = guidata (hObject);

% Get the number of the last image
num = get (hObject,'String');

% Assign it to the postpro structure
handles.postpro.plotlastimg = str2num (num);

% Update handles structure
guidata(hObject, handles);

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_ad_analyze_pb.
function GUI_ad_analyze_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_analyze_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata (hObject);

% Assign the radiobutton values to the postpro struct
handles.postpro.cellareaplot = get (handles.GUI_ad_numberofthings_rb,'Value');

% Here is where the bulk of the graphing work is done; we give it the
% postpro structure and MPM matrix to work with
ptPlotCellValues (handles.postpro, handles.MPM);

% Only if the 'speed of cells' radiobutton is selected, the speed graphs have to be done
%if get (handles.GUI_ad_speed_rb, 'Value')
%   ptPlotSpeedValues (hObject);
%end

% This is a function which was done for the cell meeting poster
% Temporary and has to be replaced by a more structured function
%ptPoster (hObject);

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
handles = guidata (hObject);

% Get the value of 'First Image'
num = get (hObject, 'String');

% Assign it to the postpro structure
handles.postpro.moviefirstimg = str2num (num);

% Update handles structure
guidata(hObject, handles);

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
handles = guidata (hObject);

% Get the value of 'Last Image'
num = get (hObject, 'String');

% Assign it to the postpro structure
handles.postpro.movielastimg = str2num (num);

% Update handles structure
guidata(hObject, handles);

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
handles = guidata (hObject);

% Get the value of 'Tracks from last...'
num = get (hObject, 'String');

% Assign it to the postpro structure
handles.postpro.dragtail = str2num (num);

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
handles = guidata (hObject);

% Get the value of the save directory
path = get (hObject, 'String');

% Assign it to the postpro structure
handles.postpro.saveallpath = path;

% Update handles structure
guidata (hObject, handles);

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_fm_universalstudios_pb.
function GUI_fm_universalstudios_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_fm_universalstudios_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Start the function that will create the dragtail movie
ptMovieMaker (hObject);

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

viewPrepFigH = figure;

h = helpdlg ('Please resize the figure and close the figure window. The selected size will be used for the movie.');

set (h, 'Position', [320.2500 272.2500 297.7500 79.5000]);

uiwait (h);

handles.postpro.figureSize = get (viewPrepFigH, 'Position');
    
close (viewPrepFigH);

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
handles = guidata (hObject);

% Get the value of the 'time to be a cluster' field
num = get (hObject,'String');

% Assign it to the postpro structure
handles.postpro.mintimeclust = str2num(num);

% Update handles structure
guidata (hObject, handles);

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

%--------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_fm_filename_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_fm_filename_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%--------------------------------------------------------------------------

function GUI_fm_filename_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_fm_filename_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_fm_filename_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_fm_filename_ed as a double
handles = guidata (hObject);

% Get the value of 'Tracks from last...'
filename = get (hObject, 'String');

% Assign it to the postpro structure
handles.postpro.dragtailfile = filename;

% Update handles structure
guidata(hObject, handles);


