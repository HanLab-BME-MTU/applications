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
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Jun 04          Cleaned up source and renamed file
% Andre Kerstens        Jul 04          Added size of movie to ptPostpro

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
                        'dragtail', 6, 'dragtailfile', 'trackmovie', 'figureSize', [], ...
                        'multFrameVelocity', 1, 'binsize', 4, 'mmpixel', 0.639, 'timeperframe', 300, ...
                        'movietype', 1);

% Assign the default postprocessing values to the GUI handle so it can be passed around
handles.defaultPostPro = defaultPostPro;

% Set the color of the gui
set(hObject,'Color',[0,0,0.627]);

% Update handles structure
guidata(hObject, handles);

% Turn the resize and int conversion (matlab 7) warnings off
iptsetpref ('TrueSizeWarning', 'off');

% For matlab 7 turn int conversion warnings off
matlabVersion = version;
if (matlabVersion(1) == '7')
  warning off all;
end

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

% Determine image file path from jobvalues path
cd (jobValPath); cd ('..');
imageFilePath = pwd;

% Change directory to the selected path
cd (jobValPath);

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
if exist ('jobvalues.mat', 'file')
   load ('jobvalues.mat');
   handles.jobvalues = jobvalues;
else
   % It might be one directory back if it is a processed MPM file
   cd ..;
   if exist ('jobvalues.mat', 'file')
      load ('jobvalues.mat');
      handles.jobvalues = jobvalues;
   else
      h = errordlg ('The file jobvalues.mat does not exist...');
      uiwait(h);          % Wait until the user presses the OK button
      return;
   end
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

% Store the size of the image
cd (imageFilePath);
tempImage = imreadnd2 (handles.jobvalues.imagename, 0, handles.jobvalues.intensityMax);
[rows, cols] = size (tempImage);
handles.postpro.rowsize = rows;
handles.postpro.colsize = cols;

% Now we have to fill up the rest of the postpro structure with
% our previously found data and parameters
handles.selectedcells = [];
handles.postpro.imagepath = imageFilePath;
handles.postpro.increment = handles.jobvalues.increment;
handles.postpro.firstimg = handles.jobvalues.firstimage;
handles.postpro.lastimg = handles.jobvalues.lastimage;
handles.postpro.maxdistpostpro = handles.jobvalues.maxsearch;
handles.postpro.plotfirstimg = handles.jobvalues.firstimage;
handles.postpro.plotlastimg = handles.jobvalues.lastimage;
handles.postpro.selectedcells = [];
handles.postpro.moviefirstimg = handles.jobvalues.firstimage;
handles.postpro.movielastimg = handles.jobvalues.lastimage;
handles.postpro.jobpath = jobValPath;
handles.postpro.imagename = handles.jobvalues.imagename;
handles.postpro.imagenameslist = handles.jobvalues.imagenameslist;
handles.postpro.intensitymax = handles.jobvalues.intensityMax;
handles.postpro.maxdistance = handles.jobvalues.maxsearch;

% These are new additions which won't be in the older jobs
if ~isempty (handles.jobvalues.timeperframe)
   handles.postpro.timeperframe = handles.jobvalues.timeperframe;
end

if ~isempty (handles.jobvalues.mmpixel)
   handles.postpro.mmpixel = handles.jobvalues.mmpixel;
end

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
set (handles.multFrameVelocity, 'String', handles.postpro.multFrameVelocity);
set (handles.GUI_ad_binsize_ed, 'String', handles.postpro.binsize);
set (handles.pp_firstframe, 'String', handles.postpro.firstimg);
set (handles.pp_lastframe, 'String', handles.postpro.lastimg);
set (handles.pp_increment, 'String', handles.postpro.increment);
set (handles.GUI_mmpixel_ed, 'String', handles.postpro.mmpixel);
set (handles.GUI_frameinterval_ed, 'String', handles.postpro.timeperframe);
set (handles.GUI_movietype_avi_rb, 'Value', 1);
set (handles.GUI_movietype_qt_rb, 'Value', 0);

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

% Check that the plot frame values provided are in range
if (handles.postpro.plotfirstimg < handles.postpro.firstimg) | ...
   (handles.postpro.plotlastimg > handles.postpro.lastimg)
   h = errordlg('Plot start and end frame value are out of range. Please reenter values...');
   uiwait(h);          % Wait until the user presses the OK button
   return;
end

% Assign the radiobutton values to the postpro struct
handles.postpro.cellclusterplot = get (handles.checkbox_clustercellstats,'Value');
handles.postpro.areaplot = get (handles.checkbox_areastats,'Value');
handles.postpro.perimeterplot = get(handles.checkbox_perimeter,'Value');
handles.postpro.speedplot = get(handles.checkbox_speed,'Value');
handles.postpro.cellcelldistplot = get(handles.checkbox_cellcelldisthist,'Value');

if (~handles.postpro.cellclusterplot & ~handles.postpro.areaplot & ...
    ~handles.postpro.perimeterplot & ~handles.postpro.speedplot & ...
    ~handles.postpro.cellcelldistplot)
   h = errordlg ('No plots selected. Please select a plot first...');
   uiwait(h);          % Wait until the user presses the OK button
   return;
else
%    if handles.postpro.cellclusterplot
%       % Generate single cell and cluster plots if the users requested these
%       ptPlotCellClusterStats (imageName, savePath, xAxis, cellAmount, clusterAmount, cellsPerCluster, ...
%                               singleCellAmount, percentageSingleCells, percentageClusteredCells);
%    end   
%    
%    if handles.postpro.areaplot
%       % Generate area plots if the users requested these
%       ptPlotAreaStats (imageName, savePath, xAxis, areaPerSingleCell, areaPerCluster);
%    end
% 
%    if handles.postpro.perimeterplot
%       % Generate perimater plots if the users requested these
%       ptPlotPerimeterStats (imageName, savePath, xAxis, perimeterLength, perimeterDivArea);
%    end

   % Here is where the bulk of the graphing work is done; we give it the
   % postpro structure and MPM matrix to work with
   % First do the area and perimeter plots
   if handles.postpro.cellclusterplot | handles.postpro.areaplot | handles.postpro.perimeterplot
      ptPlotCellValues (handles.postpro);
   end
   
   % Then do the cell-cell distance histograms
   if handles.postpro.cellcelldistplot
      ptPlotHistValues (handles.postpro);
   end
   
   % Do the speed plots as well if the user wants it
   if handles.postpro.speedplot
      ptPlotSpeedValues (handles.postpro, handles.MPM);
   end
end

%----------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_fm_movieimgone_ed_CreateFcn (hObject, eventdata, handles)
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

% % --- Executes on button press in GUI_fm_inclbody_rb.
% function GUI_fm_inclbody_rb_Callback(hObject, eventdata, handles)
% % hObject    handle to GUI_fm_inclbody_rb (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of GUI_fm_inclbody_rb

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
% Check that the plot frame values provided are in range
handles = guidata (hObject);

if (handles.postpro.moviefirstimg < handles.postpro.firstimg) | ...
   (handles.postpro.movielastimg > handles.postpro.lastimg)
   h = errordlg ('Movie start and end frame value are out of range. Please re-enter values...');
   uiwait(h);          % Wait until the user presses the OK button
   return;
end

% Start the function that will create the dragtail movie
ptMovieMaker (handles.postpro, handles.MPM);

% Update handles structure
guidata (hObject, handles);

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
handles = guidata (hObject);

% Show an empty figure
viewPrepFigH = figure;

% Show the user a dialog that explains what to do
h = helpdlg ('Please resize the figure and close the figure window. The selected size will be used for the movie.');

% Set the dialog window on a certain screen position
set (h, 'Position', [320.2500 272.2500 297.7500 79.5000]);

% wait for the user to press the ok button after resizing the figure window
uiwait (h);

% Get the new figure size
handles.postpro.figureSize = get (viewPrepFigH, 'Position');

% Close the figure    
close (viewPrepFigH);

% Update handles structure
guidata (hObject, handles);

%----------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_ad_binsize_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_ad_binsize_ed (see GCBO)
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

function GUI_ad_binsize_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_binsize_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_ad_binsize_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_ad_binsize_ed as a double
handles = guidata (hObject);

strval = get (hObject,'String');
val = str2double (strval);
if isnan (val)
   h = errordlg('Sorry, this field has to contain a number.');
   uiwait(h);          % Wait until the user presses the OK button
   handles.postpro.binsize = 4;
   set (handles.GUI_ad_binsize_ed, 'String', handles.postpro.binsize);
   return
else
   if val >= 2   % binsize should be at least 2
      handles.postpro.binsize = val;
   else
      h = errordlg('Sorry, the bin size should be at least 2.');
      uiwait(h);          % Wait until the user presses the OK button
      handles.postpro.binsize = 4;
      set (handles.GUI_ad_binsize_ed, 'String', handles.postpro.binsize);
   end
end

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

% Get the value of 'movie filename'
filename = get (hObject, 'String');

% Assign it to the postpro structure
handles.postpro.dragtailfile = filename;

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------

% --- Executes on button press in checkbox_clustercellstats.
function checkbox_clustercellstats_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_clustercellstats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_clustercellstats

%--------------------------------------------------------------------------

% --- Executes on button press in checkbox_perimeter.
function checkbox_perimeter_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_perimeter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_perimeter


% --- Executes on button press in checkbox_speed.
function checkbox_speed_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_speed

%--------------------------------------------------------------------------

% --- Executes on button press in checkbox_cellcelldisthist.
function checkbox_cellcelldisthist_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_cellcelldisthist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_cellcelldisthist

%--------------------------------------------------------------------------

% --- Executes on button press in checkbox_areastats.
function checkbox_areastats_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_areastats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_areastats

%--------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function multFrameVelocity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to multFrameVelocity (see GCBO)
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

function multFrameVelocity_Callback(hObject, eventdata, handles)
% hObject    handle to multFrameVelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of multFrameVelocity as text
%        str2double(get(hObject,'String')) returns contents of multFrameVelocity as a double
handles = guidata(hObject);

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.postpro.multFrameVelocity = 1;
    set (handles.multFrameVelocity, 'String', handles.postpro.multFrameVelocity);
    return
else
    handles.postpro.multFrameVelocity = val;
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pp_firstframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pp_firstframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function pp_firstframe_Callback(hObject, eventdata, handles)
% hObject    handle to pp_firstframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pp_firstframe as text
%        str2double(get(hObject,'String')) returns contents of pp_firstframe as a double


% --- Executes during object creation, after setting all properties.
function pp_lastframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pp_lastframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function pp_lastframe_Callback(hObject, eventdata, handles)
% hObject    handle to pp_lastframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pp_lastframe as text
%        str2double(get(hObject,'String')) returns contents of pp_lastframe as a double


% --- Executes during object creation, after setting all properties.
function pp_increment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pp_increment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function pp_increment_Callback(hObject, eventdata, handles)
% hObject    handle to pp_increment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pp_increment as text
%        str2double(get(hObject,'String')) returns contents of pp_increment as a double



function GUI_mmpixel_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_mmpixel_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_mmpixel_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_mmpixel_ed as a double
handles = guidata(hObject);

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.postpro.mmpixel = handles.postpro.default.mmpixel;
    set (handles.GUI_mmpixel_ed, 'String', handles.postpro.mmpixel);
    return
else
    handles.postpro.mmpixel = val;
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function GUI_mmpixel_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_mmpixel_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function GUI_frameinterval_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_frameinterval_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_frameinterval_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_frameinterval_ed as a double
handles = guidata(hObject);

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.postpro.timeperframe = handles.postpro.default.timeperframe;
    set (handles.GUI_frameinterval_ed, 'String', handles.postpro.timeperframe);
    return
else
    handles.postpro.timeperframe = val;
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function GUI_frameinterval_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_frameinterval_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function GUI_movietype_avi_rb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_movietype_avi_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on button press in GUI_movietype_avi_rb.
function GUI_movietype_avi_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_movietype_avi_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of GUI_movietype_avi_rb
handles = guidata(hObject);

val = get (hObject,'Value');

if val == 1
   % Button was clicked so should become 1
   handles.postpro.movietype = 1;   % avi
   set (handles.GUI_movietype_avi_rb, 'Value', 1);
   
   % The QT button should automatically become 0
   set (handles.GUI_movietype_qt_rb, 'Value', 0);
else
   % Button was already selected, but make sure that movietype is set
   % correctly
   handles.postpro.movietype = 1;   % avi
end

% Update handles structure
guidata(hObject, handles);   
   
%---------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_movietype_qt_rb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_movietype_qt_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on button press in GUI_movietype_qt_rb.
function GUI_movietype_qt_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_movietype_qt_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of GUI_movietype_qt_rb
handles = guidata(hObject);

val = get (hObject,'Value');

if val == 1
   % Button was clicked so should become 2
   handles.postpro.movietype = 2;   % qt
   set (handles.GUI_movietype_qt_rb, 'Value', 1);
   
   % The QT button should automatically become 0
   set (handles.GUI_movietype_avi_rb, 'Value', 0);
else
   % Button was already selected, but make sure that movietype is set
   % correctly
   handles.postpro.movietype = 2;   % qt
end

% Update handles structure
guidata(hObject, handles);   
