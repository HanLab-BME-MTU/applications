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
% Andre Kerstens        Jul 04          Added area/convex-hull-area plot
% Andre Kerstens        Aug 04          Complete redesign to be able to handle multiple movies

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
defaultPostPro = struct('minusframes', 5, 'plusframes', 2, 'minimaltrack', 5, 'maxsearch', 81, ...
                        'dragtail', 6, 'dragtailfile', 'trackMovie', 'figureSize', [], ...
                        'multframevelocity', 1, 'binsize', 13, 'mmpixel', 0.639, 'timeperframe', 300, ...
                        'movietype', 1, 'nrtrajectories', 5, 'neighbourdist', 81, 'windowsize', 5, ...
                        'maxcellcelldist', 81, 'ripconfint', 85, 'fulltracks', 0, 'dragtracks', 1, ...
                        'drugtimepoint', 30);

% Assign the default postprocessing values to the GUI handle so it can be passed around
handles.defaultPostPro = defaultPostPro;

% Set the HOME env variable if not set already
home = getenv('HOME');
if isempty (home)
   if ispc
      home = 'C:';
   else
      home = '/tmp';
   end
   fprintf (1, 'HOME environment variable not set. Setting default: %s\n', home);
end

% Get environment variable POLYDATA if it exists
polyDataDirectory = getenv('POLYDATA');
if isempty(polyDataDirectory)
    if ispc
        polyDataDirectory = 'C:';
    else  % Unix
        polyDataDirectory = '/tmp';
    end
    fprintf (1, 'POLYDATA environment variable not set. Setting default: %s\n', polyDataDirectory);    
end

% Create results save dir if it doesn't exist yet
if ~exist([polyDataDirectory filesep 'ptData'],'dir')
    try
        mkdir(polyDataDirectory,'ptData');
    catch
        fprintf (1, 'Error: cannot create ptData directory under %s. Please create it manually.\n', polyDataDirectory);
    end
end

% Create movies save dir if it doesn't exist yet
if ~exist([polyDataDirectory filesep 'ptMovies'],'dir')
    try
        mkdir(polyDataDirectory,'ptMovies');
    catch
        fprintf (1, 'Error: cannot create ptMovies directory under %s. Please create it manually.\n', polyDataDirectory);
    end
end

% Update biodata and tmp gui text boxes
set (handles.text_polydatadir_ptpp,'String',polyDataDirectory);

% Assign biodata and tmp to the handles struct
handles.polyDataDirectory = polyDataDirectory;

% Set settings path
handles.guiData.savesettingpath = [polyDataDirectory filesep 'fileInfoPP.mat'];
set (handles.GUI_savesettingpath_ed, 'String', handles.guiData.savesettingpath);

% Set save path
handles.guiData.savedatapath = [polyDataDirectory filesep 'ptData'];
set (handles.GUI_fm_saveallpath_ed, 'String', handles.guiData.savedatapath);

% Set movie file path
handles.guiData.dragtailfile = [polyDataDirectory filesep 'ptMovies' filesep 'trackMovie'];
set (handles.GUI_fm_filename_ed, 'String', handles.guiData.dragtailfile);

% Set binsize
handles.guiData.binsize = defaultPostPro.binsize;
set (handles.GUI_ad_binsize_ed, 'String', handles.guiData.binsize);

% Set relink distance
handles.guiData.relinkdistance = defaultPostPro.maxsearch;
set (handles.GUI_app_relinkdist_ed, 'String', handles.guiData.relinkdistance);

% Set plusframes
handles.guiData.plusframes = defaultPostPro.plusframes;
set (handles.GUI_app_plusframes_ed, 'String', handles.guiData.plusframes);

% Set minusframes
handles.guiData.minusframes = defaultPostPro.minusframes;
set (handles.GUI_app_minusframes_ed, 'String', handles.guiData.minusframes);

% Set trackdistance
handles.guiData.mintrackdistance = defaultPostPro.minimaltrack;
set (handles.GUI_app_minimaltrack_ed, 'String', handles.guiData.mintrackdistance);

% Set dragtail length
handles.guiData.dragtail = defaultPostPro.dragtail;
set (handles.GUI_fm_tracksince_ed, 'String', handles.guiData.dragtail);

% Set the movie type
set (handles.GUI_movietype_avi_rb, 'Value', 0);
set (handles.GUI_movietype_qt_rb, 'Value', 1);

% Set track type
set (handles.GUI_fm_incltracks_rb, 'Value', 1);
set (handles.GUI_fm_fulltracks_cb, 'Value', 0);

% Set window size
handles.guiData.windowsize = defaultPostPro.windowsize;
set (handles.GUI_windowsize_ed, 'String', handles.guiData.windowsize);

% Set multiple frame velocity
handles.guiData.multframevelocity = defaultPostPro.multframevelocity;
set (handles.multFrameVelocity, 'String', handles.guiData.multframevelocity);

% Set neighbourhood traj length
handles.guiData.nrtrajectories = defaultPostPro.nrtrajectories;
set (handles.nr_traj_ed, 'String', handles.guiData.nrtrajectories);

% Set neighbourhood max dist
handles.guiData.maxneighbourdist = defaultPostPro.neighbourdist;
set (handles.neighbour_dist_ed, 'String', handles.guiData.maxneighbourdist);

% Set max cellcell dist
handles.guiData.maxcellcelldist = defaultPostPro.maxcellcelldist;
set (handles.GUI_maxcellcelldist_ed, 'String', handles.guiData.maxcellcelldist);

% Set ripley confidence interval
handles.guiData.ripleyconfint = defaultPostPro.ripconfint;
set (handles.GUI_ripconfint_ed, 'String', handles.guiData.ripleyconfint);

% Set drug application time point
handles.guiData.drugtimepoint = defaultPostPro.drugtimepoint;
set (handles.GUI_drugtimepoint_ed, 'String', handles.guiData.drugtimepoint);

% Set the color of the gui
set(hObject,'Color',[0,0,0.627]);

% Update handles structure
guidata(hObject, handles);

% Turn the resize and int conversion (matlab 7) warnings off
%iptsetpref ('TrueSizeWarning', 'off');

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
function GUI_savesettingpath_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_savesettingpath_ed (see GCBO)
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

function GUI_savesettingpath_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_savesettingpath_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_savesettingpath_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_savesettingpath_ed as a double

handles = guidata(hObject);

% Get the job path and name
saveDir = get(hObject,'String');

if ~exist(saveDir, 'file')
   h=errordlg('The selected path does not exist. Please select another path... ');
   uiwait(h);
   return
end

% Save in the handles struct
handles.savesettingpath = saveDir;

% Update handles structure
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
if exist ('jobValPath', 'dir')
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

% Check that files have been selected before
if ~isfield (handles, 'allMPM')  
    errorStr = ['Jobs should be selected first by using the Select button!'];
    h = errordlg(errorStr);
    uiwait(h);          % Wait until the user presses the OK button  
    return;
end

% Make sure only 1 job is selected
if length(handles.allMPM) > 1
    errorStr = ['Manual processing can only be done for one job at a time. Please select only one job from the list.'];
    h = errordlg(errorStr);
    uiwait(h);          % Wait until the user presses the OK button  
    return;
end

% Let's make sure images are available to make the movie
if handles.jobData(1).imagesavailable == 1
    % This is to signal that the ptManualPostProcessJob function is called by the
    % GUI manual processing button
    handles.whichcallback = 1;

    % Update handles structure
    guidata(hObject, handles);

    % Do the manual postprocessing 
    ptManualPostProcessJob (hObject);
else
    errorStr = ['No images available for this job to do manual postprocessing.'];
    h = errordlg(errorStr);
    uiwait(h);          % Wait until the user presses the OK button  
    return;
end

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
handles.guiData.minusframes = str2num (num);

% Update handles structure
guidata (hObject, handles);

% Update the field on the GUI
set (handles.GUI_app_minusframes_ed, 'String', handles.guiData.minusframes);

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
handles.guiData.plusframes = str2num (num);

% Update handles structure
guidata (hObject, handles);

% Update the field on the GUI
set (handles.GUI_app_plusframes_ed, 'String', handles.guiData.plusframes);

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
handles.guiData.maxdistpostpro = str2num (num);

% Update handles structure
guidata (hObject, handles);

% Update the field on the GUI
set (handles.GUI_app_relinkdist_ed, 'String', handles.guiData.maxdistpostpro);

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
handles.guiData.minimaltrack = str2num (num);

% Update handles structure
guidata (hObject, handles);

% Update the field on the GUI
set (handles.GUI_app_minimaltrack_ed, 'String', handles.guiData.minimaltrack);

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_app_autopostpro_pb.
function GUI_app_autopostpro_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_app_autopostpro_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Check that files have been selected before
if ~isfield (handles, 'allMPM')  
    errorStr = ['Jobs should be selected first by using the Select button!'];
    h = errordlg(errorStr);
    uiwait(h);          % Wait until the user presses the OK button  
    return;
end

% Assign the radiobutton values to the radioButtons struct
radioButtons = getRadiobuttonValues (handles);

% Get the needed GUI data
plusFrames = handles.guiData.plusframes;
minusFrames = handles.guiData.minusframes;
maxDistance = handles.guiData.relinkdistance;
minTrackDistance = handles.guiData.mintrackdistance;
saveDir = handles.guiData.savedatapath;

% Process all the MPMs in the list
for iCount = 1 : length (handles.allMPM)
   updatedMPM = ptTrackFilter(handles.allMPM{iCount}, plusFrames, minusFrames, maxDistance, minTrackDistance, saveDir);
   handles.allMPM{iCount} = updatedMPM;
end

% Show a message telling the user we've finished
msgbox ('Finished auto-processing MPMs. Press OK to continue...');

% Update handles structure
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

% Check that files have been selected before
if ~isfield (handles, 'allMPM')  
    errorStr = ['Jobs should be selected first by using the Select button!'];
    h = errordlg(errorStr);
    uiwait(h);          % Wait until the user presses the OK button  
    return;
end

% Make sure only 1 job is selected
if length(handles.allMPM) > 1
    errorStr = ['Selecting cells can only be done for one job at a time. Please select only one job from the list.'];
    h = errordlg(errorStr);
    uiwait(h);          % Wait until the user presses the OK button  
    return;
end

% Ask the user which frame he wants to see
frameNr = str2double(inputdlg('Provide the frame number:'));
if frameNr == NaN
    errorStr = ['Error: Please enter a numeric value for the frame.'];
    h = errordlg(errorStr);
    uiwait(h);          % Wait until the user presses the OK button  
    return;
end

% Check that it is a valid number
if frameNr < handles.jobData(1).firstimg | frameNr > handles.jobData(1).lastimg
    errorStr = ['Error: Please enter a number between ' num2str(handles.jobData(1).firstimg) ...
                ' and ' num2str(handles.jobData(1).lastimg) '.'];
    h = errordlg(errorStr);
    uiwait(h);          % Wait until the user presses the OK button  
    return;
end

% Let the user select their cells
selectedCells = ptUserCellSelection(frameNr,handles);

% Fill the selected cells matrix in the jobdata handles
% If it's empty we have to initialize it as well
if isempty(handles.jobData(1).selectedcells)
    handles.jobData(1).selectedcells = fillSelectedCellsMatrix (selectedCells, handles.allMPM{1}, frameNr, 1, []);
else
    handles.jobData(1).selectedcells = fillSelectedCellsMatrix (selectedCells, handles.allMPM{1}, frameNr, 0, handles.jobData(1).selectedcells);
end

% Update handles structure
guidata(hObject, handles);

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
handles.guiData.plotfirstimg = str2num (num);

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
handles.guiData.plotlastimg = str2num (num);

% Update handles structure
guidata(hObject, handles);

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_ad_plot_pb.
function GUI_ad_plot_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_ad_plot_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata (hObject);

guiData = ptRetrieveGUIData(handles);

% Check that files have been selected before
if ~isfield (handles, 'allMPM')
    errorStr = ['Jobs should be selected first by using the Select button!'];
    h = errordlg(errorStr);
    uiwait(h);          % Wait until the user presses the OK button  
    return;
end

% Get the path where plot data will be saved
saveDir = guiData.savedatapath;

% If it doesn't exist yet, create it in the results directory
if ~exist (saveDir, 'dir')
   mkdir (saveDir);
end

% Test if the save directory already contains files
[answer, empty] = directoryEmpty (saveDir);

% If it does and the user doesn't want to overwrite, do nothing and return
if strcmp(answer,'') & ~empty
   return
end

% Name used for plots is the current date and time in format YYYYMMDDThhmmss
plotName = datestr(now,30);

% Get window size for running average
windowSize = guiData.windowsize;

% Get drug timepoint
drugTimepoint = guiData.drugtimepoint;

% Assign the radiobutton values to the radioButtons struct
radioButtons = getRadiobuttonValues (handles);

if (~radioButtons.cellclusterplot & ~radioButtons.areaplot & ...
    ~radioButtons.perimeterplot & ~radioButtons.speedplot & ...
    ~radioButtons.cellcelldistplot & ~radioButtons.neighbourplot & ...
    ~radioButtons.ripleyplot)
   h = errordlg ('No plots selected. Please select a plot first...');
   uiwait(h);          % Wait until the user presses the OK button
   return;
else
   
   % These calculations can take a while so set the mouse pointer to busy
   set(gcbf,'Pointer','watch');
      
   % Do the cell/cluster, area and perimeter plots
   if radioButtons.cellclusterplot | radioButtons.areaplot | radioButtons.perimeterplot
                
      % Run one iteration of the calculation
      [cellClusterStats, areaStats, perimeterStats, xAxis] = ptCalculatePlotValues (handles);
               
      % Here's where the plotting itself starts
      if radioButtons.cellclusterplot
         % Generate single cell and cluster plots if the users requested these
         ptPlotCellClusterStats (radioButtons, plotName, saveDir, xAxis, cellClusterStats, windowSize);
      end   
      if radioButtons.areaplot
         % Generate area plots if the users requested these
         ptPlotAreaStats (radioButtons, plotName, saveDir, xAxis, areaStats, windowSize);
      end

      if radioButtons.perimeterplot
         % Generate perimater plots if the users requested these
         ptPlotPerimeterStats (radioButtons, plotName, saveDir, xAxis, perimeterStats, windowSize);
      end

      % For all the figures we want to keep the xAxis as well 
      save ([saveDir filesep plotName '_xAxis-CellStats.mat'],'xAxis');
   end
   
   % Then do the cell-cell distance plots
   if radioButtons.cellcelldistplot
       
      % Run one iteration of the calculation
      [cellCellDistStats, xAxis] = ptCalculateCellCellDist (handles);
               
      % Create the plots
      ptPlotCellCellDist (radioButtons, plotName, saveDir, xAxis, cellCellDistStats, windowSize);
   end
   
   % Initialize some tmp vars
   xAxisLengthPrev = 0;
  
   % Do the speed plots if the user requested these
   if radioButtons.speedplot
       
      % Run the calculation for the velocity stats
      [avgVelocityStats, velocitySingleStats, velocityVarStats, velocityHistStats, xAxis] = ...
                               ptCalculateSpeedValues (handles);
      
      % Here's where the plotting itself starts
      if radioButtons.speedplot_2
         % Generate avg velocity plots if the users requested these
          ptPlotSpeedStats (radioButtons, plotName, saveDir, xAxis, avgVelocityStats, windowSize, ...
                            drugTimepoint);
      end   
      if radioButtons.speedplot_1
         % Generate vel. single cell plots if the users requested these
          ptPlotSingleSpeedStats (radioButtons, plotName, saveDir, xAxis, velocitySingleStats, windowSize);
      end

      if radioButtons.speedplot_3
         % Generate velocity variance plots if the users requested these
          ptPlotSpeedVarStats (radioButtons, plotName, saveDir, xAxis, velocityVarStats, windowSize);
      end

      if radioButtons.speedplot_4 | radioButtons.allcellshist | radioButtons.singlecellshist | ...
         radioButtons.clusteredcellshist
         % Generate velocity histogram plots if the users requested these
          ptPlotVelocityHist (radioButtons, plotName, saveDir, xAxis, velocityHistStats);
      end      
      
      % For all the figures we want to keep the xAxis as well 
      save ([saveDir filesep plotName '_xAxis-Velocity.mat'],'xAxis');
   end

   % Do the neighbourhood plots if the user requested these
   if radioButtons.neighbourplot            
      if radioButtons.neighbourplot_1
                    
         % Run one iteration of the calculation
         [neighTrajStats, xAxis] = ptCalculateNeighbourTraj (handles); 
            
         % Do the plots   
         ptPlotNeighbourTraj (radioButtons, plotName, saveDir, xAxis, neighTrajStats, windowSize, ...
                              drugTimepoint);
      end
            
      if radioButtons.neighbourplot_2
          
         % Run one iteration of the calculation
         [neighChangeStats, xAxis] = ptCalculateNeighbourChanges (handles); 
            
         % Do the plots   
         ptPlotNeighbourChanges (radioButtons, plotName, saveDir, xAxis, neighChangeStats, windowSize, ...
                                 drugTimepoint);
      end
      
      % For all the figures we want to keep the xAxis as well 
      save ([saveDir filesep plotName '_xAxis-Neighbours.mat'],'xAxis');
   end
   
   % Only do chaos stats calculations when image (and indirectly image
   % size) are available
   if radioButtons.ripleyplot         
      if radioButtons.ripleyplot_1
          
          % Run the calculation
          radioButtons = getRadiobuttonValues (handles);
          [chaosStats, xAxis] = ptCalculateChaosStats (handles, radioButtons);
          
          % Do the plots
          ptPlotChaosStats (radioButtons, plotName, saveDir, xAxis, chaosStats, windowSize, ...
                            drugTimepoint);
      end
      
      % For these figures we want to keep the xAxis as well 
      save ([saveDir filesep plotName '_xAxis-Chaos.mat'],'xAxis');
   end
   
   % Set the mouse pointer to normal again
   set(gcbf,'Pointer','arrow');
   
   % Show a message telling the user we've finished
   msgbox ('Finished generating plots and histograms. Press OK to continue...');
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

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get (hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    if isfield(handles.jobData(1),'moviefirstimg')
        handles.guiData.moviefirstimg = handles.jobData(1).firstimg;
        set (handles.GUI_fm_movieimgone_ed, 'String', num2str(handles.guiData.moviefirstimg));  % Revert the value back
    else
        handles.guiData.moviefirstimg = 1;
        set (handles.GUI_fm_movieimgone_ed, 'String', num2str(handles.guiData.moviefirstimg));  % Revert the value back
    end
    return
else
    handles.guiData.moviefirstimg = val;
end

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

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get (hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    if isfield(handles.jobData(1),'movielastimg')
        handles.guiData.movielastimg = handles.jobData(1).lastimg;
        set (handles.GUI_fm_movieimgone_ed, 'String', num2str(handles.guiData.movielastimg));  % Revert the value back
    else
        handles.guiData.movielastimg = 1;
        set (handles.GUI_fm_movieimgend_ed, 'String', num2str(handles.guiData.movielastimg));  % Revert the value back
    end
    return
else
    handles.guiData.movielastimg = val;
end

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------

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
handles = guidata(hObject);

val = get (hObject,'Value');

if val == 1
   % Button was clicked so should become 1
   handles.guiData.dragtracks = 1;
   set (handles.GUI_fm_incltracks_rb, 'Value', 1);
   
   % The QT button should automatically become 0
   set (handles.GUI_fm_fulltracks_cb, 'Value', 0);
   handles.guiData.fulltracks = 0;
else
   % Button was already selected, but make sure that movietype is set
   % correctly
   %handles.guiData.fulltracks = 0;
   handles.guiData.dragtracks = 0;
end

% Update handles structure
guidata(hObject, handles); 

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

strval = get (hObject,'String');
val = str2double (strval);
if isnan (val)
   h = errordlg('Sorry, this field has to contain a number.');
   uiwait(h);          % Wait until the user presses the OK button
   handles.guiData.dragtail = 6;
   set (handles.GUI_fm_tracksince_ed, 'String', handles.guiData.dragtail);
   return
else
   handles.guiData.dragtail = val;
end

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

% If the path doesn't exist yet create it
if ~exist (path, 'dir')
   msgStr = ['This directory does not exist yet. Do you want to create it?'];
   answer = questdlg(msgStr, 'Create Directory', 'Yes', 'No', 'Yes');
   if strcmp(answer,'Yes')
      mkdir (path);
   else
      path = handles.guiData.savedatapath;
      set(handles.GUI_fm_saveallpath_ed, 'String', path);
   end
end

% Assign it to the handles structure
handles.saveallpath = path;

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

% Check that files have been selected before
if ~isfield (handles, 'allMPM')  
    errorStr = ['Jobs should be selected first by using the Select button!'];
    h = errordlg(errorStr);
    uiwait(h);          % Wait until the user presses the OK button  
    return;
end

% Make sure only 1 job is selected
if length(handles.allMPM) > 1
    errorStr = ['A movie can only be generated for one job at a time. Please select only one job from the list.'];
    h = errordlg(errorStr);
    uiwait(h);          % Wait until the user presses the OK button  
    return;
end

% Assign the radiobutton values to the radioButtons struct
radioButtons = getRadiobuttonValues (handles);

% Let's make sure images are available to make the movie
if handles.jobData(1).imagesavailable == 1

    % Start the function that will create the dragtail movie
    result = ptMovieMaker (radioButtons, handles);

    % Show a message telling the user we've finished
    if result == 0
       msgbox ('Finished generating movie. Press OK to continue...');
    end
else
    errorStr = ['No images available for this job to create movie.'];
    h = errordlg(errorStr);
    uiwait(h);          % Wait until the user presses the OK button  
    return;
end

% Update handles structure
guidata (hObject, handles);

%----------------------------------------------------------------------------

% --- Executes on button press in GUI_fm_browsesaveallpath_pb.
function GUI_fm_browsesaveallpath_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_fm_browsesaveallpath_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Let the user browse for an image file and path
saveDirectory = uigetdir('','Select Save Directory');

% And store this in the handles struct
if exist(saveDirectory, 'dir') == 7
   handles.saveallpath = saveDirectory;
   set (handles.GUI_fm_saveallpath_ed, 'String', handles.saveallpath);
else
   h = errordlg('This save directory does not exist. Please select another directory...');
   uiwait(h);          % Wait until the user presses the OK button
   return;
end

% Update handles structure
guidata(hObject, handles);

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
handles.guiData.figureSize = get (viewPrepFigH, 'Position');

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
   handles.guiData.binsize = 4;
   set (handles.GUI_ad_binsize_ed, 'String', handles.guiData.binsize);
   return
else
   if val >= 2   % binsize should be at least 2
      handles.guiData.binsize = val;
   else
      h = errordlg('Sorry, the bin size should be at least 2.');
      uiwait(h);          % Wait until the user presses the OK button
      handles.guiData.binsize = 4;
      set (handles.GUI_ad_binsize_ed, 'String', handles.guiData.binsize);
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
movieDirectory = get (hObject, 'String');
[pathString, filename, ext, version] = fileparts (movieDirectory);

% If the path doesn't exist yet create it
if ~exist (pathString, 'dir')
   msgStr = ['This directory does not exist yet. Do you want to create it?'];
   answer = questdlg(msgStr, 'Create Directory', 'Yes', 'No', 'Yes');
   if strcmp(answer,'Yes')
      mkdir (pathString);
   else
      moviePath = handles.guiData.dragtailfile;
      set(handles.GUI_fm_filename_ed, 'String', moviePath);
   end
end

% Assign it to the guiData structure
handles.guiData.dragtailfile = [pathString filesep filename];

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------

% --- Executes on button press in checkbox_clustercellstats.
function checkbox_clustercellstats_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_clustercellstats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_clustercellstats
handles = guidata(hObject);

val = get (hObject,'Value');

if val == 1
   % Checkbox is selected, so select all children
   set (handles.checkbox_amount_cells, 'Value', 1);
   set (handles.checkbox_percentage_cells, 'Value', 1);
else  % val == 0
   % Checkbox was unselected so unselect all the children
   set (handles.checkbox_amount_cells, 'Value', 0);
   set (handles.checkbox_percentage_cells, 'Value', 0);
end

% Update handles structure
guidata(hObject, handles);   

%--------------------------------------------------------------------------

% --- Executes on button press in checkbox_perimeter.
function checkbox_perimeter_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_perimeter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_perimeter
% handles = guidata(hObject);
% 
% val = get (hObject,'Value');
% 
% if val == 1
%    % Checkbox is selected, so select all children
%    set (handles.checkbox_amount_cells, 'Value', 1);
%    set (handles.checkbox_percentage_cells, 'Value', 1);
% else  % val == 0
%    % Checkbox was unselected so unselect all the children
%    set (handles.checkbox_amount_cells, 'Value', 0);
%    set (handles.checkbox_percentage_cells, 'Value', 0);
% end
% 
% % Update handles structure
% guidata(hObject, handles); 

%--------------------------------------------------------------------------

% --- Executes on button press in checkbox_speed.
function checkbox_speed_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_speed
handles = guidata(hObject);

val = get (hObject,'Value');

if val == 1
   % Checkbox is selected, so select all children
   set (handles.checkbox_all_to_single_speed, 'Value', 1);
   set (handles.checkbox_average_speed, 'Value', 1);
   set (handles.checkbox_speed_variance, 'Value', 1);
   set (handles.checkbox_speed_histogram, 'Value', 1);
else  % val == 0
   % Checkbox was unselected so unselect all the children
   set (handles.checkbox_all_to_single_speed, 'Value', 0);
   set (handles.checkbox_average_speed, 'Value', 0);
   set (handles.checkbox_speed_variance, 'Value', 0);   
   set (handles.checkbox_speed_histogram, 'Value', 0);
end

% Update handles structure
guidata(hObject, handles); 

%--------------------------------------------------------------------------

% --- Executes on button press in checkbox_cellcelldisthist.
function checkbox_cellcelldisthist_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_cellcelldisthist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_cellcelldisthist
handles = guidata(hObject);

val = get (hObject,'Value');

if val == 1
   % Checkbox is selected, so select all children
   set (handles.checkbox_avg_distance_cells, 'Value', 1);
else  % val == 0
   % Checkbox was unselected so unselect all the children
   set (handles.checkbox_avg_distance_cells, 'Value', 0);
end

% Update handles structure
guidata(hObject, handles); 

%--------------------------------------------------------------------------

% --- Executes on button press in checkbox_areastats.
function checkbox_areastats_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_areastats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_areastats
handles = guidata(hObject);

val = get (hObject,'Value');

if val == 1
   % Checkbox is selected, so select all children
   set (handles.checkbox_total_area, 'Value', 1);
   set (handles.checkbox_avg_convex_hull_area, 'Value', 1);
   set (handles.checkbox_single_cluster_area, 'Value', 1);
else  % val == 0
   % Checkbox was unselected so unselect all the children
   set (handles.checkbox_total_area, 'Value', 0);
   set (handles.checkbox_avg_convex_hull_area, 'Value', 0);
   set (handles.checkbox_single_cluster_area, 'Value', 0);
end

% Update handles structure
guidata(hObject, handles); 

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
    handles.guiData.multframevelocity = 1;
    set (handles.multframevelocity, 'String', handles.guiData.multframevelocity);
    return
else
    handles.guiData.multframevelocity = val;
end

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------

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

%--------------------------------------------------------------------------

function pp_firstframe_Callback(hObject, eventdata, handles)
% hObject    handle to pp_firstframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pp_firstframe as text
%        str2double(get(hObject,'String')) returns contents of pp_firstframe as a double

%--------------------------------------------------------------------------

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

%--------------------------------------------------------------------------

function pp_lastframe_Callback(hObject, eventdata, handles)
% hObject    handle to pp_lastframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pp_lastframe as text
%        str2double(get(hObject,'String')) returns contents of pp_lastframe as a double

%--------------------------------------------------------------------------

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

%--------------------------------------------------------------------------

function pp_increment_Callback(hObject, eventdata, handles)
% hObject    handle to pp_increment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pp_increment as text
%        str2double(get(hObject,'String')) returns contents of pp_increment as a double

%--------------------------------------------------------------------------

function GUI_mmpixel_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_mmpixel_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_mmpixel_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_mmpixel_ed as a double
handles = guidata(hObject);

% Get the selected jobs, since we want to change this value for all of them
filesSelected = get(handles.GUI_filelist_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.guiData.mmpixel = handles.defaultPostpro.mmpixel;
    set (handles.GUI_mmpixel_ed, 'String', handles.guiData.mmpixel);
    return
else
    for iCount = 1:length(filesSelected)
       handles.jobData(iCount).mmpixel = val;
    end
end

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------

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

%--------------------------------------------------------------------------

function GUI_frameinterval_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_frameinterval_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_frameinterval_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_frameinterval_ed as a double
handles = guidata(hObject);

% Get the selected jobs, since we want to change this value for all of them
filesSelected = get(handles.GUI_filelist_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.guiData.timeperframe = handles.defaultPostpro.timeperframe;
    set (handles.GUI_frameinterval_ed, 'String', handles.guiData.timeperframe);
    return
else
    for iCount = 1:length(filesSelected)
       handles.jobData(iCount).timeperframe = val;
    end
end

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------

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

%--------------------------------------------------------------------------

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

%--------------------------------------------------------------------------

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
   handles.guiData.movietype = 1;   % avi
   set (handles.GUI_movietype_avi_rb, 'Value', 1);
   
   % The QT button should automatically become 0
   set (handles.GUI_movietype_qt_rb, 'Value', 0);
else
   % Button was already selected, but make sure that movietype is set
   % correctly
   handles.guiData.movietype = 1;   % avi
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

%---------------------------------------------------------------------

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
   handles.guiData.movietype = 2;   % qt
   set (handles.GUI_movietype_qt_rb, 'Value', 1);
   
   % The QT button should automatically become 0
   set (handles.GUI_movietype_avi_rb, 'Value', 0);
else
   % Button was already selected, but make sure that movietype is set
   % correctly
   handles.guiData.movietype = 2;   % qt
end

% Update handles structure
guidata(hObject, handles);   


%---------------------------------------------------------------------

% --- Executes on button press in checkbox_speed_variance.
function checkbox_speed_variance_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_speed_variance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_speed_variance

%--------------------------------------------------------------------------

% --- Executes on button press in checkbox_average_speed.
function checkbox_average_speed_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_average_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_average_speed

%---------------------------------------------------------------------

% --- Executes on button press in checkbox_all_to_single_speed.
function checkbox_all_to_single_speed_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_all_to_single_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_all_to_single_speed

%---------------------------------------------------------------------

% --- Executes on button press in checkbox_avg_distance_cells.
function checkbox_avg_distance_cells_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_avg_distance_cells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_avg_distance_cells

%---------------------------------------------------------------------

% --- Executes on button press in checkbox_avg_convex_hull_area.
function checkbox_avg_convex_hull_area_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_avg_convex_hull_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_avg_convex_hull_area

%---------------------------------------------------------------------

% --- Executes on button press in checkbox_total_area.
function checkbox_total_area_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_total_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_total_area

%---------------------------------------------------------------------

% --- Executes on button press in checkbox_single_cluster_area.
function checkbox_single_cluster_area_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_single_cluster_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_single_cluster_area

%---------------------------------------------------------------------

% --- Executes on button press in checkbox_amount_cells.
function checkbox_amount_cells_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_amount_cells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_amount_cells

%---------------------------------------------------------------------

% --- Executes on button press in checkbox_percentage_cells.
function checkbox_percentage_cells_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_percentage_cells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_percentage_cells

%---------------------------------------------------------------------

% --- Executes on button press in checkbox_neighbourhood.
function checkbox_neighbourhood_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_neighbourhood (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_neighbourhood
handles = guidata(hObject);

val = get (hObject,'Value');

if val == 1
   % Checkbox is selected, so select all children
   set (handles.checkbox_nb_trajectories, 'Value', 1);
   set (handles.checkbox_nb_interact, 'Value', 1);
else  % val == 0
   % Checkbox was unselected so unselect all the children
   set (handles.checkbox_nb_trajectories, 'Value', 0);
   set (handles.checkbox_nb_interact, 'Value', 0);
end

% Update handles structure
guidata(hObject, handles); 

%---------------------------------------------------------------------

% --- Executes on button press in checkbox_nb_trajectories.
function checkbox_nb_trajectories_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_nb_trajectories (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_nb_trajectories

%---------------------------------------------------------------------

function nr_traj_ed_Callback(hObject, eventdata, handles)
% hObject    handle to nr_traj_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nr_traj_ed as text
%        str2double(get(hObject,'String')) returns contents of nr_traj_ed as a double
handles = guidata (hObject);

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get (hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.guiData.nrtrajectories = 5;
    set (handles.nr_traj_ed, 'Value', handles.guiData.nrtrajectories);  % Revert the value back
    return
else
    handles.guiData.nrtrajectories = val;
end

% Update handles structure
guidata(hObject, handles);

%---------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function nr_traj_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nr_traj_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%---------------------------------------------------------------------

% --- Executes on button press in checkbox_nb_interact.
function checkbox_nb_interact_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_nb_interact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_nb_interact

%---------------------------------------------------------------------

function neighbour_dist_ed_Callback(hObject, eventdata, handles)
% hObject    handle to neighbour_dist_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of neighbour_dist_ed as text
%        str2double(get(hObject,'String')) returns contents of neighbour_dist_ed as a double
handles = guidata (hObject);

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get (hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.guiData.neighbourdist = 81;
    set (handles.neighbour_dist_ed, 'Value', handles.guiData.neighbourdist);  % Revert the value back
    return
else
    handles.guiData.neighbourdist = val;
end

% Update handles structure
guidata(hObject, handles);

%---------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function neighbour_dist_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to neighbour_dist_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%--------------------------------------------------------------------

% --- Executes on button press in GUI_vel_all_cells_cb.
function GUI_vel_all_cells_cb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_vel_all_cells_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_vel_all_cells_cb

%--------------------------------------------------------------------

% --- Executes on button press in GUI_vel_single_cells_cb.
function GUI_vel_single_cells_cb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_vel_single_cells_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_vel_single_cells_cb

%--------------------------------------------------------------------

% --- Executes on button press in GUI_vel_clust_cells_cb.
function GUI_vel_clust_cells_cb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_vel_clust_cells_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_vel_clust_cells_cb

% --------------------------------------------------------------------

function GUI_average_data_menu_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_average_data_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata (hObject);

% Call Averaging program
ptAverageData;

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------

function GUI_movie_menu_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_movie_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------

function GUI_make_hist_movie_menu_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_make_hist_movie_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata (hObject);

% Call Averaging program
ptMakeHistogramMovies;

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------

% --- Executes on selection change in GUI_filelist_lb.
function GUI_filelist_lb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_filelist_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GUI_filelist_lb contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GUI_filelist_lb
handles = guidata(hObject);

% Get the job list and the current job
fileList = get(handles.GUI_filelist_lb,'String');
filesSelected = get(handles.GUI_filelist_lb,'Value');

% GUI values can only be updated for 1 job at the time, so in case we have
% more, do nothing
if length(filesSelected) > 1
    return;
end

% Get the values from the job
[allMPM, allCellProps, allClusterProps, allFrameProps, allValidFrames, jobData, result] = ptRetrieveJobData (fileList, 'all');

% Check the result value (0 is good)
if result > 0
   h=errordlg ('An error occured while fetching data for the selected files (ptRetrieveJobData).');
   uiwait (h); 
   return;
end

% Assign to handles struct
handles.jobData = jobData;

% Get radioButton values
alwaysCountFrom1 = get (handles.GUI_alwayscount1_cb, 'Value');

% Default plot and movie ranges will be made similar to the frame range
% (which can be found in validFrames: 1st and last entry)
if ~alwaysCountFrom1
   %handles.guiData.plotfirstimg = jobData(filesSelected).firstimg;
   %handles.guiData.plotlastimg = jobData(filesSelected).lastimg;
   handles.guiData.plotfirstimg = allValidFrames{filesSelected}(1,1);
   handles.guiData.plotlastimg = allValidFrames{filesSelected}(1,end);
else
   handles.guiData.plotfirstimg = 1;
   %handles.guiData.plotlastimg = jobData(filesSelected).lastimg - jobData(filesSelected).firstimg + 1;
   %handles.guiData.plotlastimg = length (allValidFrames{filesSelected}(1,:));
   handles.guiData.plotlastimg = allValidFrames{filesSelected}(1,end) - allValidFrames{filesSelected}(1,1) + 1;
end

% Set movie frame start and end
%handles.guiData.moviefirstimg = jobData(filesSelected).firstimg;
%handles.guiData.movielastimg = jobData(filesSelected).lastimg;
handles.guiData.moviefirstimg = allValidFrames{filesSelected}(1,1);
handles.guiData.movielastimg = allValidFrames{filesSelected}(1,end);
   
% Set values on the GUI
handles = ptSetPostproGUIValues (handles, filesSelected);

% Update GUI handles struct
guidata (hObject,handles);

%--------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_filelist_lb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_filelist_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%--------------------------------------------------------------------

% --- Executes on button press in GUI_add_pb.
function GUI_add_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_add_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Start at the polytrack directory (POLYDATA)
currentDir = pwd;
cd (handles.polyDataDirectory);

% Select an image filename or a file called 'jobvalues.mat' from a user selected directory
[filename, directory] = uigetfile ({'MPM.mat', 'MPM-Files'; '*.*', 'all files'}, ...
                                      'Please select an MPM file');

% Do nothing in case the user doesn't select anything
if filename == 0
   cd (currentDir);
   return;
else
   % Convert filename to lowercase
   fileLower = lower (filename);
   
   % Check whether it is a MPM file (ext .mat)
   if ~strcmp(fileLower, 'mpm.mat')
      errormsg = ['File ' filename ' is not a MPM file. Please choose a file named MPM.mat.'];
      h = errordlg (errormsg);
      uiwait (h);
      return
   end
   
   % Get the current file list
   fileList = get (handles.GUI_filelist_lb, 'String');
   
   % Cat together the directory and filename
   filePath = [directory filename];
      
   % If fileList consists of more than one entry, append the new file path after the
   % last one else just put it in the fileList as the first entry
   if ~iscell(fileList)
      fileList = cell(1);
      fileList(1) = cellstr(filePath);
   else   % The list already had some files in it 
       
      % Add the filename to the list
      fileList(end+1) = cellstr(filePath);
   end
       
   % Show the values for this job on the GUI
   filePath = fileList{end};

   % Get the values from the job
   [allMPM, allCellProps, allClusterProps, allFrameProps, allValidFrames, jobData, result] = ptRetrieveJobData (fileList, 'all');

   % Check the result value (0 is good)
   if result > 0
      h=errordlg ('An error occured while fetching data for the selected files (ptRetrieveJobData).');
      uiwait (h); 
      return;
   end
   
   % Get the values from the GUI
   [guiData] = ptRetrieveGUIData (handles);
   
   % Get radioButton values
   alwaysCountFrom1 = get (handles.GUI_alwayscount1_cb, 'Value');

   % Default plot and movie ranges will be made similar to the frame range
   % (which can be found in validFrames: 1st and last entry)
   if ~alwaysCountFrom1
      %guiData.plotfirstimg = jobData(filesSelected).firstimg;
      %guiData.plotlastimg = jobData(filesSelected).lastimg;
      guiData.plotfirstimg = allValidFrames{end}(1,1);
      guiData.plotlastimg = allValidFrames{end}(1,end);
   else
      guiData.plotfirstimg = 1;
      %guiData.plotlastimg = jobData(end).lastimg - jobData(end).firstimg + 1;
      guiData.plotlastimg = length (allValidFrames{end}(1,:));
   end

   % Set movie frame start and end
   %guiData.moviefirstimg = jobData(end).firstimg;
   %guiData.movielastimg = jobData(end).lastimg;
   guiData.moviefirstimg = allValidFrames{end}(1,1);
   guiData.movielastimg = allValidFrames{end}(1,end);

   % Assign all values to the handles struct
   handles.allMPM = allMPM;
   handles.allCellProps = allCellProps;
   handles.allClusterProps = allClusterProps;
   handles.allFrameProps = allFrameProps;
   handles.allValidFrames = allValidFrames;
   handles.jobData = jobData;
   handles.guiData = guiData;

   % Set values on the GUI
   handles = ptSetPostproGUIValues (handles, length(fileList));

   % Go to the selected directory: user comfort for next file
   cd (directory);
end
      
% Store the modified job list back in the GUI handle
set(handles.GUI_filelist_lb,'String',fileList);
set (handles.GUI_filelist_lb, 'Value', length(fileList));

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------

% --- Executes on button press in GUI_remove_pb.
function GUI_remove_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_remove_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Get the job list and the current job
fileList = get(handles.GUI_filelist_lb,'String');
fileNumber = get(handles.GUI_filelist_lb,'Value');

% fileList will only then be a cell, if there are results in it.
% Otherwise it is a string (No csv's loaded)
if ~iscell(fileList)
    % Show an error dialog with an appropriate message and wait
    % for the user to press a button
    h=errordlg('Sorry, there are no files to remove.');
    uiwait(h);
    return;
end

% Remove files from the list (multiple files are possible) and empty the
% data from the cell structs
if length(fileNumber) > 0
   fileList(fileNumber) = [];
   
%    handles.jobData(fileNumber) = [];
%    handles.allMPM(fileNumber) = [];
%    handles.allCellProps(fileNumber) = [];
%    handles.allClusterProps(fileNumber) = [];
%    handles.allFrameProps(fileNumber) = [];
end

% Put a standard string in the file list window if there are no files to show
% and delete the currently selected file from the gui
if length(fileList) == 0
   fileList = char('No files loaded.');
end

% Set the list to the first file to be on the safe side
% And store new jobList in gui handle
set (handles.GUI_filelist_lb, 'Value', 1);
set (handles.GUI_filelist_lb, 'String', fileList);

% Update GUI handles struct
guidata (hObject,handles);

%--------------------------------------------------------------------

% --- Executes on button press in GUI_load_pb.
function GUI_load_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_load_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Retrieve the directory and filename where to save the result
saveDirectory = get(handles.GUI_savesettingpath_ed,'String');
[pathString, filename, ext, version] = fileparts (saveDirectory);

if ~exist (pathString, 'dir')
   % Show an error dialog with an appropriate message and wait
   % for the user to press a button
   h=errordlg ('The save directory does not exist. Please choose another directory.');
   uiwait (h);
   return
end

% Select an setting file from a user selected directory
[filename, directory] = uigetfile ({[handles.polyDataDirectory filesep '*.mat'], 'Setting Files'; '*.*', 'all files'}, ...
                                      'Please select a PP setting file');

% Do nothing in case the user doesn't select anything
if filename == 0
   cd (currentDir);
   return
else
    
   % Load previous list from disk
   load ([directory filename]);

   % Check that a fileInfoPP struct exist in the variable space
   if ~exist ('fileInfoPP', 'var')
      h=errordlg('This is not a valid setting file. Please select another one.');
      uiwait(h);
      return;
   else
       
      % Get file info from struct
      fileList = fileInfoPP.fileList;
      filesSelected = fileInfoPP.filesSelected;
      saveDir = fileInfoPP.saveDir;

      % Set the GUI values for the first one
      filePath = fileList{1};

      if ~exist(filePath)
         h=errordlg('These job directories do not exist. Please load another set.');
         uiwait(h);
         return;
      end
      
      % Get the values from the job
      [allMPM, allCellProps, allClusterProps, allFrameProps, allValidFrames, jobData, result] = ptRetrieveJobData (fileList, 'all');

      % Check the result value (0 is good)
      if result == 1
         return;
      end

      % Get the values from the GUI
      [guiData] = ptRetrieveGUIData (handles);
   
      % Get radioButton values
      alwaysCountFrom1 = get (handles.GUI_alwayscount1_cb, 'Value');

      % Default plot and movie ranges will be made similar to the frame range
      if ~alwaysCountFrom1
         handles.guiData.plotfirstimg = jobData(filesSelected).firstimg;
         handles.guiData.plotlastimg = jobData(filesSelected).lastimg;
      else
         handles.guiData.plotfirstimg = 1;
         %handles.guiData.plotlastimg = jobData(filesSelected).lastimg - jobData(filesSelected).firstimg + 1;
         handles.guiData.plotlastimg = allValidFrames{filesSelected(end)}(1,end) - allValidFrames{filesSelected(end)}(1,1) + 1;
      end

      % Set movie frame start and end
      handles.guiData.moviefirstimg = jobData(filesSelected(1)).firstimg;
      handles.guiData.movielastimg = jobData(filesSelected(1)).lastimg;
      
      % Assign all values to the handles struct
      handles.allMPM = allMPM;
      handles.allCellProps = allCellProps;
      handles.allClusterProps = allClusterProps;
      handles.allFrameProps = allFrameProps;
      handles.allValidFrames = allValidFrames;
      handles.jobData = jobData;
      handles.guiData = guiData;

      % Update the GUI with the job list and the current job
      set (handles.GUI_filelist_lb, 'String', fileList);
      set (handles.GUI_filelist_lb, 'Value', filesSelected);
      
      % Set values on the GUI
      handles = ptSetPostproGUIValues (handles, length(fileList));
      
      % Set the savepath as well
      set (handles.GUI_fm_saveallpath_ed, 'String', saveDir);
   end
end

% Update GUI handles struct
guidata (hObject,handles);

%--------------------------------------------------------------------

% --- Executes on button press in GUI_save_pb.
function GUI_save_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_save_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Retrieve the directory and filename where to save the result
saveDirectory = get(handles.GUI_savesettingpath_ed,'String');
[pathString, filename, ext, version] = fileparts (saveDirectory);

if ~exist (pathString, 'dir')
   % Show an error dialog with an appropriate message and wait
   % for the user to press a button
   h=errordlg ('The save directory does not exist. Please choose another directory.');
   uiwait (h);
   return
end

% Get the job list and the current job
fileList = get(handles.GUI_filelist_lb,'String');
filesSelected = get(handles.GUI_filelist_lb,'Value');
saveDir = get(handles.GUI_fm_saveallpath_ed,'String');

% Store this info in a struct so that we can save it easily
fileInfoPP.fileList = fileList;
fileInfoPP.filesSelected = filesSelected;
fileInfoPP.saveDir = saveDir;

% Ask the user where to save the file
[filename,path] = uiputfile(saveDirectory, 'Save settings as');

% If user presses cancel, don't do anything
if filename ~= 0
    % Save to disk
    save ([path filename], 'fileInfoPP');

    % Modify the GUI and handles as well
    handles.savepath = [path filename];
    set (handles.GUI_savesettingpath_ed, 'String', handles.savepath);
end

% Update GUI handles struct
guidata (hObject,handles);

%--------------------------------------------------------------------

% --- Executes on button press in GUI_moviebrowse_pb.
function GUI_moviebrowse_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_moviebrowse_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Retrieve the directory and filename where to save the result
movieDirectory = get(handles.GUI_fm_filename_ed,'String');
[pathString, filename, ext, version] = fileparts (movieDirectory);

% Select an image filename or a file called 'jobvalues.mat' from a user selected directory
%[filename, directory] = uigetfile ({movieDirectory, 'movie files'; '*.*', 'all files'}, ...
%                        'Please select a filename');
directory = uigetdir (pathString, 'Please select a directory');

% Do nothing in case the user doesn't select anything
if directory == 0
   return
else   
   % Get the job list and the current job
   set (handles.GUI_fm_filename_ed, 'String', [directory filesep filename]);
end

% Update GUI handles struct
guidata (hObject,handles);

%--------------------------------------------------------------------

function checkbox_speed_histogram_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_moviebrowse_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

%--------------------------------------------------------------------

function [mpmNr, mpmLength] = max_MPM_length (allMPM)
% Returns the length of the longest MPM in the list of MPM's

% Initialize vars
prevLength = 0;
mpmNr = 0;
mpmLength = 0;
      
% Go throught the list of MPMs
for iCount = 1 : length(allMPM)
   
   % Test for length and keep if longer
   curLength = size(allMPM{iCount},2);
   if curLength > prevLength
      prevLength = curLength;
  
      mpmNr = iCount;
      mpmLength = curLength;
   end
end

%--------------------------------------------------------------------

% --- Executes on button press in GUI_select_pb.
function GUI_select_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_select_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Get the job list and the current job
fileList = get(handles.GUI_filelist_lb,'String');
filesSelected = get(handles.GUI_filelist_lb,'Value');

if ~strcmp(fileList,'No files loaded.')
    % Retrieve the data for the selected jobs
    [allMPM, allCellProps, allClusterProps, allFrameProps, allValidFrames, jobData, result] = ptRetrieveJobData (fileList, filesSelected);

    % Check the result value (0 is good)
    if result > 0
       h=errordlg ('An error occured while fetching data for the selected files (ptRetrieveJobData).');
       uiwait (h); 
       return;
    end

    % Get the data from the GUI
    [guiData] = ptRetrieveGUIData (handles);

    % Get radioButton values
    alwaysCountFrom1 = get (handles.GUI_alwayscount1_cb, 'Value');

    % Check that the job and gui data fit together over movies
    for jobCount = 1 : length(allMPM)

        % Get first and last frame numbers and increment
        startFrame = jobData(jobCount).firstimg;
        endFrame = jobData(jobCount).lastimg;
        increment = jobData(jobCount).increment;

        % Default plot ranges should be tested if not started from 1
        if ~alwaysCountFrom1
            % Make sure the start and end frames fit to the selected plot frames
            if (guiData.plotfirstimg < startFrame)
                errorStr = ['The selected plot start frame (' num2str(guiData.plotfirstimg) ') does not fit with the job start frame (' num2str(startFrame) ')'];
                h = errordlg(errorStr);
                uiwait(h);          % Wait until the user presses the OK button  
                return;
            end
            if (guiData.plotlastimg > endFrame)
                errorStr = ['The selected plot end frame (' num2str(guiData.plotlastimg) ') does not fit with the job end frame (' num2str(endFrame) ')'];
                h = errordlg(errorStr);
                uiwait(h);          % Wait until the user presses the OK button  
                return;
            end

            % Make sure the start and end frames fit to the selected movie frames
            if (guiData.moviefirstimg < startFrame)
                errorStr = ['The selected movie start frame (' num2str(guiData.moviefirstimg) ') does not fit with the job start frame (' num2str(startFrame) ')'];
                h = errordlg(errorStr);
                uiwait(h);          % Wait until the user presses the OK button  
                return;
            end
            if (guiData.movielastimg > endFrame)
                errorStr = ['The selected movie end frame (' num2str(guiData.movielastimg) ') does not fit with the job end frame (' num2str(endFrame) ')'];
                h = errordlg(errorStr);
                uiwait(h);          % Wait until the user presses the OK button  
                return;
            end
        end

        % Make sure increment is consistent over movies
        if jobCount == 1
            prevIncrement = increment;
        else
            if increment ~= prevIncrement
                errorStr = ['The increment value is different between the ' length(allMPM) ' movies.'];
                h = errordlg(errorStr);
                uiwait(h);          % Wait until the user presses the OK button  
                return;
            end
        end  % if jobCount == 1

        % Get frame interval and pixel length
        frameInterval = round (jobData(jobCount).timeperframe / 60);    % In minutes
        pixelLength = jobData(jobCount).mmpixel;

        % Check that values are consistent over movies
        if jobCount == 1
            prevFrameInterval = frameInterval;
            prevPixelLength = pixelLength;
        else
            if frameInterval ~= prevFrameInterval
                errorStr = ['The frame rate is different between the ' num2str(length(allMPM)) ' movies.'];
                h = errordlg(errorStr);
                uiwait(h);          % Wait until the user presses the OK button  
                return;
            end
            if pixelLength ~= prevPixelLength
                errorStr = ['The pixel length is different between the ' num2str(length(allMPM)) ' movies.'];
                h = errordlg(errorStr);
                uiwait(h);          % Wait until the user presses the OK button  
                return;
            end
        end  % if jobCount == 1

        % Get row and colsizes
        rowSize = jobData(jobCount).rowsize;
        colSize = jobData(jobCount).colsize;

        % This is only needed for the convex hull calculations, so if these are
        % not needed, do not test this
        convexHullPlotNeeded = get (handles.checkbox_avg_convex_hull_area,'Value');

        % Make sure row and colsize is consistent over movies
        if convexHullPlotNeeded
            if jobCount == 1
                prevRowSize = rowSize;
                prevColSize = colSize;
            else
                if rowSize ~= prevRowSize
                    errorStr = ['The image row size is different between the ' num2str(length(allMPM)) ' movies.'];
                    h = errordlg(errorStr);
                    uiwait(h);          % Wait until the user presses the OK button  
                    return;
                end
                if colSize ~= prevColSize
                    errorStr = ['The image column size is different between the ' num2str(length(allMPM)) ' movies.'];
                    h = errordlg(errorStr);
                    uiwait(h);          % Wait until the user presses the OK button  
                    return;
                end
            end  % if jobCount == 1
        end   % if convexHullPlotNeeded
    end  % for jobCount = 1 : length(allMPM)

    % Assign all values to the handles struct
    handles.allMPM = allMPM;
    handles.allCellProps = allCellProps;
    handles.allClusterProps = allClusterProps;
    handles.allFrameProps = allFrameProps;
    handles.allValidFrames = allValidFrames;
    handles.jobData = jobData;
    handles.guiData = guiData;

    % Update GUI handles struct
    guidata (hObject,handles);
end

%--------------------------------------------------------------------

function [answer, empty] = directoryEmpty (path)
% Test if the save directory already contains files and in case it does,
% warn the user that they are going to be overwritten. Windows and
% linux/unix pc's have to be treated differently in this regard

% Get directory listing
dirList = dir ([path filesep '*']);

% Check whether directory is empty
if ~ispc   % linux/unix pc
    if (length (dirList) > 2) & (dirList(1).name == '.') & (dirList(2).name == '..')
       empty = 0;
    else
       empty = 1;
    end
else   % windows pc
    if (length (dirList) > 0)
       empty = 0;
    else
       empty = 1;
    end
end

% If not empty ask whether the user wants to empty it
answer = '';
if ~empty
   answer = questdlg(['Directory ' path ' is not empty. All existing files will be overwritten. Continue?']);
   if strcmp(answer,'Yes')
      delete ([path filesep '*']);
      empty = 1;
   else
      answer = '';
      empty = 0;
   end
end

%--------------------------------------------------------------------

function radioButtons = getRadiobuttonValues (handles)
% Assign the radiobutton values to the radioButtons struct

% Gett cell/cluster stats values
radioButtons.cellclusterplot = get (handles.checkbox_clustercellstats,'Value');
   radioButtons.cellclusterplot_1 = get (handles.checkbox_amount_cells,'Value');
   radioButtons.cellclusterplot_2 = get (handles.checkbox_percentage_cells,'Value');
   
% Get area stats values
radioButtons.areaplot = get (handles.checkbox_areastats,'Value');
   radioButtons.areaplot_1 = get (handles.checkbox_total_area,'Value');
   radioButtons.areaplot_2 = get (handles.checkbox_single_cluster_area,'Value');
   radioButtons.areaplot_3 = get (handles.checkbox_avg_convex_hull_area,'Value');
  
% Get perimeter stats values
radioButtons.perimeterplot = get(handles.checkbox_perimeter,'Value');

% Get velocity values
radioButtons.speedplot = get(handles.checkbox_speed,'Value');
   radioButtons.speedplot_1 = get(handles.checkbox_all_to_single_speed,'Value');
   radioButtons.speedplot_2 = get(handles.checkbox_average_speed,'Value');
   radioButtons.speedplot_3 = get(handles.checkbox_speed_variance,'Value');
   radioButtons.speedplot_4 = get(handles.checkbox_speed_histogram,'Value');
   
% Get cell/cell distance values
radioButtons.cellcelldistplot = get(handles.checkbox_cellcelldisthist,'Value');
   radioButtons.cellcelldistplot_1 = get(handles.checkbox_avg_distance_cells,'Value');
   
% Get neighbourhood stats values
radioButtons.neighbourplot = get(handles.checkbox_neighbourhood,'Value');
   radioButtons.neighbourplot_1 = get(handles.checkbox_nb_trajectories,'Value');
   radioButtons.neighbourplot_2 = get(handles.checkbox_nb_interact,'Value');
   
% Get Ripley choas values
radioButtons.ripleyplot = get(handles.GUI_ripley_cb,'Value');
   radioButtons.ripleyplot_1 = get(handles.GUI_chaosstats_cb,'Value');
   
% Get histogram radiobutton values
radioButtons.allcellshist = get (handles.GUI_vel_all_cells_cb,'Value');
radioButtons.singlecellshist = get (handles.GUI_vel_single_cells_cb,'Value');
radioButtons.clusteredcellshist = get (handles.GUI_vel_clust_cells_cb,'Value');

% Get button for show/not show plots
radioButtons.donotshowplots = get (handles.GUI_notshowplots_cb,'Value');

% Get button for running average
radioButtons.runningaverage = get (handles.GUI_running_average_cb,'Value');

% Get movie buttons
radioButtons.movieinclnuclei = get (handles.GUI_fm_inclcentromers_rb,'Value');
radioButtons.movieincltracks = get (handles.GUI_fm_incltracks_rb,'Value');
radioButtons.moviesaveastiff = get (handles.GUI_fm_saveastiff_rb,'Value'); 
radioButtons.movietypeavi = get (handles.GUI_movietype_avi_rb,'Value');
radioButtons.movietypeqt = get (handles.GUI_movietype_qt_rb,'Value');

% Get button that shows whether we should start counting from 1
radioButtons.alwayscount1 = get (handles.GUI_alwayscount1_cb,'Value');

% Get button for plotting the estimated function
radioButtons.plotestimate = get (handles.GUI_plotestimate_cb,'Value');

%--------------------------------------------------------------------

% --- Executes on button press in GUI_notshowplots_cb.
function GUI_notshowplots_cb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_notshowplots_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_notshowplots_cb

%--------------------------------------------------------------------

% --- Executes on button press in GUI_running_average_cb.
function GUI_running_average_cb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_running_average_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_running_average_cb

%--------------------------------------------------------------------

function GUI_windowsize_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_windowsize_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_windowsize_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_windowsize_ed as a double
handles = guidata (hObject);

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get (hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.guiData.windowsize = 5;
    set (handles.GUI_windowsize_ed, 'Value', handles.guiData.windowsize);  % Revert the value back
    return
else
    handles.guiData.windowsize = val;
end

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_windowsize_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_windowsize_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%--------------------------------------------------------------------

function GUI_maxcellcelldist_ed_Callback(hObject, eventdata, handles)
handles = guidata (hObject);

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get (hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.guiData.maxcellcelldist = handles.defaultPostPro.maxcellcelldist;
    set (handles.GUI_maxcellcelldist_ed, 'String', num2str(handles.guiData.maxcellcelldist));  % Revert the value back
    return
else
    handles.guiData.maxcellcelldist = val;
end

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------

function GUI_maxcellcelldist_ed_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------

function GUI_chaosstats_cb_Callback(hObject, eventdata, handles)
handles = guidata(hObject);

val = get (hObject,'Value');

if val == 1
   % Checkbox is selected, so select all children
   set (handles.GUI_ripley_cb, 'Value', 1);
else  % val == 0
   % Checkbox was unselected so unselect all the children
   set (handles.GUI_ripley_cb, 'Value', 0);
end

% Update handles structure
guidata(hObject, handles); 

%--------------------------------------------------------------------

function GUI_ripley_cb_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------

function GUI_alwayscount1_cb_Callback(hObject, eventdata, handles)
handles = guidata(hObject);

% Get the job list and the current job
fileList = get(handles.GUI_filelist_lb,'String');
filesSelected = get(handles.GUI_filelist_lb,'Value');

% GUI values can only be updated for 1 job at the time, so in case we have
% more, do nothing
if ischar(fileList) | length(filesSelected) > 1
    return;
end

% Retrieve the data for the selected jobs
[allMPM, allCellProps, allClusterProps, allFrameProps, allValidFrames, jobData, result] = ptRetrieveJobData (fileList, 'all');

% Get radioButton values
alwaysCountFrom1 = get (hObject,'Value');

% Default plot ranges will be made similar to the frame range
if ~alwaysCountFrom1
  handles.guiData.plotfirstimg = allValidFrames{filesSelected}(1,1);
  handles.guiData.plotlastimg = allValidFrames{filesSelected}(1,end);
else
  handles.guiData.plotfirstimg = 1;
  handles.guiData.plotlastimg = allValidFrames{filesSelected}(1,end) - allValidFrames{filesSelected}(1,1) + 1;
end

% Set movie frame start and end
handles.guiData.moviefirstimg = allValidFrames{filesSelected}(1,1);
handles.guiData.movielastimg = allValidFrames{filesSelected}(1,end);

% Set values on the GUI
%ptSetPostproGUIValues (handles, filesSelected);
handles = ptSetPostproGUIValues (handles, 1);

% Update handles structure
guidata(hObject, handles); 

%--------------------------------------------------------------------

function GUI_chaossim_cb_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------

function GUI_ripconfint_ed_Callback(hObject, eventdata, handles)
handles = guidata (hObject);

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get (hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.guiData.ripleyconfint = handles.defaultPostPro.ripconfint;
    set (handles.GUI_ripconfint_ed, 'String', num2str(handles.guiData.ripleyconfint));  % Revert the value back
    return
else
    handles.guiData.ripleyconfint = val;
end

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------

function GUI_ripconfint_ed_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------

function GUI_plotestimate_cb_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------------

function pp_bad_frames_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------

function pp_bad_frames_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------

function GUI_fm_saveastiff_rb_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------

function GUI_selectallcells_pb_Callback(hObject, eventdata, handles)
handles = guidata (hObject);

% Check that files have been selected before
if ~isfield (handles, 'allMPM')  
    errorStr = ['Jobs should be selected first by using the Select button!'];
    h = errordlg(errorStr);
    uiwait(h);          % Wait until the user presses the OK button  
    return;
end

% Make sure only 1 job is selected
if length(handles.allMPM) > 1
    errorStr = ['Selecting cells can only be done for one job at a time. Please select only one job from the list.'];
    h = errordlg(errorStr);
    uiwait(h);          % Wait until the user presses the OK button  
    return;
end

% To select all cells we have to set the selectedcells field to []
handles.jobData(1).selectedcells = [];

% Let the user know
msgbox('All cells have been selected.');

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------

function selectedCellsMatrix = fillSelectedCellsMatrix (selectedCells, MPM, frameNr, initMatrix, prevSelectedCells)
% This function filles the matrix for selected cells based on the cells
% that have been selected from one frame in the movie

% Check whether a new matrix has to be created
if initMatrix
    selectedCellsMatrix = zeros(size(MPM,1),size(MPM,2)/2);
else
    selectedCellsMatrix = prevSelectedCells;
end

% Start adding cells to selectedCellsMatrix: for this we have to find
% the first free row
freeCount = 1;
freeRow = find(selectedCellsMatrix(freeCount,:));
while ~isempty(freeRow)
    freeCount = freeCount + 1;
    freeRow = find(selectedCellsMatrix(freeCount,:));
end

% Walk through the selected cells and find tracks in all frames
for iCount = freeCount : freeCount+length(selectedCells)-1

    % Find out which parts of the MPM we can throw away (we only want the
    % tracks of the selected cells
    %firstRow = find(MPM(frameNr*2)) + selectedCells(iCount-freeCount+1) - 1;

    % Cut out that part of the MPM
    %tempMPM = MPM(firstRow,:);
    tempMPM = MPM(selectedCells,:);

    % Figure out where the first whole block of coordinates is
    zeroIndx = find(tempMPM == 0);
    nonZeroIndx = find(tempMPM);
    zeroIndx = zeroIndx(find(zeroIndx > nonZeroIndx(1)));
    %tempMPM = tempMPM(nonZeroIndx(2):2:zeroIndx(1)-1);

    % Loop through this row of coordinates and fill the selectedCellsMatrix
    if ~isempty(zeroIndx)
        for jCount = nonZeroIndx(2) : 2 : zeroIndx(1)-1
            selectedCellsMatrix(iCount,jCount/2) = selectedCells(iCount-freeCount+1);
        end
    else
        for jCount = nonZeroIndx(2) : 2 : nonZeroIndx(end)
            selectedCellsMatrix(iCount,jCount/2) = selectedCells(iCount-freeCount+1);
        end
    end
end

%--------------------------------------------------------------------

function GUI_fm_fulltracks_cb_Callback(hObject, eventdata, handles)
handles = guidata(hObject);

val = get (hObject,'Value');

if val == 1
   % Button was clicked so should become 1
   handles.guiData.fulltracks = 1;
   set (handles.GUI_fm_fulltracks_cb, 'Value', 1);
   
   % The QT button should automatically become 0
   set (handles.GUI_fm_incltracks_rb, 'Value', 0);
   handles.guiData.dragtracks = 0;
else
   % Button was already selected, but make sure that movietype is set
   % correctly
   handles.guiData.fulltracks = 0;
   %handles.guiData.dragtracks = 0;
end

% Update handles structure
guidata(hObject, handles); 

%--------------------------------------------------------------------

function text_polydatadir_ptpp_Callback(hObject, eventdata, handles)
handles = guidata(hObject);

% Get input from the gui and assign it to the handle; if directory does not
% exist, ask user whether to create it
polyDataDirectory = get(hObject,'String');

% If a filesep exist at the end of the path, remove it
if polyDataDirectory(end) == '/'
    polyDataDirectory(end) = '';
end

% If the path doesn't exist ask user if it should be created
if ~exist(polyDataDirectory, 'dir')
   msgStr = ['This directory does not exist yet. Do you want to create it?'];
   answer = questdlg(msgStr, 'Create Directory', 'Yes', 'No', 'Yes');
   if strcmp(answer,'Yes')
      mkdir (polyDataDirectory);
   else
      polyDataDirectory = handles.polyDataDirectory;
      set(handles.text_polydatadir_pt,'String',polyDataDirectory);
   end
end

% Assign it to the guiData structure
handles.polyDataDirectory = polyDataDirectory;

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------

function text_polydatadir_ptpp_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------

function pb_polydatabrowse_ptpp_Callback(hObject, eventdata, handles)
handles = guidata(hObject);

% Retrieve the directory and filename where to save the result
polyDataDirectory = get(handles.text_polydatadir_ptpp,'String');

directory = uigetdir(polyDataDirectory,'Please select a new POLYDATA directory');

% Do nothing in case the user doesn't select anything
if directory == 0
   return
end

% If a filesep exist at the end of the path, remove it
if directory(end) == '/'
    directory(end) = '';
end

% % Check the directory structure under BIODATA and create projects subdir 
% % if not existent
% if ~exist([directory filesep 'projects'],'dir')
%     try
%         mkdir(directory,'projects');
%     catch
%         fprintf (1, 'Error: cannot create projects directory under %s. Please create it manually.\n', polyDataDirectory);
%     end
% end
% 
% % Do the same for the experiments subdir
% if ~exist([directory filesep 'experiments'],'dir')
%     try
%         mkdir(directory,'experiments');
%     catch
%         fprintf (1, 'Error: cannot create experiments directory under %s. Please create it manually.\n', polyDataDirectory);
%     end
% end

% Update the biodata text field
set(handles.text_polydatadir_ptpp, 'String',directory);

% Update polyDataDirectory in the handles
handles.polyDataDirectory = directory;

% Update GUI handles struct
guidata (hObject,handles);

% --------------------------------------------------------------------

function GUI_drugtimepoint_ed_Callback(hObject, eventdata, handles)
handles = guidata (hObject);

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get (hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.guiData.ripleyconfint = handles.defaultPostPro.drugtimepoint;
    set (handles.GUI_drugtimepoint_ed, 'String', num2str(handles.guiData.drugtimepoint));  % Revert the value back
    return
else
    handles.guiData.drugTimepoint = val;
end

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------

function GUI_drugtimepoint_ed_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


