function varargout = PolyTrack(varargin)
% PolyTrack M-file for PolyTrack.fig
%
% PolyTrack  contains the callback functions for the Polytrack GUI. Next to
%            Matlab standard initialization code and functions, a number of
%            callbacks have been implemented that control the processing of
%            the data files
%
% SYNOPSIS   varargout = PolyTrack(varargin)
%
% INPUT      varargin (optional)
%
% OUTPUT     varargout
%
% Last Modified by A. Kerstens 01-Mar-2004

% This is matlab stuff we should not touch.
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_start_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_start_OutputFcn, ...
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

%-----------------------------------------------------------------------------

% --- Executes just before GUI_start is made visible.
function GUI_start_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_start (see VARARGIN)

% Choose default command line output for GUI_start
handles.output = hObject;

% Create a default job structure that defines all the fields used in the program
defaultjob = struct('imagedirectory', [], 'imagename', [], 'firstimage', 1, 'lastimage', [],...
                   'increment', 1, 'savedirectory', [], 'fi_background', [], 'fi_nucleus', [],...
                   'la_background', [], 'la_nucleus', [], 'maxsearch', 48, 'mmpixel', [], 'mmpixel_index', 1,...
                   'minsize',300, 'maxsize', 1500, 'minsdist', 30, 'fi_halolevel', [], 'la_halolevel', [],...
                   'minedge', 10, 'sizetemplate', 41, 'boxsize', 141, 'noiseparameter', 0.15,...
                   'mincorrqualtempl', 0.2, 'leveladjust', 0.7, 'timestepslide', 5, 'mintrackcorrqual', 0.5,...
                   'coordinatespicone', [], 'intensityMax', 4095, 'bitdepth', 12, 'bitdepth_index', 3, 'bodyname', [],...
                   'imagenameslist', [], 'timeperframe', [], 'clustering', 1, 'minmaxthresh', 0,...
                   'timestepslide_index', 2);

% Assign the default job values to the GUI handle so it can be passed around
handles.defaultjob = defaultjob;

% Set the colors of the gui
set(hObject,'Color',[0,0,0.627]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_start wait for user response (see UIRESUME)
% uiwait(handles.polyTrack_mainwindow);

%-----------------------------------------------------------------------------

% --- Outputs from this function are returned to the command line.
function varargout = GUI_start_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Here is were we output stuff to the matlab command line

% Get default command line output from handles structure
varargout{1} = handles.output;

%-----------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_job_lb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_job_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

% ispc is how we test whether we run on a windows machine
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%------------------------------------------------------------------------------

% --- Executes on selection change in GUI_st_job_lb.
function GUI_st_job_lb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_job_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GUI_st_job_lb contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GUI_st_job_lb
handles = guidata(hObject);

% Get the number of the currently selected project in the list
projNum = get(hObject,'Value');

% Use the values of this project to select the correct job and fill the
% text fields of the GUI
fillFields(handles,handles.jobs(projNum))

%-------------------------------------------------------------------------------

% --- Executes on button press in GUI_st_addjob_pb.
function GUI_st_addjob_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_addjob_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Make sure the handles struct can be used in this function
handles = guidata(hObject);

% Get a handle to the GUI job list
listhandle = handles.GUI_st_job_lb;

% Select an image filename or a file called 'jobvalues.mat' from a user selected directory
[filename,imagedirectory] = uigetfile({'*.tif','TIFF-files';'jobvalues.mat','Saved Values';'*.*','all files'},...
                                      'Please select a TIFF image file or jobvalues.mat file');

% Do nothing in case the user doesn't select anything
if filename == 0
    return
end

% Else go on and check whether a file called 'jobvalues.mat' has been selected
gotvals = 0;
if strcmp(filename,'jobvalues.mat')
    cd (imagedirectory);
    load jobvalues.mat;
    filename = jobvalues.imagename;
    gotvals = 1;
end

% Get the current job list
jobList = get(handles.GUI_st_job_lb,'String');

% If jobList consists of more than one entry append the new filename after the
% last one else just put it in the jobList as the first entry
if ~iscell(jobList)
    jobList = cellstr(filename);
else
    jobList(end+1) = cellstr(filename);
end

% Store the modified job list back in the GUI handle
set(handles.GUI_st_job_lb,'String',jobList);

% Get the last job in the job list
projNum = length(jobList); 

% In case a jobvalues.mat file was read, store this in the handle
if gotvals == 1
    handles.jobs(projNum) = jobvalues;
    clear jobvalues
% else start using the default job values and do some more
else
    handles.jobs(projNum) = handles.defaultjob;
    handles.jobs(projNum).imagedirectory = imagedirectory;
    handles.jobs(projNum).imagename = filename;
	    
    % Now we have to do the following:
    % Find out what part of the filename describes the images and which part
    % is just counting them through.
    % We beginn by starting from the rear end of filename (after cutting off
    % the extension) and converting the last few elements (starting with
    % only one and then adding on one each loop) from a string to a number.
    % As soon as the answer of the conversion is NaN (not a number) we have struck the
    % first letter. Furthermore we say that max three digits are considered
    % as numbering.
    
    number = 0;
    countNum = 0;
    while ~isnan(number) & (countNum <3)
         countNum = countNum+1;
         number = str2num(filename(end-(4+countNum):end-4));
    end

    % Extract the body of the filename and store in handles struct
    handles.jobs(projNum).bodyname = filename(1:(end-(4+countNum)));
    bodyname = handles.jobs(projNum).bodyname;
     
    % Select the current project
    set(handles.GUI_st_job_lb, 'Value', projNum);
	
    % Create a list of files present in the image directory selected by the user
    dirList = dir(imagedirectory);
    dirList = struct2cell(dirList);
    dirList = dirList(1,:);
    
    % Find all files within this directory with the same name as the selected filename
    ind = strmatch(handles.jobs(projNum).bodyname, dirList);
    dirList = dirList(ind)';
    handles.jobs(projNum).lastimage = length(dirList);
      
    % Sort the images by successive numbers:
    % First we get all numbers and write them into a vector
    for jRearange = 1:length(dirList)
        tmpName = char(dirList(jRearange));
        imageNum(jRearange) = str2num(tmpName(length(handles.jobs(projNum).bodyname)+1:end-4));
    end
    
    % Then we sort that vector and sort the dirList accordingly
    [junk,indVec] = sort(imageNum);
    handles.jobs(projNum).imagenameslist = dirList(indVec);
        
    % Create a directory to save the details and results of this job
    % Note: we call the directory results + bodyname + seq number
    % number will be lowest unoccupied number for this specific directory name
    cd (imagedirectory)
    done = 0;
    counter = 1;
    while done == 0
        newdirname = [];
        newdirname = ['results', bodyname, num2str(counter)];
           
        % Loop on untill we find an unoccupied new dirname
        if exist (newdirname, 'dir') == 0
            mkdir (imagedirectory, newdirname);
            tempname = [imagedirectory, newdirname];
            mkdir(tempname,'body');
             
            handles.jobs(projNum).savedirectory = [imagedirectory, newdirname];
            done = 1;
        end
        counter = counter + 1;
    end
end

% Update GUI handle struct
guidata(hObject, handles);
%handles = guidata(hObject);

% Last but not least make sure the text field on the GUI show the latest% values
fillFields(handles, handles.jobs(projNum))

%-------------------------------------------------------------------------------

% --- Executes on button press in GUI_st_deletejob_pb.
function GUI_st_deletejob_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_deletejob_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

listnotfilled = [];
handles = guidata(hObject);

% Get the job list and the current job
jobList = get(handles.GUI_st_job_lb,'String');
projNum = get(handles.GUI_st_job_lb,'Value');

% Joblist will only then be a cell, if there are jobs in it.
% Otherwise it is a string (No project loaded)
if ~iscell(jobList)
    % Show an error dialog with an appropriate message and wait
    % for the user to press a button
    h=errordlg('Sorry, there are no jobs to delete.');
    uiwait(h);
    return
end

% Put a standard string in the job list window if there is no project to show
% and delete the currently selected job from the gui
if length(jobList) == 1
    jobList = char('No project loaded');
    listnotfilled = 1;
else
    jobList(projNum) = [];
end

% Set the list to the first project to be on the safe side
% And store new jobList in gui handle
set(handles.GUI_st_job_lb,'Value',1);
set(handles.GUI_st_job_lb,'String',jobList);

% Store job data
handles.jobs(projNum) = [];
guidata(hObject,handles);

% Show the data of the first job in the list, or if no job is present, 
% show the data of the defaultjob
if isempty (listnotfilled)
    fillFields (handles, handles.jobs(1))
else
    fillFields (handles, handles.defaultjob)
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_path_imagedirectory_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_imagedirectory_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% On Windows computers the background color is usually set to white while
% on any other platform the default is taken.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_path_imagedirectory_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_imagedirectory_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_path_imagedirectory_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_path_imagedirectory_ed as a double

% Retrieve handle data from the gui
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Get the image directory name
imagedir = get(hObject,'String');

% Retrieve the current job number and assign it the image directory value
projNum = get(handles.GUI_st_job_lb,'Value');
handles.jobs(projNum).imagedirectory =  imagedir;

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_path_imagename_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_imagename_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% On Windows computers the background color is usually set to white while
% on any other platform the default is taken.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_path_imagename_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_imagename_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_path_imagename_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_path_imagename_ed as a double

% handles = guidata(hObject);
% 
% 
% imagename = get(hObject,'String');
% 
% % Select the current job
% projNum = get(handles.GUI_st_job_lb,'Value');
% 
% handles.jobs(projNum).imagename =  imagename;
% 
% 
% % Update handles structure
% guidata(hObject, handles);
% 
% %%%%%%%%save altered values to disk%%%%%%%%%%%%
% cd (handles.jobs(projNum).savedirectory)
% jobvalues = handles.jobs(projNum);
% save ('jobvalues','jobvalues')
% clear jobvalues
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_path_firstimage_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_firstimage_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_path_firstimage_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_firstimage_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_path_firstimage_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_path_firstimage_ed as a double
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Get the currently selected project
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).firstimage = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).firstimage = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_path_lastimage_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_lastimage_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_path_lastimage_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_lastimage_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_path_lastimage_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_path_lastimage_ed as a double
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Get the currently selected project
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).lastimage = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).lastimage = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes on button press in GUI_st_iq_set_pb.
function GUI_st_iq_set_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_set_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if iscell(jobList)
   nrofjobs = length(jobList); 
else
   h=errordlg('At least one job should be loaded first, before thresholds can be determined.');
   uiwait(h);
   return
end 

% Determine the different thresholds of the images interactively
leveldeterminer(hObject);

% Get the currently selected job
handles = guidata(hObject);
projNum = get(handles.GUI_st_job_lb,'Value');

% Fill all the threshold input fields with the newly found values
fillFields(handles,handles.jobs(projNum))

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_iq_fi_background_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_fi_background_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_iq_fi_background_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_fi_background_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_iq_fi_background_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_iq_fi_background_ed as a double
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).fi_background = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).fi_background = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_iq_fi_nucleus_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_fi_nucleus_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_iq_fi_nucleus_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_fi_nucleus_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_iq_fi_nucleus_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_iq_fi_nucleus_ed as a double
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).fi_nucleus = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).fi_nucleus = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_iq_la_background_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_la_background_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_iq_la_background_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_la_background_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_iq_la_background_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_iq_la_background_ed as a double
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).la_background = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).la_background = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_iq_la_nucleus_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_la_nucleus_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_iq_la_nucleus_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_la_nucleus_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_iq_la_nucleus_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_iq_la_nucleus_ed as a double
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).la_nucleus = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).la_nucleus = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_bp_maxsearch_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_maxsearch_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor', get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_bp_maxsearch_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_maxsearch_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_bp_maxsearch_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_bp_maxsearch_ed as a double
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).maxsearch = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).maxsearch = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_eo_timestepslide_pm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_timestepslide_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on
% Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor', get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

% --- Executes on selection change in GUI_st_eo_timestepslide_pm.
function GUI_st_eo_timestepslide_pm_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_timestepslide_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GUI_st_eo_timestepslide_pm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GUI_st_eo_timestepslide_pm
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get the value and index from the timestepslide popup menu and assign to handles struct
val = get(hObject, 'Value');
list = get(hObject, 'String');
selected_val = list{val};
handles.jobs(projNum).timestepslide = str2double(selected_val);
handles.jobs(projNum).timestepslide_index = val;

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_bp_mmpixel_pm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_mmpixel_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

% --- Executes on selection change in GUI_st_bp_mmpixel_pm.
function GUI_st_bp_mmpixel_pm_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_mmpixel_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GUI_st_bp_mmpixel_pm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GUI_st_bp_mmpixel_pm
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get the value and index from the mm-pixel popup menu and assign to handles struct
val = get(hObject, 'Value');
list = get(hObject, 'String');
selected_val = list{val};
handles.jobs(projNum).mmpixel = str2double(selected_val);
handles.jobs(projNum).mmpixel_index = val;

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_bp_minsize_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_minsize_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_bp_minsize_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_minsize_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_bp_minsize_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_bp_minsize_ed as a double
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).minsize = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).minsize = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_bp_maxsize_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_maxsize_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_bp_maxsize_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_maxsize_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_bp_maxsize_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_bp_maxsize_ed as a double
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).maxsize = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).maxsize = val;
end

% Update handles structure
guidata(hObject, handles);


% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_bp_minsdist_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_minsdist_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_bp_minsdist_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_minsdist_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_bp_minsdist_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_bp_minsdist_ed as a double
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).minsdist = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).minsdist = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes on button press in GUI_st_bp_setminsize_pb.
function GUI_st_bp_setminsize_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_setminsize_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

SetCellValues(hObject,1);
 
%-------------------------------------------------------------------------------

% --- Executes on button press in GUI_st_bp_setmaxsize_pb.
function GUI_st_bp_setmaxsize_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_setmaxsize_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

SetCellValues(hObject,2);

%-------------------------------------------------------------------------------

% --- Executes on button press in GUI_st_bp_setmindist_pb.
function GUI_st_bp_setmindist_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_setmindist_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

SetCellValues(hObject,3);

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_eo_minedge_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_minedge_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_eo_minedge_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_minedge_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_eo_minedge_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_eo_minedge_ed as a double
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).minedge = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).minedge = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_eo_noiseparameter_pm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_noiseparameter_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

% --- Executes on selection change in GUI_st_eo_noiseparameter_pm.
function GUI_st_eo_noiseparameter_pm_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_noiseparameter_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GUI_st_eo_noiseparameter_pm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GUI_st_eo_noiseparameter_pm
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).noiseparameter = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).noiseparameter = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_eo_mincorrqualtempl_pm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_mincorrqualtempl_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

% --- Executes on selection change in GUI_st_eo_mincorrqualtempl_pm.
function GUI_st_eo_mincorrqualtempl_pm_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_mincorrqualtempl_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GUI_st_eo_mincorrqualtempl_pm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GUI_st_eo_mincorrqualtempl_pm
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).mincorrqualtempl = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).mincorrqualtempl = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes on button press in GUI_st_bp_loadsettings_pb.
function GUI_st_bp_loadsettings_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_loadsettings_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be loaded...');
   uiwait(h);
   return
else
   % Select the current job
   projNum = get(handles.GUI_st_job_lb,'Value');

   imgDir = handles.jobs(projNum).imagedirectory;
   imgNam = handles.jobs(projNum).imagename;
   firstImg = handles.jobs(projNum).firstimage;
   lastImg = handles.jobs(projNum).lastimage;
   increment = handles.jobs(projNum).increment;
   savedirectory = handles.jobs(projNum).savedirectory;

   cd (handles.jobs(1).savedirectory);

   [filename,jobValPath] = uigetfile({'*.mat','mat-files'},'Please select a jobvalues file');

   if ~strcmp(filename,'jobvalues.mat')
      h = errordlg('Select a file named jobvalues.mat please...');
      uiwait(h);          % Wait until the user presses the OK button
      return
   end

   % Load the jobvalues into the jobs struct of the current job
   cd (jobValPath)
   load ('jobvalues.mat');
   handles.jobs(projNum) = jobvalues;

   handles.jobs(projNum).imagedirectory = imgDir;
   handles.jobs(projNum).imagename = imgNam;
   handles.jobs(projNum).firstimage = firstImg;
   handles.jobs(projNum).lastimage = lastImg;
   handles.jobs(projNum).increment = increment;
   handles.jobs(projNum).savedirectory = savedirectory;

   fillFields(handles,handles.jobs(projNum))
end

% Update handles structure
guidata(hObject, handles);

%-------------------------------------------------------------------------------

% --- Executes on button press in GUI_st_bp_savesettings_pb.
function GUI_st_bp_savesettings_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_savesettings_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be saved...');
   uiwait(h);
   return
else
   % Select the current job
   projNum = get(handles.GUI_st_job_lb,'Value');

   % Store the latest data in jobvalues.mat in the specified save directory
   if ~isempty(handles.jobs(projNum).savedirectory)
      cd (handles.jobs(projNum).savedirectory)
      jobvalues = handles.jobs(projNum);
      save ('jobvalues','jobvalues')
      clear jobvalues
   end
end

% Update handles structure
guidata(hObject, handles);


%-------------------------------------------------------------------------------

% --- Executes on button press in GUI_st_bp_defaultsettings_pb.
function GUI_st_bp_defaultsettings_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bp_defaultsettings_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any default settings can be set...');
   uiwait(h);
   return
else
   % Select the current job
   projNum = get(handles.GUI_st_job_lb,'Value');

   % Store the fields that shouldn't be defaulted in temp vars
   imgDir = handles.jobs(projNum).imagedirectory;
   imgNam = handles.jobs(projNum).imagename;
   firstImg = handles.jobs(projNum).firstimage;
   lastImg = handles.jobs(projNum).lastimage;
   increment = handles.jobs(projNum).increment;
   savedirectory = handles.jobs(projNum).savedirectory;

   % Set the default values to all fields
   handles.jobs(projNum) = handles.defaultjob;

   % Retrieve the values that we temporarily stored before
   handles.jobs(projNum).imagedirectory = imgDir;
   handles.jobs(projNum).imagename = imgNam;
   handles.jobs(projNum).firstimage = firstImg;
   handles.jobs(projNum).lastimage = lastImg;
   handles.jobs(projNum).increment = increment;
   handles.jobs(projNum).savedirectory = savedirectory;

   % Update all the fields on the gui
   fillFields(handles,handles.jobs(projNum))
end 

% Update handles structure
guidata(hObject, handles);

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_eo_leveladjust_pm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_leveladjust_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

% --- Executes on selection change in GUI_st_eo_leveladjust_pm.
function GUI_st_eo_leveladjust_pm_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_leveladjust_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GUI_st_eo_leveladjust_pm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GUI_st_eo_leveladjust_pm
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).leveladjust = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).leveladjust = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes on button press in GUI_st_run_pb.
function GUI_st_run_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_run_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if iscell(jobList)
   % nrofjobs is basically equal to the last entry nr in the job list which
   % again is equal to the length of the list
   nrofjobs = length(jobList); 
else
   h=errordlg('A job should be loaded first, before a run can be started.');
   uiwait(h);
   return
end

% Loop through all the jobs in the joblist
for projNum = 1:nrofjobs

   % Get the image to image increment step
   Increment = handles.jobs(projNum).increment;

   % Get the first image number
   possibleImg = handles.jobs(projNum).firstimage;
   while (possibleImg + Increment) <= handles.jobs(projNum).lastimage
      possibleImg = possibleImg + Increment;
   end

   % Make sure the last image nr fits with the last image found in the previous loop
   % If these are not the same, adjust last image number
   if handles.jobs(projNum).lastimage > possibleImg
      handles.jobs(projNum).lastimage = possibleImg;
   end

   % Save the definite version of jobvalues
   if ~isempty(handles.jobs(projNum).savedirectory)
      cd (handles.jobs(projNum).savedirectory)
      jobvalues = handles.jobs(projNum);
      save ('jobvalues','jobvalues')
      clear jobvalues
   end
        
   % Here's where the real tracking process starts for the selected job
   % AK: the try-catch should be uncommented as soon as testing is done!!!
%   try
      trackCells (hObject,projNum);
%   catch    
%     errordlg(['job number ',num2str(projNum),' had an error and could not be completed'])
%     h=errordlg(['job number ',num2str(projNum),' had an error and could not be completed',lasterr]);
%     uiwait(h);
% %   disp(lasterr)
%   end
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_eo_sizetemplate_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_sizetemplate_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_eo_sizetemplate_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_sizetemplate_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_eo_sizetemplate_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_eo_sizetemplate_ed as a double
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).sizetemplate = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).sizetemplate = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_eo_mintrackcorrqual_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_mintrackcorrqual_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_eo_mintrackcorrqual_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_mintrackcorrqual_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_eo_mintrackcorrqual_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_eo_mintrackcorrqual_ed as a double
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).mintrackcorrqual = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).mintrackcorrqual = val;
end

% Update handles structure
guidata(hObject, handles);


% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_path_savedirectory_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_savedirectory_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_path_savedirectory_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_savedirectory_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_path_savedirectory_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_path_savedirectory_ed as a double
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Get the directory path and name
savedirectory = get(hObject,'String');

%if ~exist(savedirectory, 'file')
%  h=errordlg('The selected save directory does not exist. I will create it first...');
%  uiwait(h);
%  mkdir (savedirectory);
%end

% Select the current job and store the directory name in the struct
% AK: Shouldn't some sort of validity check be done here??
projNum = get(handles.GUI_st_job_lb,'Value');
handles.jobs(projNum).savedirectory =  savedirectory;

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes on button press in GUI_st_path_savedirectory_browse_pb.
function GUI_st_path_savedirectory_browse_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_savedirectory_browse_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Using the gui browser, select a directory name
savedirectory = uigetdir;

% If nothing has been selected, do nothing
if savedirectory == 0
    return
end

% Update the savedirectory field on the gui as well
set(handles.GUI_st_path_savedirectory_ed,'String',savedirectory);

% And store the directory in the handle struct
projNum = get(handles.GUI_st_job_lb,'Value');
handles.jobs(projNum).savedirectory =  savedirectory;

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_path_increment_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_increment_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_path_increment_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_increment_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_path_increment_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_path_increment_ed as a double
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).increment = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).increment = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes on button press in pushbutton22.
function GUI_st_test_pb_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if iscell(jobList)
   nrofjobs = length(jobList); 
else
   h=errordlg('At least one job should be loaded first, before test and initialization can be started.');
   uiwait(h);
   return
end 

% Start the test and initialization process
testbutton(hObject);

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_iq_fi_halolevel_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_fi_halolevel_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_iq_fi_halolevel_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_fi_halolevel_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_iq_fi_halolevel_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_iq_fi_halolevel_ed as a double
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).fi_halolevel = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).fi_halolevel = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_iq_la_halolevel_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_la_halolevel_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_iq_la_halolevel_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_iq_la_halolevel_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_iq_la_halolevel_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_iq_la_halolevel_ed as a double
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).la_halolevel = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).la_halolevel = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_bitdepth_pm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_bitdepth_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

% --- Executes on selection change in GUI_st_bitdepth_pm.
function GUI_st_bitdepth_pm_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_bitdepth_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GUI_st_bitdepth_pm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GUI_st_bitdepth_pm
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get the value and index from the bitdepth  popup menu and assign to handles struct
val = get(hObject, 'Value');
list = get(hObject, 'String');
selected_val = list{val};
handles.jobs(projNum).bitdepth = str2double(selected_val);
handles.jobs(projNum).bitdepth_index = val;

% Calculate the maximal value of the image, depending on it's bitdepth and
% store this info in the handles struct
bitdepth = str2double(selected_val);
handles.jobs(projNum).intensityMax =  2^bitdepth - 1;

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_path_timeperframe_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_timeperframe_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-------------------------------------------------------------------------------

function GUI_st_path_timeperframe_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_path_timeperframe_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_st_path_timeperframe_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_st_path_timeperframe_ed as a double

handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(projNum).timeperframe = [];
    fillFields(handles, handles.jobs(projNum))  % Revert the value back
    return
else
    handles.jobs(projNum).timeperframe = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes on button press in GUI_st_eo_clustering_rb.
function GUI_st_eo_clustering_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_clustering_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_st_eo_clustering_rb
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Get the value of the clustering radiobutton
val = get(hObject,'Value');

% Select the current job and store the radiobutton value
projNum = get(handles.GUI_st_job_lb,'Value');
handles.jobs(projNum).clustering =  val;

% Depending on the clustering value set the minmaxthreshold to the inverse
if val
    handles.jobs(projNum).minmaxthresh = 0;
else
    handles.jobs(projNum).minmaxthresh = 1;
end
    
% And set the value on the gui
set(handles.GUI_st_eo_minmaxthresh_rb,'Value',handles.jobs(projNum).minmaxthresh);

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes on button press in GUI_st_eo_minmaxthresh_rb.
function GUI_st_eo_minmaxthresh_rb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_eo_minmaxthresh_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GUI_st_eo_minmaxthresh_rb
handles = guidata(hObject);

% Get the list of jobs
jobList = get(handles.GUI_st_job_lb,'String');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if ~iscell(jobList)
   h=errordlg('At least one job should be loaded first, before any settings can be changed...');
   uiwait(h);
   return
end 

% Get the value of the clustering radiobutton
val = get(hObject,'Value');

% Select the current job and store value
projNum = get(handles.GUI_st_job_lb,'Value');
handles.jobs(projNum).minmaxthresh =  val;

% Depending on the minmaxthreshold value set the clustering radiobutton to the inverse
if val
    handles.jobs(projNum).clustering = 0;
else
    handles.jobs(projNum).clustering = 1;
end
    
% And set the value on the gui
set(handles.GUI_st_eo_clustering_rb,'Value',handles.jobs(projNum).clustering);

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

function File_menu_Callback(hObject, eventdata, handles)
% hObject    handle to File_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This is just the filemenu which doesn't do anything really; just keep the
% submenus like 'Exit'

%-------------------------------------------------------------------------------

function exit_menuitem_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Select the current job
projNum = get(handles.GUI_st_job_lb,'Value');

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(projNum).savedirectory)
   cd (handles.jobs(projNum).savedirectory)
   jobvalues = handles.jobs(projNum);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

% Destroy the GUI
delete(handles.polyTrack_mainwindow);

%-------------------------------------------------------------------------------
