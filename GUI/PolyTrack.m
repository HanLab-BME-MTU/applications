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
% CHANGE REVISION
%
% Name                  Date            Comment
% --------------------- ----------      -----------------------------------------------
% Colin Glass           Feb 04          Initial version
% Andre Kerstens        Jun 04          Changed default mincellsize to 300 (from 250)
% Andre Kerstens        Jul 04          Image file check is only done if it hasn't be done before.
% Andre Kerstens        Jul 04          Added frame properties matrix to ptTrackCells

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
                   'increment', 1, 'savedirectory', [], 'maxsearch', 81, 'mmpixel', 0.639,...
                   'minsize',300, 'maxsize', 4000, 'minsdist', 25, 'fi_halolevel', [], 'la_halolevel', [],...
                   'minedge', 10, 'sizetemplate', 41, 'boxsize', 141, 'noiseparameter', 0.15,...
                   'mincorrqualtempl', 0.5, 'leveladjust', 0.7, 'timestepslide', 5, 'mintracklength', 2,...
                   'coordinatespicone', [], 'intensityMax', 4095, 'bitdepth', 12, 'bitdepth_index', 3, 'bodyname', [],...
                   'imagenameslist', [], 'timeperframe', 300, 'timestepslide_index', 2);

% Assign the default job values to the GUI handle so it can be passed around
handles.defaultjob = defaultjob;

% Set the colors of the gui
set(hObject,'Color',[0,0,0.627]);

% Update handles structure
guidata(hObject, handles);

% Turn the resize and int conversion (matlab 7) warnings off
iptsetpref ('TrueSizeWarning', 'off');

% For matlab 7 turn int conversion warnings off
matlabVersion = version;
if (matlabVersion(1) == '7')
  intwarning ('off');
end

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
jobNumber = get(hObject,'Value');

% Use the values of this project to select the correct job and fill the
% text fields of the GUI
ptFillFields(handles,handles.jobs(jobNumber))

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
jobNumber = length(jobList); 

% In case a jobvalues.mat file was read, store this in the handle
if gotvals == 1
    handles.jobs(jobNumber) = jobvalues;
    clear jobvalues
% else start using the default job values and do some more
else
    handles.jobs(jobNumber) = handles.defaultjob;
    handles.jobs(jobNumber).imagedirectory = imagedirectory;
    handles.jobs(jobNumber).imagename = filename;
    
%     % Get the bitdepth from the gui just in case it is different from the
%     % default
%     val = get(handles.GUI_st_bitdepth_pm, 'Value');
%     list = get(handles.GUI_st_bitdepth_pm, 'String');
%     selected_val = list{val};
%     handles.jobs(jobNumber).bitdepth = str2double(selected_val);
	    
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
    handles.jobs(jobNumber).bodyname = filename(1:(end-(4+countNum)));
    bodyname = handles.jobs(jobNumber).bodyname;
     
    % Select the current project
    set(handles.GUI_st_job_lb, 'Value', jobNumber);
	
    % Create a list of files present in the image directory selected by the user
    dirList = dir(imagedirectory);
    dirList = struct2cell(dirList);
    dirList = dirList(1,:);
    
    % Find all files within this directory with the same name as the selected filename
    ind = strmatch(handles.jobs(jobNumber).bodyname, dirList);
    dirList = dirList(ind)';
    handles.jobs(jobNumber).lastimage = length(dirList);
      
    % Go to the directory where the images are stored
    if ~exist (imagedirectory, 'dir')
       h=errordlg(['Image directory ' imagedirectory ' does not exist. Please change path...']);
       uiwait(h);
       return
    else 
       cd (imagedirectory);
    end
   
    % Test whether the max greyvalue per frame is not more than the
    % bitdepth value specified on the gui
    % But do this only if it hasn't be done before (filesChecked matfile)
    if (exist ('filesChecked.mat', 'file') ~= 2)
       fprintf (1, 'Checking image files of job %d for correctness...\n', jobNumber);
    
       % Calculate the max posible grey value
       maxGreyValue = (2^handles.jobs(jobNumber).bitdepth);
    
       % Set the mouse pointer to busy
       set(gcf,'Pointer','watch');
    
       % Sort the images by successive numbers:
       % First we get all numbers and write them into a vector
       for jRearange = 1:length(dirList)
         tmpName = char(dirList(jRearange));
         if max (max (imread(tmpName))) > maxGreyValue
            % The frame contains a value higher than the bitdepth specified
            errormsg = ['Image file ' tmpName ' contains a grey value bigger than ' ...
                        num2str(maxGreyValue) ' (' num2str(max(max(imread(tmpName)))) ...
                        '). Please correct before loading job...'];
            h = errordlg (errormsg);
            uiwait (h);
            jobList(end) = [];
            if isempty (jobList)
              jobList = char('No project loaded');
              set(handles.GUI_st_job_lb,'Value',1);
            else
              set(handles.GUI_st_job_lb,'Value', length(jobList));
            end
            set(handles.GUI_st_job_lb,'String',jobList);
            handles.jobs(jobNumber) = [];
            guidata(hObject, handles);
            return
         else  
           % Add the job to the list
           imageNum(jRearange) = str2num(tmpName(length(handles.jobs(jobNumber).bodyname)+1:end-4));
         end
         
         % Create a file that is used to skip the test next time
         cd (imagedirectory);
         filesChecked = 1;
         save ('filesChecked.mat', 'filesChecked');
         
      end   % for jRearrange
    
      % Set the mouse pointer to normal again
      set(gcf,'Pointer','arrow');
    
      fprintf (1, 'All image files of job %d are correct!\n', jobNumber);
      
    else
      % Do the sorting without the greyvalue check
      for jRearange = 1:length(dirList)
        tmpName = char(dirList(jRearange));
        imageNum(jRearange) = str2num(tmpName(length(handles.jobs(jobNumber).bodyname)+1:end-4));
      end  % for jRearange
    end  % if exist ('filesChecked.mat', 'file')
    
    % Then we sort that vector and sort the dirList accordingly
    [junk,indVec] = sort(imageNum);
    handles.jobs(jobNumber).imagenameslist = dirList(indVec);
        
    % Create a directory to save the details and results of this job
    % Note: we call the directory results + bodyname + seq number
    % number will be lowest unoccupied number for this specific directory name
    cd (imagedirectory)
    done = 0;
    counter = 1;
    while done == 0
        newdirname = [];
        newdirname = ['results', bodyname, num2str(counter)];
           
        % Loop on until we find an unoccupied new dirname
        if exist (newdirname, 'dir') == 0
            mkdir (imagedirectory, newdirname);
            tempname = [imagedirectory, newdirname];
            mkdir(tempname,'body');
             
            handles.jobs(jobNumber).savedirectory = [imagedirectory, newdirname];
            done = 1;
        end
        counter = counter + 1;
    end
end

% Store the modified job list back in the GUI handle
set(handles.GUI_st_job_lb,'String',jobList);

% Update GUI handle struct
guidata(hObject, handles);

% Last but not least make sure the text field on the GUI show the latest values
ptFillFields(handles, handles.jobs(jobNumber))

%-------------------------------------------------------------------------------

% --- Executes on button press in GUI_st_deletejob_pb.
function GUI_st_deletejob_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_deletejob_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Initialize list
listnotfilled = [];

% Get the job list and the current job
jobList = get(handles.GUI_st_job_lb,'String');
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Joblist will only then be a cell, if there are jobs in it.
% Otherwise it is a string (No project loaded)
if ~iscell(jobList)
    % Show an error dialog with an appropriate message and wait
    % for the user to press a button
    h=errordlg('Sorry, there are no jobs to delete.');
    uiwait(h);
    return;
end

% Put a standard string in the job list window if there is no project to show
% and delete the currently selected job from the gui
if length(jobList) == 1
    jobList = char('No project loaded');
    listnotfilled = 1;
else
    jobList(jobNumber) = [];
end

% Set the list to the first project to be on the safe side
% And store new jobList in gui handle
set(handles.GUI_st_job_lb,'Value',1);
set(handles.GUI_st_job_lb,'String',jobList);

% Store job data
handles.jobs(jobNumber) = [];
guidata(hObject,handles);

% Show the data of the first job in the list, or if no job is present, 
% show the data of the defaultjob
if isempty (listnotfilled)
    ptFillFields (handles, handles.jobs(1))
else
    ptFillFields (handles, handles.defaultjob)
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

if ~exist(imagedir, 'file')
  h=errordlg('The selected image directory does not exist. Please select another directory.');
  uiwait(h);
  return
end

% Retrieve the current job number and assign it the image directory value
jobNumber = get(handles.GUI_st_job_lb,'Value');
handles.jobs(jobNumber).imagedirectory =  imagedir;

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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
% jobNumber = get(handles.GUI_st_job_lb,'Value');
% 
% handles.jobs(jobNumber).imagename =  imagename;
% 
% 
% % Update handles structure
% guidata(hObject, handles);
% 
% %%%%%%%%save altered values to disk%%%%%%%%%%%%
% cd (handles.jobs(jobNumber).savedirectory)
% jobvalues = handles.jobs(jobNumber);
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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(jobNumber).firstimage = [];
    ptFillFields(handles, handles.jobs(jobNumber))  % Revert the value back
    return
else
    handles.jobs(jobNumber).firstimage = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(jobNumber).lastimage = [];
    ptFillFields(handles, handles.jobs(jobNumber))  % Revert the value back
    return
else
    handles.jobs(jobNumber).lastimage = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(jobNumber).maxsearch = [];
    ptFillFields(handles, handles.jobs(jobNumber))  % Revert the value back
    return
else
    handles.jobs(jobNumber).maxsearch = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get the value and index from the timestepslide popup menu and assign to handles struct
val = get(hObject, 'Value');
list = get(hObject, 'String');
selected_val = list{val};
handles.jobs(jobNumber).timestepslide = str2double(selected_val);
handles.jobs(jobNumber).timestepslide_index = val;

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs (jobNumber).minsize = [];
    ptFillFields (handles, handles.jobs (jobNumber))  % Revert the value back
    return
else
    handles.jobs(jobNumber).mmpixel = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(jobNumber).minsize = [];
    ptFillFields(handles, handles.jobs(jobNumber))  % Revert the value back
    return
else
    handles.jobs(jobNumber).minsize = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(jobNumber).maxsize = [];
    ptFillFields(handles, handles.jobs(jobNumber))  % Revert the value back
    return
else
    handles.jobs(jobNumber).maxsize = val;
end

% Update handles structure
guidata(hObject, handles);


% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(jobNumber).minsdist = [];
    ptFillFields(handles, handles.jobs(jobNumber))  % Revert the value back
    return
else
    handles.jobs(jobNumber).minsdist = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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

ptSetCellValues (hObject,1);
 
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

ptSetCellValues (hObject,2);

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

ptSetCellValues (hObject,3);

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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(jobNumber).minedge = [];
    ptFillFields(handles, handles.jobs(jobNumber))  % Revert the value back
    return
else
    handles.jobs(jobNumber).minedge = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(jobNumber).noiseparameter = [];
    ptFillFields(handles, handles.jobs(jobNumber))  % Revert the value back
    return
else
    handles.jobs(jobNumber).noiseparameter = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(jobNumber).mincorrqualtempl = [];
    ptFillFields(handles, handles.jobs(jobNumber))  % Revert the value back
    return
else
    handles.jobs(jobNumber).mincorrqualtempl = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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
   jobNumber = get(handles.GUI_st_job_lb,'Value');

   imgDir = handles.jobs(jobNumber).imagedirectory;
   imgNam = handles.jobs(jobNumber).imagename;
   firstImg = handles.jobs(jobNumber).firstimage;
   lastImg = handles.jobs(jobNumber).lastimage;
   increment = handles.jobs(jobNumber).increment;
   savedirectory = handles.jobs(jobNumber).savedirectory;

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
   handles.jobs(jobNumber) = jobvalues;

   handles.jobs(jobNumber).imagedirectory = imgDir;
   handles.jobs(jobNumber).imagename = imgNam;
   handles.jobs(jobNumber).firstimage = firstImg;
   handles.jobs(jobNumber).lastimage = lastImg;
   handles.jobs(jobNumber).increment = increment;
   handles.jobs(jobNumber).savedirectory = savedirectory;

   ptFillFields(handles,handles.jobs(jobNumber))
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
   jobNumber = get(handles.GUI_st_job_lb,'Value');

   % Store the latest data in jobvalues.mat in the specified save directory
   if ~isempty(handles.jobs(jobNumber).savedirectory)
      cd (handles.jobs(jobNumber).savedirectory)
      jobvalues = handles.jobs(jobNumber);
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
   jobNumber = get(handles.GUI_st_job_lb,'Value');

   % Store the fields that shouldn't be defaulted in temp vars
   imgDir = handles.jobs(jobNumber).imagedirectory;
   imgNam = handles.jobs(jobNumber).imagename;
   firstImg = handles.jobs(jobNumber).firstimage;
   lastImg = handles.jobs(jobNumber).lastimage;
   increment = handles.jobs(jobNumber).increment;
   savedirectory = handles.jobs(jobNumber).savedirectory;

   % Set the default values to all fields
   handles.jobs(jobNumber) = handles.defaultjob;

   % Retrieve the values that we temporarily stored before
   handles.jobs(jobNumber).imagedirectory = imgDir;
   handles.jobs(jobNumber).imagename = imgNam;
   handles.jobs(jobNumber).firstimage = firstImg;
   handles.jobs(jobNumber).lastimage = lastImg;
   handles.jobs(jobNumber).increment = increment;
   handles.jobs(jobNumber).savedirectory = savedirectory;

   % Update all the fields on the gui
   ptFillFields(handles,handles.jobs(jobNumber))
end 

% Update handles structure
guidata(hObject, handles);

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
   % nrOfJobs is basically equal to the last entry nr in the job list which
   % again is equal to the length of the list
   nrOfJobs = length(jobList); 
else
   h=errordlg('A job should be loaded first, before a run can be started.');
   uiwait(h);
   return
end

% Loop through all the jobs in the joblist
for jobNumber = 1 : nrOfJobs

   % Get the image to image increment step
   Increment = handles.jobs(jobNumber).increment;

   % Get the first image number
   possibleImg = handles.jobs(jobNumber).firstimage;
   while (possibleImg + Increment) <= handles.jobs(jobNumber).lastimage
      possibleImg = possibleImg + Increment;
   end

   % Make sure the last image nr fits with the last image found in the previous loop
   % If these are not the same, adjust last image number
   if handles.jobs(jobNumber).lastimage > possibleImg
      handles.jobs(jobNumber).lastimage = possibleImg;
   end

   % Save the definite version of jobvalues
   if ~isempty(handles.jobs(jobNumber).savedirectory)
      cd (handles.jobs(jobNumber).savedirectory);
      jobvalues = handles.jobs(jobNumber);
      save ('jobvalues','jobvalues');
      clear jobvalues;
   end
        
   % Here's where the real tracking process starts for the selected job
   % AK: the try-catch should be uncommented as soon as testing is done!!!
   %try
      [M, clusterProps, cellProps, frameProps, imageCount] = ptTrackCells (handles.jobs(jobNumber), jobNumber);
   %catch    
      %fprintf (1, '\nJob number %d  had an error and could not be completed: %s\n', jobNumber, lasterr);
  
     % Save M, cluster and cell data
     %cd (handles.jobs(jobNumber).savedirectory);
     %save ('M','M');
     %save ('clusterProps','clusterProps');
     %save ('cellProps','cellProps');
   %end
   
   % Final message for the user to mark the end
   fprintf (1, '\nTracking finished...\n\n');
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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(jobNumber).sizetemplate = [];
    ptFillFields(handles, handles.jobs(jobNumber))  % Revert the value back
    return
else
    handles.jobs(jobNumber).sizetemplate = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

%-------------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_mintracklength_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_mintracklength_ed (see GCBO)
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

function GUI_mintracklength_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_mintracklength_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_mintracklength_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_mintracklength_ed as a double
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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(jobNumber).mintracklength = [];
    ptFillFields(handles, handles.jobs(jobNumber))  % Revert the value back
    return
else
    handles.jobs(jobNumber).mintracklength = val;
end

% Update handles structure
guidata(hObject, handles);


% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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

%if ~exist(savedirectory, 'dir')
%  mkdir (savedirectory);
%end

% Select the current job and store the directory name in the struct
% AK: Shouldn't some sort of validity check be done here??
jobNumber = get(handles.GUI_st_job_lb,'Value');
handles.jobs(jobNumber).savedirectory =  savedirectory;

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(savedirectory)
   cd (savedirectory)
   jobvalues = handles.jobs(jobNumber);
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
jobNumber = get(handles.GUI_st_job_lb,'Value');
handles.jobs(jobNumber).savedirectory =  savedirectory;

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(jobNumber).increment = [];
    ptFillFields(handles, handles.jobs(jobNumber))  % Revert the value back
    return
else
    handles.jobs(jobNumber).increment = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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

% Select the current job
jobNumber = get(handles.GUI_st_job_lb,'Value');

% If the joblist has entries get the number of entries else return,
% because there is really nothing to do
if iscell(jobList)
   nrOfJobs = length(jobList); 
else
   h=errordlg('At least one job should be loaded first, before test and initialization can be started.');
   uiwait(h);
   return
end 

% Start the test and initialization process
newCoord = ptInitializeJob (handles.jobs(jobNumber), jobNumber);

% Store the newly found coordinates in the handles struct
handles.jobs(jobNumber).coordinatespicone = newCoord;

% Update handles structure
guidata (hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(jobNumber).fi_halolevel = [];
    ptFillFields(handles, handles.jobs(jobNumber))  % Revert the value back
    return
else
    handles.jobs(jobNumber).fi_halolevel = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(jobNumber).la_halolevel = [];
    ptFillFields(handles, handles.jobs(jobNumber))  % Revert the value back
    return
else
    handles.jobs(jobNumber).la_halolevel = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get the value and index from the bitdepth  popup menu and assign to handles struct
val = get(hObject, 'Value');
list = get(hObject, 'String');
selected_val = list{val};
handles.jobs(jobNumber).bitdepth = str2double(selected_val);
handles.jobs(jobNumber).bitdepth_index = val;

% Calculate the maximal value of the image, depending on it's bitdepth and
% store this info in the handles struct
bitdepth = str2double(selected_val);
handles.jobs(jobNumber).intensityMax =  2^bitdepth - 1;

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Get number from the gui, convert it to a number and assign it to the handle;
% If it is not an number, throw and error dialog and revert to the old number
strval = get(hObject,'String');
val = str2double(strval);
if isnan (val)
    h = errordlg('Sorry, this field has to contain a number.');
    uiwait(h);          % Wait until the user presses the OK button
    handles.jobs(jobNumber).timeperframe = [];
    ptFillFields(handles, handles.jobs(jobNumber))  % Revert the value back
    return
else
    handles.jobs(jobNumber).timeperframe = val;
end

% Update handles structure
guidata(hObject, handles);

% Store the latest data in jobvalues.mat in the specified save directory
if ~isempty(handles.jobs(jobNumber).savedirectory)
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
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
jobNumber = get(handles.GUI_st_job_lb,'Value');

% Store the latest data in jobvalues.mat in the specified save directory
if exist ('handles.jobs(jobNumber).savedirectory')
   cd (handles.jobs(jobNumber).savedirectory)
   jobvalues = handles.jobs(jobNumber);
   save ('jobvalues','jobvalues')
   clear jobvalues
end

% Destroy the GUI
delete(handles.polyTrack_mainwindow);

%-------------------------------------------------------------------------------
