function varargout = ptAverageData(varargin)
%PTAVERAGEDATA M-file for ptAverageData.fig
%      PTAVERAGEDATA, by itself, creates a new PTAVERAGEDATA or raises the existing
%      singleton*.
%
%      H = PTAVERAGEDATA returns the handle to a new PTAVERAGEDATA or the handle to
%      the existing singleton*.
%
%      PTAVERAGEDATA('Property','Value',...) creates a new PTAVERAGEDATA using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ptAverageData_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      PTAVERAGEDATA('CALLBACK') and PTAVERAGEDATA('CALLBACK',hObject,...) call the
%      local function named CALLBACK in PTAVERAGEDATA.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES


% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ptAverageData_OpeningFcn, ...
                   'gui_OutputFcn',  @ptAverageData_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before ptAverageData is made visible.
function ptAverageData_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for ptAverageData
handles.output = hObject;

% Create a default job structure that defines all the fields used in the program
%data = struct('savedirectory', []);

% Set the HOME env variable if not set already
home = getenv('HOME');
if isempty (home)
   if ispc
      home = 'H:';
   else
      home = '/tmp';
   end
   fprintf (1, 'HOME environment variable not set. Setting default: %s\n', home);
end

handles.savedirectory = [home filesep 'polytrack_average.csv'];
set (handles.GUI_savepath_ed, 'String', handles.savedirectory);

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

% --- Outputs from this function are returned to the command line.
function varargout = ptAverageData_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%--------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_st_result_lb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_st_result_lb (see GCBO)
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
function GUI_st_result_lb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_result_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GUI_st_job_lb contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GUI_st_job_lb
handles = guidata(hObject);

% Get the number of the currently selected project in the list
fileNumber = get (hObject,'Value');

%-------------------------------------------------------------------------------

% --- Executes on button press in GUI_st_addjob_pb.
function GUI_st_add_result_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_add_result_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Make sure the handles struct can be used in this function
handles = guidata(hObject);

% Assign the radiobutton value to the handles struct
alwaysReadCsvFile = get (handles.GUI_always_read_csv_cb, 'Value');

% Select an csv file from a user selected directory
[filename, directory] = uigetfile ({'*.csv', 'CSV-Files'; '*.*', 'all files'}, ...
                                      'Please select a CSV file with Polytrack Results');

% Do nothing in case the user doesn't select anything
if filename == 0
   return
else
   % Convert filename to lowercase
   fileLower = lower (filename);
   
   % Check whether it is a CSV file (ext .csv)
   if (strfind (fileLower, '.csv')) == []
      errormsg = ['File ' filename ' is not a CSV file. Please choose another file.'];
      h = errordlg (errormsg);
      uiwait (h);
      return
   end
         
   % Get the current file list
   fileList = get (handles.GUI_st_result_lb, 'String');
   
   % Cat together the directory and filename
   filePath = [directory filename];
      
   % If fileList consists of more than one entry, append the new file path after the
   % last one else just put it in the fileList as the first entry
   if ~iscell (fileList)
      fileList = cellstr (filePath); 
   else   % The list already had some files in it
       
      % Skip dimensions and x-axis test when user requests this
      if ~alwaysReadCsvFile
         % Read the csv file
         csvFile = csvread (filePath);
       
         % Read the previous csv file
         csvFilePrev = csvread (fileList{end});
      
         % Make sure the new csv file has the same dimensions as the one read before
         if size (csvFile,1) ~= size (csvFilePrev,1) | size (csvFile,2) ~= size (csvFilePrev,2)
            errormsg = ['File ' filename ' has another dimension as the files already loaded. Please choose another file.'];
            h = errordlg (errormsg);
            uiwait (h);
            return
         end
      
         
         % Also make sure their X-axis are the same (unless the user requests not to)
         if csvFile(1,:) ~= csvFilePrev(1,:)
            errormsg = ['File ' filename ' has a different X-axis structure as the files loaded before. Please choose another file.'];
            h = errordlg (errormsg);
            uiwait (h);
            return
         end
      end
      
      % Add the filename to the list
      fileList(end+1) = cellstr (filePath);
   end
       
   % Go to the selected directory: user comfort for next file
   cd (directory);
            
   % Clear the csvFile matrix for the next entry
   clear csvFile; clear csvFilePrev;
end
      
% Store the modified job list back in the GUI handle
set(handles.GUI_st_result_lb,'String',fileList);

% Update GUI handle struct
guidata(hObject, handles);

%------------------------------------------------------------------------

% --- Executes on button press in GUI_st_remove_result_pb.
function GUI_st_remove_result_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_remove_result_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Get the job list and the current job
fileList = get(handles.GUI_st_result_lb,'String');
fileNumber = get(handles.GUI_st_result_lb,'Value');

% fileList will only then be a cell, if there are results in it.
% Otherwise it is a string (No csv's loaded)
if ~iscell(fileList)
    % Show an error dialog with an appropriate message and wait
    % for the user to press a button
    h=errordlg('Sorry, there are no CSV files to delete.');
    uiwait(h);
    return;
end

% Put a standard string in the file list window if there are no files to show
% and delete the currently selected file from the gui
if length(fileList) == 1
    fileList = char('No CSV files loaded');
else
    fileList(fileNumber) = [];
end

% Set the list to the first project to be on the safe side
% And store new jobList in gui handle
set (handles.GUI_st_result_lb, 'Value', 1);
set (handles.GUI_st_result_lb, 'String', fileList);

% Update GUI handles struct
guidata (hObject,handles);

%-------------------------------------------------------------------------------

% --- Executes on button press in GUI_st_average_results_pb.
function GUI_st_average_results_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_average_results_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Assign the radiobutton value to the handles struct
alwaysReadCsvFile = get (handles.GUI_always_read_csv_cb, 'Value');

% Retrieve the directory and filename where to save the result
saveDirectory = handles.savedirectory;
[pathString, filename, ext, version] = fileparts (saveDirectory);

if ~exist (pathString, 'dir')
   % Show an error dialog with an appropriate message and wait
   % for the user to press a button
   h=errordlg ('The save directory does not exist. Please choose another directory.');
   uiwait (h);
   return
end

% Cd into the directory
cd (pathString);

% retrieve the csv file list
fileList = get(handles.GUI_st_result_lb,'String');

% fileList will only then be a cell, if there are results in it.
% Otherwise it is a string (No csv's loaded)
if ~iscell(fileList)
    % Show an error dialog with an appropriate message and wait
    % for the user to press a button
    h=errordlg ('Sorry, there are no CSV files to average.');
    uiwait (h);
    return
end

% Prepare storage for the values to be averaged
% valuesToAverage = zeros (length(fileList), size (csvread(fileList{1}),2)+1);
clear valuesToAverage;

% Set the mouse pointer to busy
set(gcbf,'Pointer','watch');

% Read all the csv files and store in matrix
for iCount = 1 : length (fileList)
   csvFile = csvread (fileList{iCount});
   
   % In case csv files where selected without limit checking, we have to work
   % with the shortest one and set the other ones to NaN at the end of the
   % vector
   if alwaysReadCsvFile
      % Find the shortest vector (the -1 after 'size' is necessary for the
      % rank we added to valuesToAverage)
      if iCount > 1
         if size(csvFile,2) > (size(valuesToAverage(iCount-1,:),2)-1)
            % Shorten the current one
            overshoot = size(csvFile,2) - (size(valuesToAverage(iCount-1,:),2)-1);
            csvFile (end-overshoot:end) = NaN;
            valuesToAverage(:, end+1:(size(csvFile,2)+1)) = NaN;
         else
            % Shorten the ones already stored
            overshoot = (size(valuesToAverage(iCount-1,:),2)-1) - size(csvFile,2);
            csvFile (:,end+1:(size(valuesToAverage(iCount-1,:),2)-1)) = NaN;
            valuesToAverage(:, end-(overshoot-1):end) = NaN;
         end
      end
   end
   
   % Add the csv row to the total
   valuesToAverage(iCount,:) = [iCount csvFile(2,:)];
end

% Get the X-axis values from the last csv file read (any csv file is okay
% since these are always the  same)
xAxisValues = [NaN csvFile(1,:)];

% Average the values
averageValues = sum (valuesToAverage) / size (valuesToAverage,1);
averageValues(1) = NaN;   % First one shoulden't be averaged since it is only a sequence number

% Put it all back together
averageForCsv = [xAxisValues ; valuesToAverage ; averageValues];

% Write the csv file to the specified save file
csvwrite (saveDirectory, averageForCsv);

% Set the mouse pointer to normal again
set(gcbf,'Pointer','arrow');

% Tell the user that we have finished
message = ['The average files have been written into ' saveDirectory];
uiwait (msgbox (message, 'Averaging finished', 'modal'));

% Update the handles structure
guidata (hObject,handles);

%------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_savepath_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_savepath_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%------------------------------------------------------------------------

function GUI_savepath_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_savepath_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_savepath_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_savepath_ed as a double
handles = guidata(hObject);

% Get the directory path and name
saveDirectory = get (hObject,'String');

% Get the directory part of the string
[pathString, filename, ext, version] = fileparts (saveDirectory);

if ~exist (pathString, 'dir')
   % Ask the user whether to create it
   question = ['The directory ' pathString ' does not exist. Shall I create it for you?'];
   answer = questdlg (question);
   switch answer 
      case 'Yes',
         % Create the directory
         mkdir (pathString);
      case 'No',
         % Put it back to what it was
         saveDirectory = handles.savedirectory;
      case 'Cancel',
         % Put it back to what it was
         saveDirectory = handles.savedirectory;   
   end  % switch
end

% Update the edit field
set (handles.GUI_savepath_ed, 'String', saveDirectory);

% Select the current job and store the directory name in the struct
handles.savedirectory =  saveDirectory;

% Update GUI handles struct
guidata (hObject,handles);

%--------------------------------------------------------------------------

function File_menu_Callback(hObject, eventdata, handles)
% hObject    handle to File_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This is just the filemenu which doesn't do anything really; just keep the
% submenus like 'Exit'

%--------------------------------------------------------------------------

function exit_menuitem_Callback(hObject, eventdata, handles)
% hObject    handle to exit_menuitem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Destroy the GUI
delete(handles.ptAverageData_mainwindow);

%--------------------------------------------------------------------------

% --- Executes on button press in GUI_always_read_csv_cb.
function GUI_always_read_csv_cb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_always_read_csv_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of GUI_always_read_csv_cb

% Nothing to be done here

%--------------------------------------------------------------------------

% --- Executes on button press in GUI_browse_savepath_pb.
function GUI_browse_savepath_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_browse_savepath_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Make sure the handles struct can be used in this function
handles = guidata(hObject);

% Get the directory path and name
saveDirectory = handles.savedirectory;

% Get the directory part of the string
[pathString, filename, ext, version] = fileparts (saveDirectory);

% Select a directory where the average csv file will be saved
directory = uigetdir (pathString, 'Select Directory');

% Do nothing in case the user doesn't select anything
if directory == 0
   return
else
   % this will be the new save directory
   newDirectory = [directory filesep filename ext];
   
   % Update the edit field
   set (handles.GUI_savepath_ed, 'String', newDirectory);

   % Select the current job and store the directory name in the struct
   handles.savedirectory =  newDirectory;
end

% Update GUI handles struct
guidata (hObject,handles);

