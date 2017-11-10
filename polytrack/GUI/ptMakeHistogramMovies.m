function varargout = ptMakeHistogramMovies(varargin)
%PTMAKEHISTOGRAMMOVIES M-file for ptMakeHistogramMovies.fig
%      PTMAKEHISTOGRAMMOVIES, by itself, creates a new PTMAKEHISTOGRAMMOVIES or raises the existing
%      singleton*.
%
%      H = PTMAKEHISTOGRAMMOVIES returns the handle to a new PTMAKEHISTOGRAMMOVIES or the handle to
%      the existing singleton*.
%
%      PTMAKEHISTOGRAMMOVIES('Property','Value',...) creates a new PTMAKEHISTOGRAMMOVIES using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ptMakeHistogramMovies_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      PTMAKEHISTOGRAMMOVIES('CALLBACK') and PTMAKEHISTOGRAMMOVIES('CALLBACK',hObject,...) call the
%      local function named CALLBACK in PTMAKEHISTOGRAMMOVIES.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ptMakeHistogramMovies

% Last Modified by GUIDE v2.5 20-Aug-2004 17:24:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ptMakeHistogramMovies_OpeningFcn, ...
                   'gui_OutputFcn',  @ptMakeHistogramMovies_OutputFcn, ...
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

%------------------------------------------------------------------------

% --- Executes just before ptMakeHistogramMovies is made visible.
function ptMakeHistogramMovies_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for ptMakeHistogramMovies
handles.output = hObject;

% Assign the default job values to the GUI handle so it can be passed around
home = getenv('HOME');
if isempty (home)
   if ispc
      home = 'H:';
   else
      home = '/tmp';
   end
   fprintf (1, 'HOME environment variable not set. Setting default: %s\n', home);
end

handles.savedirectory = [home filesep 'hist_movie.mov'];
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

%------------------------------------------------------------------------

% --- Outputs from this function are returned to the command line.
function varargout = ptMakeHistogramMovies_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%------------------------------------------------------------------------

% --- Executes on button press in GUI_st_add_hist_pb.
function GUI_st_add_hist_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_add_hist_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Select a file called '<something>.mat' from a user selected directory
[filename, directory] = uigetfile ({'*.mat', 'MAT-Files'; '*.*', 'all files'}, ...
                                      'Please select the first mat file with histogram data (of a range)');

% Do nothing in case the user doesn't select anything
if filename == 0
   return
else
   % Convert filename to lowercase
   fileLower = lower (filename);
   
   % Check whether it is a CSV file (ext .csv)
   if (strfind (fileLower, '.mat')) == []
      errormsg = ['File ' filename ' is not a MAT file. Please choose another file.'];
      h = errordlg (errormsg);
      uiwait (h);
      return
   end
         
   % Get the current file list
   fileList = get (handles.GUI_histogram_file_lb, 'String');
   
   % Cat together the directory and filename
   filePath = [directory filename];
      
   % If fileList consists of more than one entry, append the new file path after the
   % last one else just put it in the fileList as the first entry
   if ~iscell (fileList)
      fileList = cellstr (filePath); 
   else   % The list already had some files in it
       
      % Read the mat file: returns a struct of which the fieldnames are not known
      matFile = load (filePath);
      
      % Get the field out; if more than 1: stop
      field = fieldnames (matFile);
      if size (field,1) > 1
         errormsg = ['File ' filename ' contains more than a single histogram vector. Please choose another file.'];
         h = errordlg (errormsg);
         uiwait (h);
         return
      end
      eval (['histCur=matFile.', char(field), ';']);
       
      % Read the previous mat file
      matFilePrev = load (fileList{end});
      fieldPrev = fieldnames (matFilePrev);
      eval (['histPrev=matFilePrev.', char(fieldPrev), ';']);
      
      % Make sure the new mat file has the same dimensions as the one read before
      if length (histCur) ~= length (histPrev) 
         errormsg = ['Histogram ' filename ' has another dimension as the ones already loaded. Please choose another file.'];
         h = errordlg (errormsg);
         uiwait (h);
         return
      end
      
      % Add the filename to the list
      fileList(end+1) = cellstr (filePath);
   end
       
   % Go to the selected directory: user comfort for next file
   cd (directory);
            
   % Clear the matFile matrix for the next entry
   clear matFile; clear matFilePrev; clear histCur; clear histPrev;
end
      
% Store the modified job list back in the GUI handle
set(handles.GUI_histogram_file_lb,'String',fileList);

% Update GUI handle struct
guidata(hObject, handles);

%------------------------------------------------------------------------

% --- Executes on button press in GUI_st_remove_hist_pb.
function GUI_st_remove_hist_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_st_remove_hist_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Get the job list and the current job
fileList = get(handles.GUI_histogram_file_lb,'String');
fileNumber = get(handles.GUI_histogram_file_lb,'Value');

% fileList will only then be a cell, if there are results in it.
% Otherwise it is a string (No csv's loaded)
if ~iscell(fileList)
    % Show an error dialog with an appropriate message and wait
    % for the user to press a button
    h=errordlg('Sorry, there are no MAT files to delete.');
    uiwait(h);
    return;
end

% Put a standard string in the file list window if there are no files to show
% and delete the currently selected file from the gui
if length(fileList) == 1
    fileList = char('No histogram files loaded');
else
    fileList(fileNumber) = [];
end

% Set the list to the first project to be on the safe side
% And store new jobList in gui handle
set (handles.GUI_histogram_file_lb, 'Value', 1);
set (handles.GUI_histogram_file_lb, 'String', fileList);

% Update GUI handles struct
guidata (hObject,handles);

%------------------------------------------------------------------------

% --- Executes on button press in GUI_generate_movie_pb.
function GUI_generate_movie_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_generate_movie_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

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

% Cd into the directory and init the movie file
cd (pathString);
movie = avifile ([filename ext]);
%makeQTMovie ('start', [filename ext]);

% retrieve the histogram file list
fileList = get(handles.GUI_histogram_file_lb,'String');

% fileList will only then be a cell, if there are entries in it.
% Otherwise it is a string (No mat's loaded)
if iscell(fileList)
    nrOfFiles = length(fileList);
else
    % Show an error dialog with an appropriate message and wait
    % for the user to press a button
    h=errordlg ('Sorry, cannot start movie generation since no histogram files are loaded.');
    uiwait (h);
    return
end

% Set the mouse pointer to busy
set(gcbf,'Pointer','watch');

filePath = {};
basename = {};
numberOfFiles = [];

% Loop through all the jobs in the joblist
for fileNumber = 1 : nrOfFiles
    
   % Get the filename and break it in pieces
   filePathAndName = char(fileList(fileNumber));
   [pathString, filename, ext, version] = fileparts (filePathAndName);
   
   % Store the path for later use
   filePath(fileNumber) = cellstr(pathString);
   
   % Find the basename (filename without the digits)
   number = 0;
   countNum = 0;
   while ~isnan(number) & (countNum <3)
      countNum = countNum+1;
      number = str2num(filename((end-countNum):end));
   end
   basename(fileNumber) = cellstr(filename(1:((end-countNum))));

   % Create a list of files present in the image directory selected by the user
   dirList = dir(pathString);
   dirList = struct2cell(dirList);
   dirList = dirList(1,:);
    
   % Find all files within this directory with the same name as the selected filename
   ourFiles = strmatch(char(basename(fileNumber)), dirList);
   dirList = dirList(ourFiles)';
   numberOfFiles(fileNumber) = length(dirList);
   
end
    
% Load one file to see what the binsize is
matFile = load (filePathAndName);
fieldCur = fieldnames (matFile);
eval (['histCur=matFile.', char(fieldCur), ';']);
binSize = length(histCur);

% Initialize the histogram summary vector
histSum = zeros(numberOfFiles(1),binSize);

% Now that we know the names we're looking for and the number of files we
% can start looping through the hist files and find the max value (for the
% figures we show later on
for iCount = 1 : numberOfFiles(1)
    
   % Initialize some variables
   filename = {};
    
   % Get the filenames we need
   for jCount = 1 : nrOfFiles
      filename(jCount) = cellstr([char(filePath(jCount)) filesep char(basename(jCount)) num2str(iCount) '.mat']);
      
      % Read the histogram mat file
      matFile = load (char(filename(jCount)));
      fieldCur = fieldnames (matFile);
      eval (['histCur=matFile.', char(fieldCur), ';']);
      
      if jCount == 1
         histSum(iCount,:) = histCur;
      else
         histSum(iCount,:) = histSum(iCount,:) + histCur;
      end
   end
      
   % Calculate the maximum velocity value
   yMax(iCount) = max (histSum(iCount,:));
end

% Calculate the max y-axis value
ymax = max(yMax);

% Get the movie title
movieTitle = get (handles.GUI_movie_title_ed, 'String');

% Let the user know what we're doing
fprintf (1, 'Frame: ');

for movieStep = 1 : numberOfFiles(1)
    
   % Let the user know where we are
   fprintf (1, '%d ', movieStep);
      
   % Show the figure
   fig = figure ('Position', [200 300 500 400]); bar (histSum(movieStep,:)); title([movieTitle]);
    
   % Set the axis
   axis ([0 binSize+1 0 ymax])
    
   % Get a handle to the figure frame
   F = getframe (fig);
      
   % Add the frame to the movie
   movie = addframe (movie, F);
   %makeQTMovie ('addaxes', gca);
      
   % Close the figure
   close; 
end

% Close the movie file
movie = close(movie);
%makeQTMovie ('finish');

% Set the mouse pointer to normal again
set(gcbf,'Pointer','arrow');

% Tell the user that we have finished
message = ['The movie has been written into ' saveDirectory];
uiwait (msgbox (message, 'Movie generation finished', 'modal'));

% Update GUI handles struct
guidata (hObject,handles);

%------------------------------------------------------------------------

% --- Executes on button press in GUI_browse_savepath_pb.
function GUI_browse_savepath_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_browse_savepath_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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

% --- Executes on selection change in GUI_histogram_file_lb.
function GUI_histogram_file_lb_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_histogram_file_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GUI_histogram_file_lb contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GUI_histogram_file_lb
handles = guidata(hObject);

% Get the number of the currently selected histogram in the list
fileNumber = get (hObject,'Value');

%------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_histogram_file_lb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_histogram_file_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%------------------------------------------------------------------------

function GUI_movie_title_ed_Callback(hObject, eventdata, handles)
% hObject    handle to GUI_movie_title_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GUI_movie_title_ed as text
%        str2double(get(hObject,'String')) returns contents of GUI_movie_title_ed as a double
handles = guidata(hObject);

% Get the directory path and name
movieTitle = get (hObject,'String');

% Update the edit field
set (handles.GUI_movie_title_ed, 'String', movieTitle);

% Select the current job and store the directory name in the struct
handles.movieTitle =  saveDirectory;

% Update GUI handles struct
guidata (hObject,handles);

%------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function GUI_movie_title_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_movie_title_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


