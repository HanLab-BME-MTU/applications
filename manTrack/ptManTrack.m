function varargout = ptManTrack(varargin)
% PTMANTRACK M-file for ptManTrack.fig
%      PTMANTRACK, by itself, creates a new PTMANTRACK or raises the existing
%      singleton*.
%
%      H = PTMANTRACK returns the handle to a new PTMANTRACK or the handle to
%      the existing singleton*.
%
%      PTMANTRACK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PTMANTRACK.M with the given input arguments.
%
%      PTMANTRACK('Property','Value',...) creates a new PTMANTRACK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ptManTrack_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ptManTrack_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help ptManTrack

% Last Modified by GUIDE v2.5 18-Feb-2005 12:35:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ptManTrack_OpeningFcn, ...
                   'gui_OutputFcn',  @ptManTrack_OutputFcn, ...
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


% --- Executes just before ptManTrack is made visible.
function ptManTrack_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ptManTrack (see VARARGIN)

% Choose default command line output for ptManTrack
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ptManTrack wait for user response (see UIRESUME)
% uiwait(handles.ptManTrack);


% --- Outputs from this function are returned to the command line.
function varargout = ptManTrack_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pb_selectmovie_ptmt.
function pb_selectmovie_ptmt_Callback(hObject, eventdata, handles)
% hObject    handle to pb_selectmovie_ptmt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% Let the user browse for an image file and path
movieDirectory = uigetdir('','Select movie directory');

if movieDirectory == 0
    return;
end

% And store this in the handles struct
if ~exist(movieDirectory, 'dir')
   h = errordlg('This directory does not exist. Please select another directory.');
   uiwait(h);          % Wait until the user presses the OK button
   return;
end

% Update handles
handles.movieDirectory = movieDirectory;
set (handles.edit_moviepath_ptmt, 'String', movieDirectory);

% There should be at least one result dir available to get the coordinates
% from (MPM). If more are available ask the user
dirList = dir(movieDirectory);
dirList = struct2cell(dirList);
dirList = dirList(1,:);
ind = strmatch('results_', dirList);
dirList = dirList(ind)';

if isempty(dirList)
    msgStr = ['No results directory available. Please choose an image set that has been analyzed before.'];
    h = errordlg(msgStr);
    uiwait(h);
    return;
elseif length(dirList) == 1
    handles.resultsDirectory = char(dirList(1));
else
    % Let the user browse for a result directory
    resultsDirectory = uigetdir(movieDirectory,'Select results directory');
    
    if resultsDirectory == 0
        return;
    elseif isempty(findstr(resultsDirectory,'results'))
        msgStr = ['This is not a results directory. Exiting...'];
        h = errordlg(msgStr);
        uiwait(h);
        return;
    end  
    
    handles.resultsDirectory = resultsDirectory;
end

% Check for images in the directory provided
dirList = dir(movieDirectory);
dirList = struct2cell(dirList);
dirList = dirList(1,:);

% If we found some tiff's ask the user to select a specific file
if ~isempty(dirList)
    curDir = pwd; cd(directory);
    [imageFile, imageDirectory] = uigetfile('*.tif','TIFF files','Select an image from the desired sequence');
    cd(curDir);
else
    msg = 'No TIFF files can be found in this directory. Exiting...';
    h = errordlg(msg);
    uiwait(h);
    return;
end

if imageFile == 0
    % Nothing selected
    msg = 'No TIFF file has been selected. Exiting...';
    h = errordlg(msg);
    uiwait(h);
    return;
end

% Separate into parts
%[directory, filename, extension, dummy] = fileparts([imageDirectory filesep imageFile]);

% Find the body part of the filename
number = 0;
countNum = 0;
while ~isnan(number) & (countNum < 3)   
    countNum = countNum + 1;
    number = str2num(imageFile(end-(4+countNum):end-4));
end
bodyname = lower(imageFile(1:(end-(4+countNum))));

% Filter the directory list previously acquired
ind = strmatch(bodyname, dirList);
dirList = dirList(ind)';

% Store all these values
imageName      = bodyname;
firstImage     = 1;
lastImage      = length(dirList);
imageNameList  = dirList;
intensityMax   = 4095;

% Also store in the handles struct
handles.imageName      = imageName;
handles.firstImage     = firstImage;
handles.lastImage      = lastImage;
handles.imageNameList  = imageNameList;
handles.intensityMax   = intensityMax;

% Update handles structure
guidata(hObject, handles);


function edit_moviepath_ptmt_Callback(hObject, eventdata, handles)
% hObject    handle to edit_moviepath_ptmt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_moviepath_ptmt as text
%        str2double(get(hObject,'String')) returns contents of edit_moviepath_ptmt as a double
% Get input from the gui and assign it to the handle; if directory does not
% exist, ask user whether to create it
movieDirectory = get(hObject,'String');

% If a filesep exist at the end of the path, remove it
if movieDirectory(end) == '/'
    movieDirectory(end) = '';
end

% If the path doesn't exist ask user if it should be created
if ~exist(tmpDirectory, 'dir')
   msgStr = ['This directory does not exist. Please choose another directory.'];
   h = errordlg(msgStr);
   uiwait(h);
   return;
end

% Assign it to the guiData structure
handles.movieDirectory = movieDirectory;

% There should be at least one result dir available to get the coordinates
% from (MPM). If more are available ask the user
dirList = dir(directory);
dirList = struct2cell(dirList);
dirList = dirList(1,:);
ind = strmatch('results_', dirList);
dirList = dirList(ind)';

if isempty(dirList)
    msgStr = ['No results directory available. Please choose an image set that has been analyzed before.'];
    h = errordlg(msgStr);
    uiwait(h);
    return;
elseif length(dirList) == 1
    handles.resultsDirectory = char(dirList(1));
else
    hFig = figure;
    hBox = uicontrol(hFig, 'Style','listbox','String', dirList, 'Value', 1, 'Callback', 'selectResultDir');
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_moviepath_ptmt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_moviepath_ptmt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_starttrack_ptmt.
function pb_starttrack_ptmt_Callback(hObject, eventdata, handles)
% hObject    handle to pb_starttrack_ptmt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

if isfield(handles,'movieDirectory');
    % Show the slider window
    ptShowImageSequence(handles.movieDirectory);
else
    msgStr = ['No movie directory selected. Please choose a directory first'];
    h = errordlg(msgStr);
    uiwait(h);
    return;
end

% Update handles structure
guidata(hObject, handles);

