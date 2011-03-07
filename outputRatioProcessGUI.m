function varargout = outputRatioProcessGUI(varargin)
% OUTPUTRATIOPROCESSGUI M-file for outputRatioProcessGUI.fig
%      OUTPUTRATIOPROCESSGUI, by itself, creates a new OUTPUTRATIOPROCESSGUI or raises the existing
%      singleton*.
%
%      H = OUTPUTRATIOPROCESSGUI returns the handle to a new OUTPUTRATIOPROCESSGUI or the handle to
%      the existing singleton*.
%
%      OUTPUTRATIOPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OUTPUTRATIOPROCESSGUI.M with the given input arguments.
%
%      OUTPUTRATIOPROCESSGUI('Property','Value',...) creates a new OUTPUTRATIOPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before outputRatioProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to outputRatioProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help outputRatioProcessGUI

% Last Modified by GUIDE v2.5 24-Aug-2010 11:18:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @outputRatioProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @outputRatioProcessGUI_OutputFcn, ...
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


% --- Executes just before outputRatioProcessGUI is made visible.
function outputRatioProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% Available tools 
% UserData data:
%       userData.mainFig - handle of main figure
%       userData.handles_main - 'handles' of main figure
%       userData.procID - The ID of process in the current package
%       userData.crtProc - handle of current process
%       userData.crtPackage - handles of current package
%       userData.procConstr - constructor of current process
%
%       userData.questIconData - help icon image information
%       userData.colormap - color map information
%

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright)

userData = get(handles.figure1, 'UserData');
% Choose default command line output for segmentationProcessGUI
handles.output = hObject;

% Get main figure handle and process id
t = find(strcmp(varargin,'mainFig'));
userData.mainFig = varargin{t+1};
userData.procID = varargin{t+2};
userData.handles_main = guidata(userData.mainFig);

% Notify main figure that the current process setting panel is open
% setFlag = getappdata(userData.mainFig,'setFlag');
% setFlag(userData.procID) = 1;
% setappdata(userData.mainFig,'setFlag',setFlag);

% Get current package and process
userData_main = get(userData.mainFig, 'UserData');
userData.crtPackage = userData_main.crtPackage;
userData.crtProc = userData.crtPackage.processes_{userData.procID};

% Get current process constructer
eval ( [ 'userData.procConstr = @', ...
    userData.crtPackage.processClassNames_{userData.procID},';']);

% If process does not exist, create a default one in user data.
if isempty(userData.crtProc)
    userData.crtProc = userData.procConstr(userData_main.MD(userData_main.id), ...
                                userData.crtPackage.outputDirectory_);
end

% Get icon infomation
userData.questIconData = userData_main.questIconData;
userData.colormap = userData_main.colormap;

% ---------------------- Channel Setup -------------------------

funParams = userData.crtProc.funParams_;

% Set up available input channels
set(handles.listbox_input1, 'String', {userData_main.MD(userData_main.id).channels_.channelPath_},...
        'Userdata', 1: length(userData_main.MD(userData_main.id).channels_));
    
% Set up input channel (one channel)
if ~isempty(funParams.ChannelIndex)
    set(handles.edit_dir, 'String', ...
            {userData_main.MD(userData_main.id).channels_(funParams.ChannelIndex).channelPath_}, ...
            'Userdata', funParams.ChannelIndex )
else
    
    temp = cellfun(@(x)isa(x, 'RatioProcess'), userData.crtPackage.processes_);
    
    % If ratio process exist and has a numerator
    if any(temp) && ~isempty(userData.crtPackage.processes_{temp}.funParams_.ChannelIndex)
        nu = userData.crtPackage.processes_{temp}.funParams_.ChannelIndex(1);
        
        set(handles.edit_dir, 'String', ...
            {userData_main.MD(userData_main.id).channels_(nu).channelPath_}, ...
            'Userdata', nu )
        set(handles.listbox_input1, 'Value', nu)        
    end    
end
    
if ~isempty(funParams.ChannelIndex)
    set(handles.listbox_input1, 'Value', funParams.ChannelIndex(1))
end
    
% ---------------------- Parameter Setup -------------------------

set(handles.edit_factor, 'String', num2str(funParams.ScaleFactor))
set(handles.edit_path, 'String', funParams.OutputDirectory)

% ----------------------Set up help icon------------------------

% Set up help icon
set(hObject,'colormap',userData.colormap);
% Set up package help. Package icon is tagged as '0'
axes(handles.axes_help);
Img = image(userData.questIconData); 
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off','YDir','reverse');
set(Img,'ButtonDownFcn',@icon_ButtonDownFcn);

if openHelpFile
    set(Img, 'UserData', struct('class',class(userData.crtProc)))
end

% ----------------------------------------------------------------

% Update user data and GUI data
set(userData.mainFig, 'UserData', userData_main);
set(hObject, 'UserData', userData);

uicontrol(handles.pushbutton_done);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = outputRatioProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% Call back function of 'Apply' button
userData = get(handles.figure1, 'UserData');
userData_main = get(userData.mainFig, 'UserData');

% -------- Check user input --------
channelIndex = get (handles.edit_dir, 'Userdata');
outputDir = get(handles.edit_path, 'String');

if isempty(channelIndex)
   errordlg('Please select a channel as ratio channel.','Setting Error','modal') 
    return;
end
if isempty(outputDir)
   errordlg('Please select a path to save .tiff file output.','Setting Error','modal') 
    return;
end
if isnan(str2double(get(handles.edit_factor, 'String'))) ...
                || str2double(get(handles.edit_factor, 'String')) < 0
    errordlg('Please provide a valid input for ''Scale Factor''.','Setting Error','modal');
    return;
end

% -------- Process Sanity check --------
% ( only check underlying data )

try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

%---------Check if channel indexs are changed---------

funParams = userData.crtProc.funParams_;

if isempty(funParams.ChannelIndex) || funParams.ChannelIndex ~= channelIndex

    % If channel indexs are changed, set procChanged to true
    userData.crtProc.setProcChanged(true);
end

% -------- Set parameter --------


    
    % Get parameter
    
    funParams.ChannelIndex = channelIndex;
    funParams.OutputDirectory = outputDir;    
    funParams.ScaleFactor = str2double(get(handles.edit_factor, 'String'));
    
    % Set parameters
    userData.crtProc.setPara(funParams);



% --------------------------------------------------


% If this is a brand new process, attach current process to MovieData and 
% package's process list 
if isempty( userData.crtPackage.processes_{userData.procID} )
    
    % Add new process to both process lists of MovieData and current package
    userData_main.MD(userData_main.id).addProcess( userData.crtProc );
    userData.crtPackage.setProcess(userData.procID, userData.crtProc);
    
    % Set font weight of process name bold
    eval([ 'set(userData.handles_main.checkbox_',...
            num2str(userData.procID),', ''FontWeight'',''bold'')' ]);
end

% ----------------------Sanity Check (II, III check)----------------------

% Do sanity check - only check changed parameters
procEx = userData.crtPackage.sanityCheck(false,'all');

% Return user data !!!
set(userData.mainFig, 'UserData', userData_main)

% Draw some bugs on the wall 
for i = 1: length(procEx)
   if ~isempty(procEx{i})
       % Draw warning label on the i th process
       userfcn_drawIcon(userData.handles_main,'warn',i,procEx{i}(1).message, true); % user data is retrieved, updated and submitted
   end
end

% Refresh user data !!
userData_main = get(userData.mainFig, 'UserData');

% -------------------- Apply setting to all movies ------------------------

if get(handles.checkbox_applytoall, 'Value')
errorIndex = [];

for x = 1: length(userData_main.MD)
    
   if x == userData_main.id
      continue 
   end
   
   % Customize funParams to other movies 
   % ChannelIndex - all channels
   % bleedChannelIndex - the indexs of bleedthrough channels
   % bleedCoefficients - the coefficients of bleedthrough channels
   % OutputDirectory - pacakge output directory
   l = length(userData_main.MD(x).channels_);
   
   % If channel index is larger than the number of channels of current movie, report error 
   if channelIndex > l
       errorIndex = cat(2, errorIndex, x);
       continue
   else
   
       funParams.ChannelIndex = channelIndex;
   end
   
   % If process exists, don't change, if does not exist, set output dir to
   % default directory
   if isempty(userData_main.package(x).processes_{userData.procID})
       funParams.OutputDirectory  = [userData_main.package(x).outputDirectory_  filesep 'ratio_tiffs'];
   else
       funParams.OutputDirectory = userData_main.package(x).processes_{userData.procID}.funParams_.OutputDirectory;
   end

   % if new process, create a new process with funParas and add to
   % MovieData and package's process list
   if isempty(userData_main.package(x).processes_{userData.procID})

       process = userData.procConstr(userData_main.MD(x), userData_main.package(x).outputDirectory_, funParams);
       userData_main.MD(x).addProcess( process )
       userData_main.package(x).setProcess(userData.procID, process )
       
   % if process exist, replace the funParams with the new one
   else
       userData_main.package(x).processes_{userData.procID}.setPara(funParams)
   end
   
   % If current process is changed, then assume funParams are changed in
   % all movies
   if userData.crtProc.procChanged_ 
       
       userData_main.package(x).processes_{userData.procID}.setProcChanged(true);
   end
   
    % Do sanity check - only check changed parameters
    procEx = userData_main.package(x).sanityCheck(false,'all');

    % Draw some bugs on the wall 
    for i = 1: length(procEx)
       if ~isempty(procEx{i})
           % Record the icon and message to user data
           userData_main.statusM(x).IconType{userData.procID} = 'warn';
           userData_main.statusM(x).Msg{i} = procEx{i}(1).message;
       end
    end   
end

% Write error report
if ~isempty( errorIndex )
    msg = '';
    for i = errorIndex
        msg = strcat(msg, sprintf('Movie %d  ', i));
    end
    msg = strcat(msg, sprintf('\n\nThe above movie(s) do not have the input channel you have selected. Please set up step 11, the output step in above movie(s) manually.'));
    titlemsg = 'Fail to be set up Output step of the following movie(s):';
    userData_main.msgboxGUI = msgboxGUI('title',titlemsg,'text', msg);
end

% Save user data
set(userData.mainFig, 'UserData', userData_main)

end
% -------------------------------------------------------------------------


set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);



function edit_path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_path as text
%        str2double(get(hObject,'String')) returns contents of edit_path as a double


% --- Executes during object creation, after setting all properties.
function edit_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_output.
function pushbutton_output_Callback(hObject, eventdata, handles)

pathname = uigetdir(pwd);
if isnumeric(pathname)
    return;
end

set(handles.edit_path, 'String', pathname);



function edit_factor_Callback(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);

% --- Executes during object creation, after setting all properties.
function edit_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_input1.
function listbox_input1_Callback(hObject, eventdata, handles)

contents1 = get(hObject, 'String');
chanIndex = get(hObject, 'Userdata');

id = get(hObject, 'Value');

if isempty(contents1) || isempty(id)
   return;
else
    set(handles.edit_dir, 'string', contents1{id}, 'Userdata',chanIndex(id));
end


% --- Executes during object creation, after setting all properties.
function listbox_input1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_input1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dir_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dir as text
%        str2double(get(hObject,'String')) returns contents of edit_dir as a double


% --- Executes during object creation, after setting all properties.
function edit_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end


% --- Executes on button press in checkbox_applytoall.
function checkbox_applytoall_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_applytoall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_applytoall
