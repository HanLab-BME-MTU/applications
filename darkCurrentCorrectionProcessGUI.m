function varargout = darkCurrentCorrectionProcessGUI(varargin)
% DARKCURRENTCORRECTIONPROCESSGUI M-file for darkCurrentCorrectionProcessGUI.fig
%      DARKCURRENTCORRECTIONPROCESSGUI, by itself, creates a new DARKCURRENTCORRECTIONPROCESSGUI or raises the existing
%      singleton*.
%
%      H = DARKCURRENTCORRECTIONPROCESSGUI returns the handle to a new DARKCURRENTCORRECTIONPROCESSGUI or the handle to
%      the existing singleton*.
%
%      DARKCURRENTCORRECTIONPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DARKCURRENTCORRECTIONPROCESSGUI.M with the given input arguments.
%
%      DARKCURRENTCORRECTIONPROCESSGUI('Property','Value',...) creates a new DARKCURRENTCORRECTIONPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before darkCurrentCorrectionProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to darkCurrentCorrectionProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help darkCurrentCorrectionProcessGUI

% Last Modified by GUIDE v2.5 07-May-2011 12:41:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @darkCurrentCorrectionProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @darkCurrentCorrectionProcessGUI_OutputFcn, ...
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


% --- Executes just before darkCurrentCorrectionProcessGUI is made visible.
function darkCurrentCorrectionProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% Available tools 
% UserData data:
%       userData.mainFig - handle of main figure
%       userData.handles_main - 'handles' of main figure
%       userData.procID - The ID of process in the current package
%       userData.crtProc - handle of current process
%       userData.crtPackage - handles of current package
%       userData.procConstr - constructor of current process
%       userData.userDir - save the default directory when user click 'Add Channel'
%
%       userData.questIconData - help icon image information
%       userData.colormap - color map information
%

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright)

userData = get(handles.figure1, 'UserData');
% Choose default command line output for darkCurrentCorrectionProcessGUI
handles.output = hObject;

% Get main figure handle and process id
t = find(strcmp(varargin,'mainFig'));
userData.mainFig = varargin{t+1};
userData.procID = varargin{t+2};
userData.handles_main = guidata(userData.mainFig);

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

% Get default user directory when clicking 'Add channel'
userData.userDir = pwd;

% Get icon infomation
userData.questIconData = userData_main.questIconData;
userData.colormap = userData_main.colormap;

% ---------------------- Channel Setup -------------------------

funParams = userData.crtProc.funParams_;

% Set up available input channels
set(handles.listbox_1, 'String', {userData_main.MD(userData_main.id).channels_.channelPath_},...
        'Userdata', 1: length(userData_main.MD(userData_main.id).channels_));
    
% Set up selected input data channels and channel index
parentI = find( userData.crtPackage.depMatrix_(userData.procID,:) );

if isempty(parentI) || ~isempty( userData.crtPackage.processes_{userData.procID} )
    
    % If process has no dependency, or process already exists, display saved channels 
    set(handles.listbox_2, 'String', ...
        {userData_main.MD(userData_main.id).channels_(funParams.ChannelIndex).channelPath_}, ...
        'Userdata',funParams.ChannelIndex);
    
elseif isempty( userData.crtPackage.processes_{userData.procID} )
    % If new process
        empty = false;
        for i = parentI
           if isempty(userData.crtPackage.processes_{i})
               empty = true;
               break;
           end
        end
            
        if ~empty

            % If all dependent processes exist
            channelIndex = userData.crtPackage.processes_{parentI(1)}.funParams_.ChannelIndex;
            for i = 2: length(parentI)
                channelIndex = intersect(channelIndex, ...
                    userData.crtPackage.processes_{parentI(i)}.funParams_.ChannelIndex);
            end  
            
            if ~isempty(channelIndex)
                set(handles.listbox_2, 'String', ...
                    {userData_main.MD(userData_main.id).channels_(channelIndex).channelPath_}, ...
                    'Userdata',channelIndex);    
            end
        end
end

% ---------------------- Parameter Setup -------------------------
if ~all(cellfun(@isempty,userData.crtProc.inFilePaths_(2,:)));
    set(handles.listbox_3, 'String', ...
        userData.crtProc.inFilePaths_(2,funParams.ChannelIndex ));
end

if ~funParams.MedianFilter
    set(handles.checkbox_medianfilter, 'Value',0)
end
if funParams.GaussFilterSigma >= 1
    set(handles.checkbox_gaussianfilter, 'Value', 1)
    set(handles.text_filters1, 'Enable','on');
    set(handles.edit_sigma, 'Enable','on',...
            'String',num2str(funParams.GaussFilterSigma))
end
    

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

set(handles.listbox_3, 'Value',1);


% Update user data and GUI data
set(userData.mainFig, 'UserData', userData_main);
set(handles.figure1, 'UserData', userData);

uicontrol(handles.pushbutton_done)
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = darkCurrentCorrectionProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
delete(handles.figure1);


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% Call back function of 'Apply' button
userData = get(handles.figure1, 'UserData');
userData_main = get(userData.mainFig, 'UserData');


% -------- Check user input --------

if isempty(get(handles.listbox_2, 'String'))
   errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal') 
    return;
end

if length(get(handles.listbox_2, 'String')) ~= ...
                            length(get(handles.listbox_3, 'String'))
   errordlg('Please provide the same number of dark-current image channels as input channels.','Setting Error','modal') 
    return;
end

if get(handles.checkbox_gaussianfilter, 'value')

        if isnan(str2double(get(handles.edit_sigma, 'String'))) ...
                || str2double(get(handles.edit_sigma, 'String')) < 0
            errordlg('Please provide a valid input for ''Sigma''.','Setting Error','modal');
            return;
            
        elseif str2double(get(handles.edit_sigma, 'String')) < 1
            warndlg('In order to perform Gaussian filtering, ''Sigma'' should be larger than 1.','Setting Error','modal');
            return;
        end    

end

% -------- Process Sanity check --------
% ( only check underlying data )

channelIndex = get (handles.listbox_2, 'Userdata');
funParams = userData.crtProc.funParams_;

% Get input index
tempChannelIndex = funParams.ChannelIndex;
funParams.ChannelIndex = channelIndex; 
userData.crtProc.setPara(funParams);

% Get dark current image path
temp = userData.crtProc.inFilePaths_(2,:);
userData.crtProc.setCorrectionImagePath(channelIndex, get(handles.listbox_3, 'String'));

try
    userData.crtProc.sanityCheck;
catch ME

    errordlg(ME.message,'Setting Error','modal');
    if ~isempty(temp)
        userData.crtProc.setCorrectionImagePath(1:length(temp), temp)
    end
    funParams.ChannelIndex = tempChannelIndex; 
    userData.crtProc.setPara(funParams);
    return;
end

funParams.ChannelIndex = tempChannelIndex; 
userData.crtProc.setPara(funParams);


% -------- Set parameter --------

  
funParams.ChannelIndex = channelIndex;
  
% Filter parameters
if get(handles.checkbox_medianfilter, 'Value')
    funParams.MedianFilter = true;
else
    funParams.MedianFilter = false;
end
if get(handles.checkbox_gaussianfilter, 'Value')
    funParams.GaussFilterSigma = str2double(get(handles.edit_sigma, 'String'));
else
    funParams.GaussFilterSigma = 0;
end

% Set parameters
userData.crtProc.setPara(funParams);

% --------------------------------------------------


% If this is a brand new process, attach current process to MovieData and 
% package's process list 
if isempty( userData.crtPackage.processes_{userData.procID} )
    
    % Add new process to both process lists of MovieData and current package
    userData_main.MD(userData_main.id).addProcess( userData.crtProc );
    userData.crtPackage.setProcess(userData.procID, userData.crtProc);
    
    % Customized to Biosensors Package (hard coded)
%     userData.crtPackage.setDepMatrix(5, userData.procID, 1);   
    
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

for x = 1: length(userData_main.MD)
    
   if x == userData_main.id
      continue 
   end
   
   l = length(userData_main.MD(x).channels_);
   
   temp = arrayfun(@(x)(x > l),channelIndex, 'UniformOutput', true );
   
   % Customize funParams to other movies 
   % ChannelIndex - all channels
   % OutputDirectory - pacakge output directory
   funParams.ChannelIndex = channelIndex(logical(~temp));
   funParams.OutputDirectory  = [userData_main.package(x).outputDirectory_  filesep 'dark_current_corrected_images'];
   
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
   
   % Additionally, for dark current process, set correction image path
   correctionPath = get(handles.listbox_3, 'String');
   userData_main.package(x).processes_{userData.procID}.setCorrectionImagePath(funParams.ChannelIndex, correctionPath(logical(~temp)));
   
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

% Save user data
set(userData.mainFig, 'UserData', userData_main)

end
% -------------------------------------------------------------------------

% Save user data
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);

% --- Executes on button press in pushbutton_deletechannel.
function pushbutton_deletechannel_Callback(hObject, eventdata, handles)
% Call back function of 'delete' button
contents = get(handles.listbox_3,'String');
% Return if list is empty
if isempty(contents)
    return;
end
num = get(handles.listbox_3,'Value');

% Delete selected item
contents(num) = [ ];

% Refresh listbox
set(handles.listbox_3,'String',contents);
% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if (num>length(contents) && num>1)
    set(handles.listbox_3,'Value',length(contents));
end

guidata(hObject, handles);

% --- Executes on button press in pushbutton_add.
function pushbutton_add_Callback(hObject, eventdata, handles)
% Call back of 'add' button
% Directories can be the same, no need to check directory pool
userData = get(handles.figure1, 'UserData');

set(handles.listbox_3, 'Value',1)

path = uigetdir(userData.userDir, 'Add Channels ...');
if path == 0
    return;
end
% Input validation function ...
% Get current list
contents = get(handles.listbox_3,'String');

% Add current formula to the listbox
contents{end+1} = path;
set(handles.listbox_3,'string',contents);

% Set user directory
sepDir = regexp(path, filesep, 'split');
dir = sepDir{1};
for i = 2: length(sepDir)-1
    dir = [dir filesep sepDir{i}];
end
userData.userDir = dir;

set(handles.figure1, 'Userdata', userData)
guidata(hObject, handles);

% --- Executes on button press in checkbox_all.
function checkbox_all_Callback(hObject, eventdata, handles)

% Hint: get(hObject,'Value') returns toggle state of checkbox_all
contents1 = get(handles.listbox_1, 'String');

chanIndex1 = get(handles.listbox_1, 'Userdata');
chanIndex2 = get(handles.listbox_2, 'Userdata');

% Return if listbox1 is empty
if isempty(contents1)
    return;
end

switch get(handles.checkbox_all,'Value')
    case 1
        set(handles.listbox_2, 'String', contents1);
        chanIndex2 = chanIndex1;
    case 0
        set(handles.listbox_2, 'String', {}, 'Value',1);
        chanIndex2 = [ ];
end
set(handles.listbox_2, 'UserData', chanIndex2);


% --- Executes on button press in pushbutton_select.
function pushbutton_select_Callback(hObject, eventdata, handles)
% call back function of 'select' button

contents1 = get(handles.listbox_1, 'String');
contents2 = get(handles.listbox_2, 'String');
id = get(handles.listbox_1, 'Value');

% If channel has already been added, return;
chanIndex1 = get(handles.listbox_1, 'Userdata');
chanIndex2 = get(handles.listbox_2, 'Userdata');

for i = id
    if any(strcmp(contents1{i}, contents2) )
        continue;
    else
        contents2{end+1} = contents1{i};
        
        chanIndex2 = cat(2, chanIndex2, chanIndex1(i));

    end
end

set(handles.listbox_2, 'String', contents2, 'Userdata', chanIndex2);



% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)
% Call back function of 'delete' button
contents = get(handles.listbox_2,'String');
id = get(handles.listbox_2,'Value');

% Return if list is empty
if isempty(contents) || isempty(id)
    return;
end

% Delete selected item
contents(id) = [ ];

% Delete userdata
chanIndex2 = get(handles.listbox_2, 'Userdata');
chanIndex2(id) = [ ];
set(handles.listbox_2, 'Userdata', chanIndex2);

% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if (id >length(contents) && id>1)
    set(handles.listbox_2,'Value',length(contents));
end
% Refresh listbox
set(handles.listbox_2,'String',contents);


% --- Executes on button press in checkbox_gaussianfilter.
function checkbox_gaussianfilter_Callback(hObject, eventdata, handles)

switch get(hObject, 'Value')
    case 0
        set(handles.text_filters1, 'Enable', 'off');
        set(handles.edit_sigma, 'Enable', 'off');
    case 1
        set(handles.text_filters1, 'Enable', 'on');
        set(handles.edit_sigma, 'Enable', 'on');        
end


% --- Executes on button press in pushbutton_up.
function pushbutton_up_Callback(hObject, eventdata, handles)

% call back of 'Up' button
id = get(handles.listbox_3,'Value');
contents = get(handles.listbox_3,'String');

% Return if list is empty
if isempty(contents) || isempty(id) || id == 1
    return;
end

temp = contents{id};
contents{id} = contents{id-1};
contents{id-1} = temp;

set(handles.listbox_3, 'string', contents);
set(handles.listbox_3, 'value', id-1);

% --- Executes on button press in pushbutton_down.
function pushbutton_down_Callback(hObject, eventdata, handles)

% Call back of 'Down' button
id = get(handles.listbox_3,'Value');
contents = get(handles.listbox_3,'String');

% Return if list is empty
if isempty(contents) || isempty(id) || id == length(contents)
    return;
end

temp = contents{id};
contents{id} = contents{id+1};
contents{id+1} = temp;

set(handles.listbox_3, 'string', contents);
set(handles.listbox_3, 'value', id+1);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'helpFig') && ishandle(userData.helpFig)
   delete(userData.helpFig) 
end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


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
