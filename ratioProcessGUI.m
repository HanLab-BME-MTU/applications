function varargout = ratioProcessGUI(varargin)
% RATIOPROCESSGUI M-file for ratioProcessGUI.fig
%      RATIOPROCESSGUI, by itself, creates a new RATIOPROCESSGUI or raises the existing
%      singleton*.
%
%      H = RATIOPROCESSGUI returns the handle to a new RATIOPROCESSGUI or the handle to
%      the existing singleton*.
%
%      RATIOPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RATIOPROCESSGUI.M with the given input arguments.
%
%      RATIOPROCESSGUI('Property','Value',...) creates a new RATIOPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ratioProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ratioProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ratioProcessGUI

% Last Modified by GUIDE v2.5 06-Apr-2011 17:07:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ratioProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ratioProcessGUI_OutputFcn, ...
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


% --- Executes just before ratioProcessGUI is made visible.
function ratioProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
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
% Choose default command line output for ratioProcessGUI
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

% Get icon infomation
userData.questIconData = userData_main.questIconData;
userData.colormap = userData_main.colormap;

% -----------------------Channel Setup--------------------------

funParams = userData.crtProc.funParams_;

set(handles.listbox_input, 'String', {userData_main.MD(userData_main.id).channels_.channelPath_},...
        'Userdata', 1: length(userData_main.MD(userData_main.id).channels_));
set(handles.listbox_mask, 'String', {userData_main.MD(userData_main.id).channels_.channelPath_},...
        'Userdata', 1: length(userData_main.MD(userData_main.id).channels_));
    
if ~isempty(funParams.ChannelIndex)
    
    set(handles.edit_nu_input, 'String', ...
        {userData_main.MD(userData_main.id).channels_(funParams.ChannelIndex(1)).channelPath_}, ...
        'Userdata',funParams.ChannelIndex(1));    
    
    set(handles.edit_de_input, 'String', ...
        {userData_main.MD(userData_main.id).channels_(funParams.ChannelIndex(2)).channelPath_}, ...
        'Userdata',funParams.ChannelIndex(2));      
end

if ~isempty(funParams.MaskChannelIndex)
    
    set(handles.edit_nu_mask, 'String', ...
        {userData_main.MD(userData_main.id).channels_(funParams.MaskChannelIndex(1)).channelPath_}, ...
        'Userdata',funParams.MaskChannelIndex(1));    
    
    set(handles.edit_de_mask, 'String', ...
        {userData_main.MD(userData_main.id).channels_(funParams.MaskChannelIndex(2)).channelPath_}, ...
        'Userdata',funParams.MaskChannelIndex(2));      
end
    
% ---------------------- Parameter Setup -------------------------

if ~funParams.ApplyMasks
        set(handles.listbox_mask, 'enable', 'off')
        set(handles.pushbutton_nu_mask, 'enable', 'off')
        set(handles.pushbutton_de_mask, 'enable', 'off')
        set(handles.checkbox_mask, 'Value', 0)
        set(handles.checkbox_newmask, 'enable', 'off', 'Value', 0)
        set(handles.text_input6, 'enable', 'off')
        set(handles.edit_de_mask, 'enable', 'off')
        set(handles.edit_nu_mask, 'enable', 'off')
else
    if ~funParams.CreateMasks
        set(handles.checkbox_newmask, 'Value', 0);
    end
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

% Update user data and GUI data
set(userData.mainFig, 'UserData', userData_main);
set(handles.figure1, 'UserData', userData);

uicontrol(handles.pushbutton_done)
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = ratioProcessGUI_OutputFcn(hObject, eventdata, handles) 
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
if isempty(get(handles.edit_nu_input, 'String'))
   errordlg('Please select a channel as numerator ''Input channel''.','Setting Error','modal') 
    return;
else
    if strcmp(get(handles.edit_nu_input, 'String'), ...
                            get(handles.edit_de_input, 'String'))
        errordlg('Numerator and denominator cannot be the same ''Input Channel''.','Setting Error','modal') 
        return;
    end
end

if isempty(get(handles.edit_de_input, 'String'))
   errordlg('Please select a channel as denominator ''Input Channel''.','Setting Error','modal') 
    return;
end

if get(handles.checkbox_mask, 'Value')
    if isempty(get(handles.edit_nu_mask, 'String'))
        errordlg('Please select a channel as numerator ''Mask Channel''.','Setting Error','modal') 
        return;
    else
        if strcmp(get(handles.edit_nu_mask, 'String'), ...
                            get(handles.edit_de_mask, 'String'))
            errordlg('Numerator and denominator cannot be the same ''Mask Channel''.','Setting Error','modal') 
            return;
        end
    end

    if isempty(get(handles.edit_de_mask, 'String'))
        errordlg('Please select a channel as denominator ''Mask Channel''.','Setting Error','modal') 
        return;
    end
    
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

% -------- Set parameter --------
channelIndex = [get(handles.edit_nu_input, 'Userdata'), get(handles.edit_de_input, 'Userdata')];
maskChannelIndex = [ get(handles.edit_nu_mask, 'Userdata'), get(handles.edit_de_mask, 'Userdata')];

funParams = userData.crtProc.funParams_;

funParams.ChannelIndex = channelIndex;
funParams.MaskChannelIndex = maskChannelIndex;

if get(handles.checkbox_mask, 'Value')
    
    funParams.ApplyMasks = true;
    
    if get(handles.checkbox_newmask, 'Value')
        funParams.CreateMasks = true;
    else
        funParams.CreateMasks = false;
    end
else
    funParams.ApplyMasks = false;
    funParams.CreateMasks = false;
    
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
   if any( channelIndex > l ) || any( maskChannelIndex > l )
       errorIndex = cat(2, errorIndex, x);
       continue
   else
       funParams.ChannelIndex = channelIndex;
   end
   
   funParams.OutputDirectory  = [userData_main.package(x).outputDirectory_  filesep 'ratio_images'];
   
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
   
    % Do sanity check - only check changed parameters
    procEx = userData_main.package(x).sanityCheck(false,'all');

    % Draw some bugs on the wall 
    for i = 1: length(procEx)
       if ~isempty(procEx{i})
           % Record the icon and message to user data
           userData_main.statusM(x).IconType{i} = 'warn';
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
    msg = strcat(msg, sprintf('\n\nThe above movie(s) do not have the channel you selected as numerator or denominator. Please set up step 9, the Ratioing in above movie(s) manually.'));
    titlemsg = 'Fail to be set up Ratioing step of the following movie(s):';
    userData_main.msgboxGUI = msgboxGUI('title',titlemsg,'text', msg);
end

% Save user data
set(userData.mainFig, 'UserData', userData_main)

end
% -------------------------------------------------------------------------


set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);



% --- Executes on button press in pushbutton_nu_input.
function pushbutton_nu_input_Callback(hObject, eventdata, handles)
% call back function of 'Choose as Numerator' button in input channel

contents1 = get(handles.listbox_input, 'String');
chanIndex = get(handles.listbox_input, 'Userdata');
id = get(handles.listbox_input, 'Value');

if isempty(contents1) || isempty(id)
   return;
% elseif strcmp( get(handles.edit_de_input, 'string'),  contents1{id})
%     return;
else
    set(handles.edit_nu_input, 'string', contents1{id}, 'Userdata', chanIndex(id));
end

% If numerator of mask channel is not set up. Set the same
% channel as the numerator of input channel
if get(handles.checkbox_mask, 'Value') && ...
        isempty(get(handles.edit_nu_mask, 'String')) 
   set(handles.edit_nu_mask, 'String', contents1{id}, 'Userdata', chanIndex(id)); 
end


% --- Executes on button press in checkbox_mask.
function checkbox_mask_Callback(hObject, eventdata, handles)
% Call back of 'Use Mask Channels' checkbox
switch get(hObject, 'Value')
    
    case 0

        set(handles.listbox_mask, 'enable', 'off')
        set(handles.pushbutton_nu_mask, 'enable', 'off')
        set(handles.pushbutton_de_mask, 'enable', 'off')
        set(handles.checkbox_newmask, 'enable', 'off', 'Value', 0)
        set(handles.text_input6, 'enable', 'off')
        set(handles.edit_nu_mask, 'Enable', 'off')
        set(handles.edit_de_mask, 'Enable', 'off')
        
    case 1

        set(handles.listbox_mask, 'enable', 'on')
        set(handles.pushbutton_nu_mask, 'enable', 'on')
        set(handles.pushbutton_de_mask, 'enable', 'on')
        set(handles.checkbox_newmask, 'enable', 'on')
        set(handles.text_input6, 'enable', 'on')
        set(handles.edit_nu_mask, 'Enable', 'inactive')
        set(handles.edit_de_mask, 'Enable', 'inactive')        
        
end


% --- Executes on button press in pushbutton_de_input.
function pushbutton_de_input_Callback(hObject, eventdata, handles)
% call back function of 'Choose as Denominator' button in input channel

contents1 = get(handles.listbox_input, 'String');
chanIndex = get(handles.listbox_input, 'Userdata');

id = get(handles.listbox_input, 'Value');

if isempty(contents1) || isempty(id)
   return;
% elseif strcmp( get(handles.edit_nu_input, 'string'),  contents1{id})
%     return;
else
    set(handles.edit_de_input, 'string', contents1{id},'Userdata',chanIndex(id));
end

% If denominator of mask channel is not set up. Set the same
% channel as the denominator of input channel
if get(handles.checkbox_mask, 'Value') && ...
        isempty(get(handles.edit_de_mask, 'String')) 
   set(handles.edit_de_mask, 'String', contents1{id}, 'Userdata', chanIndex(id)); 
end


% --- Executes on button press in pushbutton_nu_mask.
function pushbutton_nu_mask_Callback(hObject, eventdata, handles)
% call back function of 'Choose as Numerator' button in mask channel

contents1 = get(handles.listbox_mask, 'String');
chanIndex = get(handles.listbox_mask, 'userdata');
id = get(handles.listbox_mask, 'Value');

if isempty(contents1) || isempty(id)
   return;
% elseif strcmp( get(handles.edit_de_mask, 'string'),  contents1{id})
%     return;
else
    set(handles.edit_nu_mask, 'string', contents1{id}, 'Userdata', chanIndex(id));
end

% --- Executes on button press in pushbutton_de_mask.
function pushbutton_de_mask_Callback(hObject, eventdata, handles)
% call back function of 'Choose as Denominator' button in mask channel

contents1 = get(handles.listbox_mask, 'String');
chanIndex = get(handles.listbox_mask, 'Userdata');

id = get(handles.listbox_mask, 'Value');

if isempty(contents1) || isempty(id)
   return;
% elseif strcmp( get(handles.edit_nu_mask, 'string'),  contents1{id})
%     return;
else
    set(handles.edit_de_mask, 'string', contents1{id}, 'Userdata',chanIndex(id));
end


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
