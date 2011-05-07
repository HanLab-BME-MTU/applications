function varargout = bleedthroughCorrectionProcessGUI(varargin)
% BLEEDTHROUGHCORRECTIONPROCESSGUI M-file for bleedthroughCorrectionProcessGUI.fig
%      BLEEDTHROUGHCORRECTIONPROCESSGUI, by itself, creates a new BLEEDTHROUGHCORRECTIONPROCESSGUI or raises the existing
%      singleton*.
%
%      H = BLEEDTHROUGHCORRECTIONPROCESSGUI returns the handle to a new BLEEDTHROUGHCORRECTIONPROCESSGUI or the handle to
%      the existing singleton*.
%
%      BLEEDTHROUGHCORRECTIONPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BLEEDTHROUGHCORRECTIONPROCESSGUI.M with the given input arguments.
%
%      BLEEDTHROUGHCORRECTIONPROCESSGUI('Property','Value',...) creates a new BLEEDTHROUGHCORRECTIONPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bleedthroughCorrectionProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bleedthroughCorrectionProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bleedthroughCorrectionProcessGUI

% Last Modified by GUIDE v2.5 24-Aug-2010 11:16:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bleedthroughCorrectionProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @bleedthroughCorrectionProcessGUI_OutputFcn, ...
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


% --- Executes just before bleedthroughCorrectionProcessGUI is made visible.
function bleedthroughCorrectionProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
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
% Choose default command line output for bleedthroughCorrectionProcessGUI
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

% ---------------------- Channel Setup -------------------------

funParams = userData.crtProc.funParams_;

% Set up available input channels
set(handles.listbox_input1, 'String', {userData_main.MD(userData_main.id).channels_.channelPath_},...
        'Userdata', 1: length(userData_main.MD(userData_main.id).channels_));
    
% Set up input channel (one channel)
if ~isempty(funParams.ChannelIndex)
    set(handles.edit_dir, 'String', ...
        {userData_main.MD(userData_main.id).channels_(funParams.ChannelIndex).channelPath_} , ...
        'Userdata', funParams.ChannelIndex )
end

if ~isempty(funParams.ChannelIndex)
    set(handles.listbox_input1, 'Value', funParams.ChannelIndex(1))
end
    
% Set up bleed channels
set(handles.listbox_mask1, 'String', {userData_main.MD(userData_main.id).channels_.channelPath_},...
        'Userdata', 1: length(userData_main.MD(userData_main.id).channels_));    

% Set up bleedthrough channels
if ~isempty(funParams.BleedChannelIndex)
    set(handles.listbox_mask2, 'String', ...
        {userData_main.MD(userData_main.id).channels_(funParams.BleedChannelIndex).channelPath_}, ...
        'Userdata', funParams.BleedChannelIndex)
end

% ---------------------- Parameter Setup -------------------------

strBleedCoef = cell(1,length(funParams.BleedCoefficients));

for i = 1:length(funParams.BleedCoefficients)
    strBleedCoef{i} = funParams.BleedCoefficients(i);
end

set(handles.listbox_coef1, 'String', strBleedCoef)

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

uicontrol(handles.pushbutton_done)
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = bleedthroughCorrectionProcessGUI_OutputFcn(hObject, eventdata, handles) 
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
channelIndex = get (handles.edit_dir, 'Userdata');

if isempty(channelIndex)
   errordlg('Please select a channel as input channel.','Setting Error','modal') 
    return;
end

bleedChannelIndex = get(handles.listbox_mask2, 'Userdata');
if isempty(bleedChannelIndex)
   errordlg('Please select at least one bleedthrough channel from ''Available Bleed Channels''.','Setting Error','modal') 
    return;
else
    for i = 1: length(bleedChannelIndex)
        if bleedChannelIndex(i) == channelIndex
            errordlg('Bleedthrough channel cannot be input channel.','Setting Error','modal') 
            return;
        end
    end
end

if isempty(get(handles.listbox_coef1, 'String'))
   errordlg('Please add bleedthrough coefficients.','Setting Error','modal') 
   return;    
elseif length(bleedChannelIndex) ~= length(get(handles.listbox_coef1, 'String'))
        
   errordlg('The coefficients must be the same number as selected bleed channels.','Setting Error','modal') 
   return;
else
    bleedCoefficients = str2double(get(handles.listbox_coef1, 'String'));
    if any( isnan(bleedCoefficients) ) || any(bleedCoefficients < 0)
        errordlg('Please provide valid bleedthrough coefficients. Coefficients cannot be a negative number.','Setting Error','modal');
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

funParams = userData.crtProc.funParams_;

funParams.ChannelIndex = channelIndex;
funParams.BleedChannelIndex = bleedChannelIndex;
funParams.BleedCoefficients = bleedCoefficients;

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
    userData.crtPackage.setDepMatrix(9, userData.procID, 1);
    
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
   
   temp = arrayfun(@(x)(x > l),bleedChannelIndex, 'UniformOutput', true );
   funParams.BleedChannelIndex = bleedChannelIndex(logical(~temp));
   funParams.BleedCoefficients = bleedCoefficients(logical(~temp));
   funParams.OutputDirectory  = [userData_main.package(x).outputDirectory_  filesep 'bleedthrough_corrected_images'];
   
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
    msg = strcat(msg, sprintf('\n\nThe above movie(s) do not have the input channel you have selected. Please set up step 8, the Bleedthrough Correction, in above movie(s) manually.'));
    titlemsg = 'Fail to be set up Bleedthrough Correction step of the following movie(s):';
    userData_main.msgboxGUI = msgboxGUI('title',titlemsg,'text', msg);
end

% Save user data
set(userData.mainFig, 'UserData', userData_main)

end
% -------------------------------------------------------------------------


set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);

% --- Executes on button press in checkbox_all.
function checkbox_all_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of checkbox_all
contents1 = get(handles.listbox_mask1, 'String');

chanIndex1 = get(handles.listbox_mask1, 'Userdata');
chanIndex2 = get(handles.listbox_mask2, 'Userdata');

% Return if listbox1 is empty
if isempty(contents1)
    return;
end

switch get(hObject,'Value')
    case 1
        set(handles.listbox_mask2, 'String', contents1);
        chanIndex2 = chanIndex1;
    case 0
        set(handles.listbox_mask2, 'String', {}, 'Value',1);
        chanIndex2 = [ ];
end
set(handles.listbox_mask2, 'UserData', chanIndex2);


% --- Executes on button press in pushbutton_mask_select.
function pushbutton_mask_select_Callback(hObject, eventdata, handles)
% call back function of 'select' button
set(handles.listbox_mask2, 'Value', 1)

contents1 = get(handles.listbox_mask1, 'String');
contents2 = get(handles.listbox_mask2, 'String');
id = get(handles.listbox_mask1, 'Value');

% If channel has already been added, return;
chanIndex1 = get(handles.listbox_mask1, 'Userdata');
chanIndex2 = get(handles.listbox_mask2, 'Userdata');

for i = id
    if any(strcmp(contents1{i}, contents2) )
        continue;
    else
        contents2{end+1} = contents1{i};
        
        chanIndex2 = cat(2, chanIndex2, chanIndex1(i));

    end
end

set(handles.listbox_mask2, 'String', contents2, 'Userdata', chanIndex2);


% --- Executes on button press in pushbutton_mask_delete.
function pushbutton_mask_delete_Callback(hObject, eventdata, handles)
% Call back function of 'delete' button
contents = get(handles.listbox_mask2,'String');
id = get(handles.listbox_mask2,'Value');

% Return if list is empty
if isempty(contents) || isempty(id)
    return;
end

% Delete selected item
contents(id) = [ ];

% Delete userdata
chanIndex2 = get(handles.listbox_mask2, 'Userdata');
chanIndex2(id) = [ ];
set(handles.listbox_mask2, 'Userdata', chanIndex2);

% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if (id >length(contents) && id>1)
    set(handles.listbox_mask2,'Value',length(contents));
end
% Refresh listbox
set(handles.listbox_mask2,'String',contents);


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


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'helpFig') && ishandle(userData.helpFig)
   delete(userData.helpFig) 
end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_add.
function pushbutton_add_Callback(hObject, eventdata, handles)


set(handles.listbox_coef1, 'Value', 1);
text = get(handles.edit_coef, 'String');
if isempty(text)
    return;
end

if isnan(str2double(text)) || str2double(text) < 0 
    errordlg('Please provide a valid coefficient. Coefficient must be positive.','Setting Error','modal');
    return;
end

contents = get(handles.listbox_coef1, 'String');
contents{end + 1} = text;
set(handles.listbox_coef1, 'String', contents)
set(handles.edit_coef, 'String', '')


% --- Executes on button press in pushbutton_coef_delete.
function pushbutton_coef_delete_Callback(hObject, eventdata, handles)
% Call back function of 'delete' button
contents = get(handles.listbox_coef1,'String');
% Return if list is empty
if isempty(contents)
    return;
end
id = get(handles.listbox_coef1,'Value');

% Delete selected item
contents(id) = [ ];

% Refresh listbox
set(handles.listbox_coef1,'String',contents);
% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if (id>length(contents) && id>1)
    set(handles.listbox_coef1,'Value',length(contents));
end


% --- Executes on button press in pushbutton_up.
function pushbutton_up_Callback(hObject, eventdata, handles)
% call back of 'Up' button
id = get(handles.listbox_coef1,'Value');
contents = get(handles.listbox_coef1,'String');

% Return if list is empty
if isempty(contents) || isempty(id) || id == 1
    return;
end

temp = contents{id};
contents{id} = contents{id-1};
contents{id-1} = temp;

set(handles.listbox_coef1, 'string', contents);
set(handles.listbox_coef1, 'value', id-1);


% --- Executes on button press in pushbutton_down.
function pushbutton_down_Callback(hObject, eventdata, handles)
% Call back of 'Down' button
id = get(handles.listbox_coef1,'Value');
contents = get(handles.listbox_coef1,'String');

% Return if list is empty
if isempty(contents) || isempty(id) || id == length(contents)
    return;
end

temp = contents{id};
contents{id} = contents{id+1};
contents{id+1} = temp;

set(handles.listbox_coef1, 'string', contents);
set(handles.listbox_coef1, 'value', id+1);



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
