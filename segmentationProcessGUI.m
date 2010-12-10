function varargout = segmentationProcessGUI(varargin)
% segmentationProcessGUI M-file for segmentationProcessGUI.fig
%      segmentationProcessGUI, by itself, creates a new segmentationProcessGUI or raises the existing
%      singleton*.
%
%      H = segmentationProcessGUI returns the handle to a new segmentationProcessGUI or the handle to
%      the existing singleton*.
%
%      segmentationProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in segmentationProcessGUI.M with the given input arguments.
%
%      segmentationProcessGUI('Property','Value',...) creates a new segmentationProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before segmentationProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to segmentationProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help segmentationProcessGUI

% Last Modified by GUIDE v2.5 26-Aug-2010 16:33:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @segmentationProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @segmentationProcessGUI_OutputFcn, ...
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


% --- Executes just before segmentationProcessGUI is made visible.
function segmentationProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
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

set(handles.text_copyright, 'String', userfcn_copyright)

userData = get(handles.figure1, 'UserData');
% Choose default command line output for segmentationProcessGUI
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
                    {userData_main.MD(userData_main.id).channels(channelIndex).channelPath_}, ...
                    'Userdata',channelIndex);    
            end
        end
end

% ---------------------- Parameter Setup -------------------------

if isempty(funParams.ThresholdValue)
    if funParams.MaxJump
       set(handles.checkbox_max, 'Value', 1);
       set(handles.edit_jump, 'Enable', 'on', 'String',...
                num2str(funParams.MaxJump));
    end
else
    set(handles.text_body3, 'Enable', 'on')
    set(handles.edit_coef, 'Enable', 'on')
    set(handles.pushbutton_add, 'Enable', 'on')
    set(handles.pushbutton_up, 'Enable', 'on')
    set(handles.pushbutton_down, 'Enable', 'on')
    set(handles.pushbutton_coef_delete, 'Enable', 'on')
    set(handles.listbox_coef1, 'Enable', 'on')
    
    set(handles.checkbox_auto, 'Value', 0)
    set(handles.checkbox_max, 'Enable', 'off')
    
    
    threshold = cell(1, length(funParams.ThresholdValue));
    
    for i = 1:length(funParams.ThresholdValue)
       threshold{i} = funParams.ThresholdValue(i); 
    end
    
    set(handles.listbox_coef1, 'String', threshold)
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
set(Img, 'UserData', userData.crtProc.getHelp(true))

% ----------------------------------------------------------------

% Update user data and GUI data
set(userData.mainFig, 'UserData', userData_main);
set(hObject, 'UserData', userData);

uicontrol(handles.pushbutton_done);
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = segmentationProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% Delete figure
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

if get(handles.checkbox_auto, 'value')
    if get(handles.checkbox_max, 'Value')
        % If both checkbox are checked
        if isnan(str2double(get(handles.edit_jump, 'String'))) ...
                || str2double(get(handles.edit_jump, 'String')) < 0
            errordlg('Please provide a valid input for ''Maximum threshold jump''.','Setting Error','modal');
            return;
        end    
    end
else
    threshold = get(handles.listbox_coef1, 'String');
    if isempty(threshold)
       errordlg('Please provide at least one threshold value.','Setting Error','modal')
       return
    elseif length(threshold) ~= 1 && length(threshold) ~= length(get(handles.listbox_2, 'String'))
       errordlg('Please provide the same number of threshold values as the input channels.','Setting Error','modal')
       return
    else
        threshold = str2double(threshold);
        if any(isnan(threshold)) || any(threshold < 0)
            errordlg('Please provide valid threshold values. Threshold cannot be a negative number.','Setting Error','modal')
            return            
        end
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

%---------Check if channel indexs are changed---------

channelIndex = get (handles.listbox_2, 'Userdata');
funParams = userData.crtProc.funParams_;

if ~isempty( setdiff(channelIndex, funParams.ChannelIndex) ) ...
    || ~isempty( setdiff(funParams.ChannelIndex, channelIndex) )

    % If channel indexs are changed, set procChanged to true
    userData.crtProc.setProcChanged(true);
end
    
% -------- Set parameter --------

if userData.crtProc.procChanged_ 
    
    % Get parameter
    
    funParams.ChannelIndex = channelIndex;

    if get(handles.checkbox_auto, 'value')
        % if automatic thresholding
        funParams.ThresholdValue = [ ];
        if get(handles.checkbox_max, 'value')
            funParams.MaxJump = str2double(get(handles.edit_jump,'String'));
        else
            funParams.MaxJump = 0;
        end
    else
        % if fixed thresholding
        if length(threshold) == 1
            funParams.ThresholdValue = repmat(threshold, [1 length(channelIndex)]);
        else
            funParams.ThresholdValue = threshold;
        end
    end
    % Set parameters
    userData.crtProc.setPara(funParams);
    

end


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
       userfcn_drawIcon(userData.handles_main,'warn',i,procEx{i}(1).message, true) % user data is retrieved, updated and submitted
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
   
   % Customize funParams to other movies 
   % ChannelIndex - all channels
   % OutputDirectory - pacakge output directory
       
   l = length(userData_main.MD(x).channels_);
   temp = arrayfun(@(x)(x > l),channelIndex, 'UniformOutput', true );
   funParams.ChannelIndex = channelIndex(logical(~temp));
   
   if get(handles.checkbox_auto, 'value')
       
       funParams.ThresholdValue = [ ];
   else
        if length(threshold) == 1
            funParams.ThresholdValue = repmat(threshold, [1 length(funParams.ChannelIndex)]);
        else
            funParams.ThresholdValue = threshold(logical(~temp));
        end
   end
   
   funParams.OutputDirectory  = [userData_main.package(x).outputDirectory_  filesep 'masks'];
   
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

% Save user data
set(userData.mainFig, 'UserData', userData_main)

end
% -------------------------------------------------------------------------

% Save user data
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);



% --- Executes on selection change in listbox_1.
function listbox_1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_1


% --- Executes during object creation, after setting all properties.
function listbox_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_2.
function listbox_2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_2


% --- Executes during object creation, after setting all properties.
function listbox_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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

switch get(hObject,'Value')
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


% --- Executes on button press in checkbox_auto.
function checkbox_auto_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of checkbox_auto
switch get(hObject, 'Value')
    case 0
        set(handles.text_body3, 'Enable', 'on')
        set(handles.edit_coef, 'Enable', 'on')
        set(handles.pushbutton_add, 'Enable', 'on')
        set(handles.pushbutton_up, 'Enable', 'on')
        set(handles.pushbutton_down, 'Enable', 'on')
        set(handles.pushbutton_coef_delete, 'Enable', 'on')
        set(handles.listbox_coef1, 'Enable', 'on')
        
        
        set(handles.checkbox_max, 'Enable', 'off', 'Value', 0);
        set(handles.edit_jump, 'Enable', 'off');
       
    case 1
        set(handles.text_body3, 'Enable', 'off')
        set(handles.edit_coef, 'Enable', 'off')
        set(handles.pushbutton_add, 'Enable', 'off')
        set(handles.pushbutton_up, 'Enable', 'off')
        set(handles.pushbutton_down, 'Enable', 'off')
        set(handles.pushbutton_coef_delete, 'Enable', 'off')
        set(handles.listbox_coef1, 'Enable', 'off')        
        
        set(handles.checkbox_max, 'Enable', 'on');
                
        
end
userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);


function edit_jump_Callback(hObject, eventdata, handles)
% Report to process object that the parameters are changed
userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);




% --- Executes during object creation, after setting all properties.
function edit_jump_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_jump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% Notify the package GUI that the setting panel is closed
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'helpFig') && ishandle(userData.helpFig)
   delete(userData.helpFig) 
end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);

% --- Executes on button press in checkbox_max.
function checkbox_max_Callback(hObject, eventdata, handles)


switch get(hObject, 'value')
    case 0
        set(handles.edit_jump, 'Enable', 'off');
    case 1
        set(handles.edit_jump, 'Enable', 'on');
end

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);




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


% --- Executes on selection change in listbox_coef1.
function listbox_coef1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_coef1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_coef1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_coef1


% --- Executes during object creation, after setting all properties.
function listbox_coef1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_coef1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);


% --- Executes on button press in pushbutton_up.
function pushbutton_up_Callback(hObject, eventdata, handles)
% % call back of 'Up' button
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

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);


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

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);


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

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);



function edit_coef_Callback(hObject, eventdata, handles)
% hObject    handle to edit_coef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_coef as text
%        str2double(get(hObject,'String')) returns contents of edit_coef as a double


% --- Executes during object creation, after setting all properties.
function edit_coef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_coef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
