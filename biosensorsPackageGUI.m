function varargout = biosensorsPackageGUI(varargin)
% biosensorsPackageGUI M-file for biosensorsPackageGUI.fig
%      biosensorsPackageGUI, by itself, creates a new biosensorsPackageGUI or raises the existing
%      singleton*.
%
%      H = biosensorsPackageGUI returns the handle to a new biosensorsPackageGUI or the handle to
%      the existing singleton*.
%
%      biosensorsPackageGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in biosensorsPackageGUI.M with the given input arguments.
%
%      biosensorsPackageGUI('Property','Value',...) creates a new biosensorsPackageGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before biosensorsPackageGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to biosensorsPackageGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help biosensorsPackageGUI

% Last Modified by GUIDE v2.5 04-May-2010 16:35:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @biosensorsPackageGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @biosensorsPackageGUI_OutputFcn, ...
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


% --- Executes just before biosensorsPackageGUI is made visible.
function biosensorsPackageGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% Useful tools
% GUI data:
%       handles.MD - the MovieData object
%       handles.dependM - dependency matrix
%       handles.passIconData - pass icon image data
%       handles.errorIconData - error icon image data
%       handles.warnIconData - warning icon image data
%       handles.setFig - array of handles of sub-windows
%
% App Data:
%       'setFlag' in figure1 - set flag of open sub window, open or close
%
userData = get(handles.figure1,'UserData');
% Choose default command line output for biosensorsPackageGUI
handles.output = hObject;
if nargin > 3
    % Pass the MovieData from setup MovieData GUI to Biosensors GUI
    userData.MD = varargin{1};
    % Handle of current package 
    userData.crtPackage = varargin{2};
    % Dependency matrix is defined in BioSensorsPackage class
    userData.dependM = userData.crtPackage.depMatrix_;
else
    % Default dependency matrix. For test reason, define a dependency 
    % matrix in here
    load movieData.mat
    userData.MD = MD;
    userData.crtPackage = MD.packages_{1};
    userData.dependM = [0 0 0 0
                       1 0 0 0
                       0 1 0 0
                       0 0 1 0];
end

% Load icon images from dialogicons.mat
load lccbGuiIcons.mat
% Save Icon data to GUI data
userData.passIconData = passIconData;
userData.errorIconData = errorIconData;
userData.warnIconData = warnIconData;
% Set figure colormap
supermap(1,:) = get(hObject,'color');
set(hObject,'colormap',supermap);
% Set up package help. Package icon is tagged as '0'
axes(handles.axes_help);
Img = image(questIconData); 
% handles.icon_help = Img;
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off','YDir','reverse');
set(Img,'ButtonDownFcn',@help_ButtonDownFcn);
set(Img,'tag','0');
% Set up process help. Process icons are tagged as '1','2' ... 'n'
for i = 1:size(userData.dependM, 1)
    eval (['axes(handles.axes_help' num2str(i) ')']);
    Img = image(questIconData);
%     eval(['handles.icon_help',num2str(i),' = Img; '])
    set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
        'visible','off','YDir','reverse');  
    set(Img,'ButtonDownFcn',@help_ButtonDownFcn);
    eval([ 'set(Img,''tag'',''',num2str(i),''');' ])
end

% set text body
set(handles.text_body1, 'string',[userData.crtPackage.name_ ' Package']);

% Set flag of sub window. Sub window open flag = 1, close flag = 0
setappdata(hObject, 'setFlag', zeros(1,size(userData.dependM,1)));

% Update handles structure
set(handles.figure1,'UserData',userData);
guidata(hObject, handles);

% all processes full sanity check
procEx = userData.crtPackage.sanityCheck(true, 'all');
k = [];
for i = 1: size(userData.dependM, 1)
   if ~isempty(procEx{i})
       l = 0;
       for j = 1: length(procEx{i})
           if strcmp(procEx{i}(j).identifier, 'lccb:set:fatal')
               if l <= 2
                   l = 2;
                   eInd = j;
               end
           elseif strcmp(procEx{i}(j).identifier, 'lccb:parachanged:warn')
               if l <= 1
                  l = 1;
                  eInd = j;
               end
           else
               if l == 0
                   eInd = j;
               end
           end
       end
       if l > 1
           userfcn_drawIcon(handles,'error',i,procEx{i}(eInd).message);
       else
           userfcn_drawIcon(handles,'warn',i,procEx{i}(eInd).message);
       end
   else
       if ~isempty(userData.crtPackage.processes_{i}) && ...
          userData.crtPackage.processes_{i}.success_ && ...
           ~userData.crtPackage.processes_{i}.procChanged_ && ...
           userData.crtPackage.processes_{i}.updated_
      
           userfcn_drawIcon(handles,'pass',i,'Current step is processed successfully') ;
           
       end
   end
   % If process's sucess = 1, release the process from GUI enable/disable
   % control
   if ~isempty(userData.crtPackage.processes_{i}) && ...
      userData.crtPackage.processes_{i}.success_ 
       k = [k, i];
       eval([ 'set(handles.pushbutton_show',num2str(i),', ''enable'', ''on'');']);
   end
end
tempDependM = userData.dependM;
tempDependM(:,k) = zeros(size(userData.dependM,1),length(k));
% Checkbox enable/disable set up
userfcn_enable(find (any(tempDependM,2)), 'off',handles);



% UIWAIT makes biosensorsPackageGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = biosensorsPackageGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox_1.
function checkbox_1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_1

% Switch lamps
userfcn_lampSwitch(1, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set1.
function pushbutton_set1_Callback(hObject, eventdata, handles)
% The process setting panel this button triggers is defined by 'procID', 
% who is the index of corresponding process in current package's process list
userData = get(handles.figure1, 'UserData');
procID = 1;
userData.setFig(procID) = thresholdProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_show1.
function pushbutton_show1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_2.
function checkbox_2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_2
userfcn_lampSwitch(2, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set2.
function pushbutton_set2_Callback(hObject, eventdata, handles)
% The process setting panel this button triggers is defined by 'procID', 
% who is the index of corresponding process in current package's process list
procID = 2;
userData = get(handles.figure1, 'UserData');
userData.setFig(procID) = ...
    backgroundMaskProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_show2.
function pushbutton_show2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_3.
function checkbox_3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_3
userfcn_lampSwitch(3, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set3.
function pushbutton_set3_Callback(hObject, eventdata, handles)
% The process setting panel this button triggers is defined by 'procID', 
% who is the index of corresponding process in current package's process list
procID = 3;
userData = get(handles.figure1, 'UserData');
userData.setFig(procID) = testProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_show3.
function pushbutton_show3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_all.
function checkbox_all_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1, 'UserData');
switch get(hObject,'value')
    case 1
        userfcn_enable(1:length(userData.dependM(:,1)),'on',handles,true);
    case 0
        k = [];
        for i = 1: size(userData.dependM, 1)
            if ~isempty(userData.crtPackage.processes_{i}) && ...
                userData.crtPackage.processes_{i}.success_ 
                k = [k, i];
            end
            eval( ['set(handles.checkbox_',num2str(i),',''value'',0)']  );
        end
        tempDependM = userData.dependM;
        tempDependM(:,k) = zeros(size(userData.dependM,1),length(k));
        userfcn_enable(find(any(tempDependM,2)),'off',handles,true);
        
end


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
MD = userData.MD;
save([userData.MD.movieDataPath_ userData.MD.movieDataFileName_], 'MD');
delete(handles.figure1);


% --- Executes on button press in pushbutton_status.
function pushbutton_status_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1, 'UserData'); 
disp('status');

% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)

userData = get(handles.figure1,'UserData');

procCheck = [ ]; % save id of checked processes 
procRun = [ ]; % save id of processes to run
nProcesses = size(userData.dependM,1);
MD = userData.MD;
for i = 1: nProcesses
    % collect the processes that are checked
    eval([ 'checked = get(handles.checkbox_',num2str(i),', ''value'');' ])
    if checked
        procCheck = horzcat(procCheck, i);
    end
end
if isempty(procCheck)
    errordlg('Please select a step to run','No Step Selected','modal');
    return
end
% Check if process exist
% Check if selected processes have alrady be successfully run
k = true;
for i = procCheck
    if isempty (userData.crtPackage.processes_{i})
        errordlg([num2str(i),' th step is not set up yet'], ...
            'Step Not Set Up','modal');
        return;
    end
    if ~( userData.crtPackage.processes_{i}.success_ ...
                && ~userData.crtPackage.processes_{i}.procChanged_ ) ...
            || ~userData.crtPackage.processes_{i}.updated_
        k = false;
        procRun = horzcat(procRun, i);
    end
end
if k
    warndlg('All selected steps have been processed successfully. Please change settings to run again', ...
        'Step Processed')
    return;
end


% Package full sanity check. Sanitycheck every checked process
procEx = userData.crtPackage.sanityCheck(true, procCheck);
k = {};
for i = procCheck
   if ~isempty(procEx{i})
       % Check if there is fatal error in exception array
       for j = 1: length(procEx{i})
           if strcmp(procEx{i}(j).identifier, 'lccb:set:fatal');
               k = horzcat(k,[num2str(i) ' ']);
%                userData.crtPackage.processes_{i}.setSuccess(false);
               userfcn_drawIcon(handles,'error',i,procEx{i}(j).message)
           end
       end
   end
end
% If setting error, stop and pop up notice
if ~isempty(k)
    errordlg(['Settings are incorrect in ' k{:} 'th step.' ...
        'Click corresponding error icon for further information'],...
             'Setting Error','modal');
    return
end

% Clear icons of selected processes
userfcn_drawIcon(handles,'clear',procRun,'');

% Set all running processes' sucess = false; 
for i = procRun
    userData.crtPackage.processes_{i}.setSuccess(false);
end

% Run the algorithms!
try
    for i = procRun
        userfcn_runProc_dfs(i,procRun,MD,handles);
    end
    
catch ME
    errordlg(ME.message,'Run Time Error','modal');
end


% --- Executes on button press in checkbox_4.
function checkbox_4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_4
userfcn_lampSwitch(4, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set4.
function pushbutton_set4_Callback(hObject, eventdata, handles)
% The process setting panel this button triggers is defined by 'procID', 
% who is the index of corresponding process in current package's process list
procID = 4;
userData = get(handles.figure1, 'userData');
userData.setFig(procID) = testProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_show4.
function pushbutton_show4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function userfcn_enable (index, onoff, handles, check)
% This is a user-defined function used to change the 'visible' property of
% uicontrols on control panel. The name of the uicontrols are pre-defined
% in the following way: 
%       checkbox: 
%               checkbox_1  checkbox_2 ...
%       pushbutton:
%               pushbutton_set1  pushbutton_set2 ...
%               pushbutton_show1  pushbutton_show2 ...
% Input: 
%       index - vector of check box index
%       onoff - enable or disable, 'on' or 'off'
%       handles - handles of control panel
%       check - (Optional) true or false. It provides a option to select/unselect 
%       the checkboxs that have been enabled/disabled.
% 
if nargin < 4
    check = false;
end

for i = 1: length(index)
    eval (['set(handles.checkbox_', num2str(index(i)),...
                                        ',''enable'',''',onoff,''')']);
    eval (['set(handles.pushbutton_set', num2str(index(i)),...
                                        ',''enable'',''',onoff,''')']);
%     eval (['set(handles.pushbutton_show', num2str(index(i)),...
%                                         ',''enable'',''',onoff,''')']);
end
if check
    switch onoff
        case 'on'
            for i = 1: length(index)
                eval( ['set(handles.checkbox_',...
                    num2str(index(i)),',''value'',1);'] );
            end
        case 'off'
            for i = 1: length(index)
                 eval( ['set(handles.checkbox_',...
                    num2str(index(i)),',''value'',0);'] );           
            end
    end

end

function userfcn_lampSwitch(index, value, handles)
% index - the index of current checkbox
% value - checked or unchecked
userData = get(handles.figure1, 'UserData');
M = userData.dependM;

if ~any(M(:,index))
   % if no follower exists, return.
        return;
else
    subindex = find(M(:,index));
    switch value
        % Checkbox is selected
        case 1
            for i = 1: length(subindex)
               parentI = find(M(subindex(i),:));
               for j = 1: length(parentI)
                   if ~eval(['get(handles.checkbox_',...
                                       num2str(parentI(j)),',''value'')'])
                       k = 1;
                       break;
                   else
                       k = 0;
                   end
               end
               if k == 1
                   continue;
               end
               % The following code will probably not be executed
               % Leave it here just in case design is changed
               % ------------------------------------------ %
               if eval(['get(handles.checkbox_', ...
                                      num2str(subindex(i)),',''value'')'])
                    userfcn_lampSwitch(subindex(i),1,handles)
               % ------------------------------------------ %
               else
                    % Turn on the subindex checkbox
                    userfcn_enable (subindex(i),'on',handles);
               end
            end
        % Checkbox is unselected
        case 0
            % If success = 1, release checkbox dependency enable/disable control
            if ~isempty(userData.crtPackage.processes_{index}) ...
                   && userData.crtPackage.processes_{index}.success_
                return;
            else
                for i =1:length(subindex)
                    % Turn off and uncheck the follower checkboxes
                    userfcn_enable(subindex(i),'off',handles,true);
                
                    userfcn_lampSwitch(subindex(i),0,handles);
                end
            end
        otherwise
            error(['User-defined error: unexpected value of ''value'' property',...
                            'in checkbox object']);
     end
            
end



function edit_notes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_notes as text
%        str2double(get(hObject,'String')) returns contents of edit_notes as a double


% --- Executes during object creation, after setting all properties.
function edit_notes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function help_ButtonDownFcn(hObject, eventdata)
% Call back function when help icon is clicked

tag = get(hObject,'tag');
handles = guidata(hObject);
userData = get(handles.figure1, 'UserData');
% get identity of help event creater
id = str2double(tag);
% For test purpose. No MovieData saved in handles
if ~isfield(userData,'MD')
    if id
        helpdlg(['This is help of process ',tag],'Help');
    else
        helpdlg(['This is a general help text for the current package'],...
                        'Help');
    end
    return;
end

% if id ==0 package help; id ~= 0 process help
if id
    % process help
    processName = userData.crtPackage.processClassNames_{id};
    % how do ppl handle similar situation in C or Java?
    eval(['helpdlg(',processName,'.getHelp)']);
%     helpdlg(text,'Help');
else
    % package help
    helpdlg(userData.crtPackage.getHelp, 'Help');
end

function icon_ButtonDownFcn(hObject, eventdata)

helpdlg( get(hObject,'UserData') ,'Warning');

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setFlag = getappdata(hObject, 'setFlag');
userData = get(handles.figure1, 'UserData');
if any(setFlag)
    index = find(setFlag);
    for i = index
        delete(userData.setFig(i));
    end
end
