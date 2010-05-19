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

% Last Modified by GUIDE v2.5 19-May-2010 10:33:03

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
% User Data:
%       userData.MD - the MovieData object
%       userData.dependM - dependency matrix
%       userData.crtPackage - current package
%
%       userData.passIconData - pass icon image data
%       userData.errorIconData - error icon image data
%       userData.warnIconData - warning icon image data
%       userData.setFig - array of handles of sub-windows (may not exist)
%       userData.setupMovieDataFig - handle of setupMovieData figure (may not exist)
%       userData.overviewFig - handles of overviewMovieDataGUI figure (may not exist)
%       
% App Data:
%       'setFlag' in figure1 - set flag of open sub window, open or close
%       'setupMovieDataFlag' in figure1 - flag of setupMovieData figure
%       'overviewFlag' in figure1 - flag of view status figure
%

handles.output = hObject;
userData = get(handles.figure1,'UserData');

% Call package GUI error

if nargin == 3
    error('User-defined: Please call biosensors control panel with a MovieData object. E.g. biosensorsPackageGUI(movieDataObject)');
%     userData.setupMovieDataFig = setupMovieDataGUI;
    
%     set(handles.figure1,'UserData',userData);
end

% ----------------------------- Load MovieData ----------------------------

MD = varargin{1};

% I. Before loading MovieData, firstly check if the current package exists
    packageExist = false;

    for i = 1: length(MD.packages_)
        if strcmp(class(MD.packages_{i}), 'BioSensorsPackage')
            % TODO: dlg box: previous package exists; Save package and
            % proceeds  
            userData.crtPackage = MD.packages_{i};
            packageExist = true;
            break;
        end
    end
    
    if ~packageExist
        % No same package is found.. create a new package object
        MD.addPackage( BioSensorsPackage(MD, MD.outputDirectory_) )
        userData.crtPackage = MD.packages_{end};
    end


% II. Check if existing processes can be recycled. (No use in first release of software)

if ~packageExist && ~isempty(MD.processes_)
    for i = 1: length(userData.crtPackage.processClassNames_)
        % ith process in package's process list
        % 'sameProcID' save the ID of existing processes in MovieData that 
        % could be recycled by current package. 

        sameProcID = [];
        sameProcName = {};
        
        for j = 1: length(MD.processes_)
            % jth available process in MovieData's process pool
           if strcmp(class(MD.processes_{j}), ...
                        userData.crtPackage.processClassNames_{i}) 
                sameProcID = horzcat (sameProcID, j);
                sameProcName = horzcat(sameProcName, ...
                                            {MD.processes_{j}.name_});
           end
        end
        
        if ~isempty(sameProcID)
            [select, ok] = listdlg('ListString', sameProcName, ...
                 'SelectionMode','single','ListSize',[420 200], ...
                 'Name','Select an existing process','PromptString', ...
                 {['There are existing ',sameProcName{1},' processes saved during previous sessions. You can choose to:'],...
                 '',...
                 '1. Use one of the following existing processes by click ''Select''',...
                 '2. click ''New'' to start over this process'} , ...
                 'OKString','Select','CancelString','New'); 
             if ok 
                 % Add process to current package
                 userData.crtPackage.setProcess(i,MD.processes_{sameProcID(select)});
             end
        end
    end
end
save([MD.movieDataPath_ MD.movieDataFileName_], 'MD');

    userData.MD = MD;

    % Dependency matrix is defined in BioSensorsPackage class
    userData.dependM = userData.crtPackage.depMatrix_;

% -----------------------Load and set up icons----------------------------

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
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off','YDir','reverse');
set(Img,'ButtonDownFcn',@help_ButtonDownFcn);
set(Img,'tag','0');

% Set up process help. Process icons are tagged as '1','2' ... 'n'
for i = 1:size(userData.dependM, 1)
    eval (['axes(handles.axes_help' num2str(i) ')']);
    Img = image(questIconData);
    set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
        'visible','off','YDir','reverse');  
    set(Img,'ButtonDownFcn',@help_ButtonDownFcn);
    eval([ 'set(Img,''tag'',''',num2str(i),''');' ])
end

% --------------------------Other GUI settings-----------------------------

% set text body
set(handles.text_body1, 'string',[userData.crtPackage.name_ ' Package']);

% Set flag of sub window. Sub window open flag = 1, close flag = 0
setappdata(hObject, 'setFlag', zeros(1,size(userData.dependM,1)));
setappdata(hObject, 'setupMovieDataFlag', 0);
setappdata(hObject,'overviewFlag',0);

% Update handles structure
set(handles.figure1,'UserData',userData);
guidata(hObject, handles);

% --------------------------Package Sanity Check---------------------------

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

% --------------------------Set up uicontrols------------------------------

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
userData.setFig(procID) = segmentationProcessGUI('mainFig',handles.figure1,procID);
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
userData = get(handles.figure1, 'UserData');
procID = 2;
userData.setFig(procID) = backgroundMasksProcessGUI('mainFig',handles.figure1,procID);
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
save([MD.movieDataPath_ MD.movieDataFileName_], 'MD');
delete(handles.figure1);


% --- Executes on button press in pushbutton_status.
function pushbutton_status_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

userData.overviewFig = overviewMovieDataGUI('mainFig',handles.figure1);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


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
for i = procCheck
    if isempty (userData.crtPackage.processes_{i})
        errordlg([num2str(i),' th step is not set up yet'], ...
            'Step Not Set Up','modal');
        return;
    end    
end

% Check if selected processes have alrady be successfully run
% If force run, re-run every process that is checked
if ~get(handles.checkbox_forcerun, 'Value')

    k = true;
    for i = procCheck

        if ~( userData.crtPackage.processes_{i}.success_ ...
                && ~userData.crtPackage.processes_{i}.procChanged_ ) ...
                || ~userData.crtPackage.processes_{i}.updated_
            k = false;
            procRun = horzcat(procRun, i);
        end
    end
    if k
        warndlg('All selected steps have been processed successfully. Please check ''Force Run'' to re-run the steps.', ...
        'Step Processed')
        return;
    end
else
    procRun = procCheck;
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
userData = get(handles.figure1, 'UserData');
setFlag = getappdata(hObject, 'setFlag');

% Delete setting figures
if any(setFlag)
    index = find(setFlag);
    for i = index
        delete(userData.setFig(i));
    end
end

% If open, delete setupMovieData GUI figure
if getappdata(hObject,'setupMovieDataFlag');
    delete(userData.setupMovieDataFig);
end

% If open, delete setupMovieData GUI figure
if getappdata(hObject,'overviewFlag');
    delete(userData.overviewFig);
end


% --- Executes on button press in checkbox_5.
function checkbox_5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_5
userfcn_lampSwitch(5, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set5.
function pushbutton_set5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_set5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_show5.
function pushbutton_show5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_6.
function checkbox_6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_6
userfcn_lampSwitch(6, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set6.
function pushbutton_set6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_set6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_show6.
function pushbutton_show6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_7.
function checkbox_7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_7
userfcn_lampSwitch(7, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set7.
function pushbutton_set7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_set7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_show7.
function pushbutton_show7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_8.
function checkbox_8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_8
userfcn_lampSwitch(8, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set8.
function pushbutton_set8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_set8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_show8.
function pushbutton_show8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_9.
function checkbox_9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_9
userfcn_lampSwitch(9, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set9.
function pushbutton_set9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_set9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_show9.
function pushbutton_show9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_10.
function checkbox_10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_10
userfcn_lampSwitch(10, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set10.
function pushbutton_set10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_set10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_show10.
function pushbutton_show10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_help_lccb_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help_lccb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_about_lccb_Callback(hObject, eventdata, handles)
% hObject    handle to menu_about_lccb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_open_Callback(hObject, eventdata, handles)
% Call back function of 'New' in menu bar
userData = get(handles.figure1, 'UserData');

% Open setupMovieData, notify new MovieData object is needed.
userData.setupMovieDataFig = ...
    setupMovieDataGUI('mainFig',handles.figure1);
set(handles.figure1, 'UserData', userData);



% --------------------------------------------------------------------
function menu_file_exit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);


% --- Executes on button press in checkbox_forcerun.
function checkbox_forcerun_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_forcerun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_forcerun


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1,'Userdata');
MD = userData.MD;

user_response = questdlg(['Do you want to save the current progress to ',MD.movieDataFileName_,'?'], ...
    'BioSensors Package Control Panel');
switch lower(user_response)
    case 'yes'
        save([MD.movieDataPath_ MD.movieDataFileName_], 'MD');
        delete(handles.figure1);
    case 'no'
        delete(handles.figure1);
    case 'cancel'
end
