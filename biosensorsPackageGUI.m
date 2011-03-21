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

% Last Modified by GUIDE v2.5 17-Nov-2010 14:23:31

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
%
% biosensorsPackageGUI(MD)   MD: MovieData object
%
% Useful tools
%
% User Data:
%
%       userData.MD - array of MovieData object
%       userData.package - array of package (same length with userData.MD)
%       userData.crtPackage - the package of current MD
%       userData.id - the id of current MD on board
%
%       userData.dependM - dependency matrix
%       userdata.statusM - GUI status matrix
%       userData.optProcID - optional process id
%
%       userData.passIconData - pass icon image data
%       userData.errorIconData - error icon image data
%       userData.warnIconData - warning icon image data
%       userData.questIconData - help icon image data
%       userData.colormap - color map
%
%       userData.setFig - array of handles of (multiple) setting figures (may not exist)
%       userData.resultFig - array of handles of (multiple) result figures (may not exist)
%       userData.setupMovieDataFig - handle of (single) setupMovieData figure (may not exist)
%       userData.overviewFig - handles of (single) overviewMovieDataGUI figure (may not exist)
%       userData.packageHelpFig - handle of (single) help figure (may not exist)
%       userData.iconHelpFig - handle of (single) help figures (may not exist)
%       userData.processHelpFig - handle of (multiple) help figures (may not exist) 
%       userData.msgboxGUI - handle of message box
%       
%
% NOTE:
%   
%   userData.statusM - 1 x m stucture array, m is the number of Movie Data 
%                      this user data is used to save the status of movies
%                      when GUI is switching between different movie(s)
%                   
%   	fields: IconType - the type of status icons, 'pass', 'warn', 'error'
%               Msg - the message displayed when clicking status icons
%               Checked - 1 x n logical array, n is the number of processes
%                         used to save value of check box of each process
%               Visited - logical true or false, if the movie has been
%                         loaded to GUI before 
%
% Load movie data and recycle processes
userfcn_iniPackage_commonCode;


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

% In case the GUI has been called without argument
if (isfield(handles,'startMovieSelectorGUI') && handles.startMovieSelectorGUI)
    menu_file_open_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in checkbox_1.
function checkbox_1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_1

% Switch lamps
userfcn_checkAllMovies(1, get(hObject,'value'), handles);
userfcn_lampSwitch(1, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set_1.
function pushbutton_set_1_Callback(hObject, eventdata, handles)
% The process setting panel this button triggers is defined by 'procID', 
% who is the index of corresponding process in current package's process list
userData = get(handles.figure1, 'UserData');
procID = 1;
userData.setFig(procID) = segmentationProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_show_1.
function pushbutton_show_1_Callback(hObject, eventdata, handles)

procID = 1;
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'resultFig') && ishandle(userData.resultFig)
    
    delete(userData.resultFig)
end
%     userData.resultFig = userData.crtPackage.processes_{procID}.showResult;
    userData.resultFig = userData.crtPackage.processes_{procID}.resultDisplay;


set(handles.figure1, 'UserData', userData);


% --- Executes on button press in checkbox_2.
function checkbox_2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_2
userfcn_checkAllMovies(2, get(hObject,'value'), handles);
userfcn_lampSwitch(2, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set_2.
function pushbutton_set_2_Callback(hObject, eventdata, handles)
% The process setting panel this button triggers is defined by 'procID', 
% who is the index of corresponding process in current package's process list
userData = get(handles.figure1, 'UserData');
procID = 2;
userData.setFig(procID) = backgroundMasksProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_show_2.
function pushbutton_show_2_Callback(hObject, eventdata, handles)
procID = 2;
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'resultFig') && ishandle(userData.resultFig)
    delete(userData.resultFig)
end
%     userData.resultFig = userData.crtPackage.processes_{procID}.showResult;
    userData.resultFig = userData.crtPackage.processes_{procID}.resultDisplay;


set(handles.figure1, 'UserData', userData);


% --- Executes on button press in checkbox_3.
function checkbox_3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_3
userfcn_checkAllMovies(3, get(hObject,'value'), handles);
userfcn_lampSwitch(3, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set_3.
function pushbutton_set_3_Callback(hObject, eventdata, handles)
% The process setting panel this button triggers is defined by 'procID', 
% who is the index of corresponding process in current package's process list
userData = get(handles.figure1, 'UserData');
procID = 3;
userData.setFig(procID) = maskRefinementProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_show_3.
function pushbutton_show_3_Callback(hObject, eventdata, handles)
procID = 3;
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'resultFig') && ishandle(userData.resultFig)
    delete(userData.resultFig)
end
%     userData.resultFig = userData.crtPackage.processes_{procID}.showResult;
    userData.resultFig = userData.crtPackage.processes_{procID}.resultDisplay;


set(handles.figure1, 'UserData', userData);


% --- Executes on button press in checkbox_all.
function checkbox_all_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% userData = get(handles.figure1, 'UserData');
% l = size(userData.dependM, 1);

% switch get(hObject,'value')
%     case 1
%         userfcn_enable(1:length(userData.dependM(:,1)),'on',handles,true);
        
        % track checked status
%         for x = 1: length(userData.MD)
%             userData.statusM(x).Checked = ones(1, l);
%         end
        
%     case 0
%         k = [];
%         for i = 1: size(userData.dependM, 1)
%             if ~isempty(userData.crtPackage.processes_{i}) && ...
%                 userData.crtPackage.processes_{i}.success_ 
%                 k = [k, i];
%             end
%             eval( ['set(handles.checkbox_',num2str(i),',''value'',0)']  );
%         end
%         tempDependM = userData.dependM;
%         tempDependM(:,k) = zeros(size(userData.dependM,1),length(k));
%         userfcn_enable(find(any(tempDependM,2)),'off',handles,true);
        
        % track checked status
%         for x = 1: length(userData.MD)
%             userData.statusM(x).Checked = zeros(1, l);
%         end        
        
% end

% set(handles.figure1, 'UserData', userData)

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
for i = 1: length(userData.MD)
    userData.MD(i).saveMovieData
end
delete(handles.figure1);


% --- Executes on button press in pushbutton_status.
function pushbutton_status_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

% if newMovieDataGUI exist
if isfield(userData, 'overviewFig') && ishandle(userData.overviewFig)
    delete(userData.overviewFig)
end

userData.overviewFig = newMovieDataGUI('mainFig',handles.figure1, 'overview', userData.MD(userData.id));
set(handles.figure1, 'UserData', userData);



% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)

userfcn_pushbutton_run_common



% --- Executes on button press in checkbox_4.
function checkbox_4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_4
userfcn_checkAllMovies(4, get(hObject,'value'), handles);
userfcn_lampSwitch(4, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set_4.
function pushbutton_set_4_Callback(hObject, eventdata, handles)
% The process setting panel this button triggers is defined by 'procID', 
% who is the index of corresponding process in current package's process list
userData = get(handles.figure1, 'UserData');
procID = 4;
userData.setFig(procID) = darkCurrentCorrectionProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_show_4.
function pushbutton_show_4_Callback(hObject, eventdata, handles)
procID = 4;
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'resultFig') && ishandle(userData.resultFig)
    delete(userData.resultFig)
end
%     userData.resultFig = userData.crtPackage.processes_{procID}.showResult;
    userData.resultFig = userData.crtPackage.processes_{procID}.resultDisplay;


set(handles.figure1, 'UserData', userData);




% function userfcn_refreshGUI(handles)
% 
% userData = get(handles.figure1, 'UserData');
% l = size(userData.dependM,1);
% k = zeros(1,l);
% 
% % reset GUI
% userfcn_resetGUI(handles)
% 
% % Set movie data path
% set(handles.text_path, 'String', ...
%     [userData.MD(userData.id).movieDataPath_ userData.MD(userData.id).movieDataFileName_ ])
% 
% for i = 1: l
%    % Draw icons
%    if ~isempty(userData.statusM(userData.id).IconType{i})
%         userfcn_drawIcon(handles, userData.statusM(userData.id).IconType{i}, i, userData.statusM(userData.id).Msg{i}, false);
%    end
%    
%    % Bold the Name of Existing Processes
%    if ~isempty(userData.crtPackage.processes_{i})
%        eval([ 'set(handles.checkbox_',num2str(i),', ''FontWeight'',''bold'')' ])
%    end   
%    
%    % If process sucess = 1 or is checked, release the process from GUI
%    % enable/disable control
%    if ~isempty(userData.crtPackage.processes_{i}) && ...
%       userData.crtPackage.processes_{i}.success_ 
%       
%        k(i) = 1;
%        eval([ 'set(handles.pushbutton_show',num2str(i),', ''enable'', ''on'');'])
%    end
%    
%    % If process is checked, check and enable the process and enable decendent
%    % processes
%    if userData.statusM(userData.id).Checked(i)
%        k(i) = 1;
%        eval([ 'set(handles.checkbox_',num2str(i),', ''Value'',1, ''Enable'', ''on'')' ])
%        userfcn_lampSwitch(i, 1, handles)
%    end   
%    
% end
% 
% tempDependM = userData.dependM;
% tempDependM(:, logical(k)) = zeros(l, nnz(k));
% 
% % Checkbox enable/disable set-up
% userfcn_enable(find(any(tempDependM, 2)), 'off', handles);





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


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1, 'UserData');
% setFlag = getappdata(hObject, 'setFlag');

% Delete setting figures (multiple)
if isfield(userData, 'setFig')
    for i = 1: length(userData.setFig)
        if userData.setFig(i)~=0 && ishandle(userData.setFig(i))
            delete(userData.setFig(i))
        end
    end
end

% Close result figure
if isfield(userData, 'resultFig') && ishandle(userData.resultFig)
   delete(userData.resultFig) 
end

% If open, delete MovieData Overview GUI figure (single)
if isfield(userData, 'setupMovieDataFig') && ishandle(userData.setupMovieDataFig)
   delete(userData.setupMovieDataFig) 
end

% If open, delete MovieData Overview GUI figure (single)
if isfield(userData, 'overviewFig') && ishandle(userData.overviewFig)
   delete(userData.overviewFig) 
end

% Delete pre-defined package help dialog (single)
if isfield(userData, 'packageHelpFig') && ishandle(userData.packageHelpFig)
   delete(userData.packageHelpFig) 
end

% Delete pre-defined icon help dialog (single)
if isfield(userData, 'iconHelpFig') && ishandle(userData.iconHelpFig)
   delete(userData.iconHelpFig) 
end

% Delete pre-defined process help dialogssetting figures (multiple)
if isfield(userData, 'processHelpFig')
    for i = 1: length(userData.processHelpFig)
        if userData.processHelpFig(i)~=0 && ishandle(userData.processHelpFig(i))
            delete(userData.processHelpFig(i))
        end
    end
end

% msgboxGUI used for error reports
if isfield(userData, 'msgboxGUI') && ishandle(userData.msgboxGUI)
   delete(userData.msgboxGUI) 
end


% --- Executes on button press in checkbox_5.
function checkbox_5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_5
userfcn_checkAllMovies(5, get(hObject,'value'), handles);
userfcn_lampSwitch(5, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set_5.
function pushbutton_set_5_Callback(hObject, eventdata, handles)
% The process setting panel this button triggers is defined by 'procID', 
% who is the index of corresponding process in current package's process list
userData = get(handles.figure1, 'UserData');
procID = 5;
userData.setFig(procID) = shadeCorrectionProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_show_5.
function pushbutton_show_5_Callback(hObject, eventdata, handles)
procID = 5;
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'resultFig') && ishandle(userData.resultFig)
    delete(userData.resultFig)
end
%     userData.resultFig = userData.crtPackage.processes_{procID}.showResult;
    userData.resultFig = userData.crtPackage.processes_{procID}.resultDisplay;


set(handles.figure1, 'UserData', userData);

% --- Executes on button press in checkbox_6.
function checkbox_6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_6
userfcn_checkAllMovies(6, get(hObject,'value'), handles);
userfcn_lampSwitch(6, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set_6.
function pushbutton_set_6_Callback(hObject, eventdata, handles)
% The process setting panel this button triggers is defined by 'procID', 
% who is the index of corresponding process in current package's process list
userData = get(handles.figure1, 'UserData');
procID = 6;
userData.setFig(procID) = backgroundSubtractionProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_show_6.
function pushbutton_show_6_Callback(hObject, eventdata, handles)
procID = 6;
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'resultFig') && ishandle(userData.resultFig)
    delete(userData.resultFig)
end
%     userData.resultFig = userData.crtPackage.processes_{procID}.showResult;
    userData.resultFig = userData.crtPackage.processes_{procID}.resultDisplay;


set(handles.figure1, 'UserData', userData);


% --- Executes on button press in checkbox_7.
function checkbox_7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_7
userfcn_checkAllMovies(7, get(hObject,'value'), handles);
userfcn_lampSwitch(7, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set_7.
function pushbutton_set_7_Callback(hObject, eventdata, handles)
% The process setting panel this button triggers is defined by 'procID', 
% who is the index of corresponding process in current package's process
% list
userData = get(handles.figure1, 'UserData');
procID = 7;
userData.setFig(procID) = transformationProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_show_7.
function pushbutton_show_7_Callback(hObject, eventdata, handles)
procID = 7;
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'resultFig') && ishandle(userData.resultFig)
    delete(userData.resultFig)
end
%     userData.resultFig = userData.crtPackage.processes_{procID}.showResult;
    userData.resultFig = userData.crtPackage.processes_{procID}.resultDisplay;


set(handles.figure1, 'UserData', userData);


% --- Executes on button press in checkbox_8.
function checkbox_8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_8
userfcn_checkAllMovies(8, get(hObject,'value'), handles);
userfcn_lampSwitch(8, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set_8.
function pushbutton_set_8_Callback(hObject, eventdata, handles)
% The process setting panel this button triggers is defined by 'procID', 
% who is the index of corresponding process in current package's process list
userData = get(handles.figure1, 'UserData');
procID = 8;
userData.setFig(procID) = bleedthroughCorrectionProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_show_8.
function pushbutton_show_8_Callback(hObject, eventdata, handles)
procID = 8;
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'resultFig') && ishandle(userData.resultFig)
    delete(userData.resultFig)
end
%     userData.resultFig = userData.crtPackage.processes_{procID}.showResult;
    userData.resultFig = userData.crtPackage.processes_{procID}.resultDisplay;


set(handles.figure1, 'UserData', userData);


% --- Executes on button press in checkbox_9.
function checkbox_9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_9
userfcn_checkAllMovies(9, get(hObject,'value'), handles);
userfcn_lampSwitch(9, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set_9.
function pushbutton_set_9_Callback(hObject, eventdata, handles)
% The process setting panel this button triggers is defined by 'procID', 
% who is the index of corresponding process in current package's process list
userData = get(handles.figure1, 'UserData');
procID = 9;
userData.setFig(procID) = ratioingProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_show_9.
function pushbutton_show_9_Callback(hObject, eventdata, handles)
procID = 9;
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'resultFig') && ishandle(userData.resultFig)
    delete(userData.resultFig)
end

% userData.resultFig = userData.crtPackage.processes_{procID}.showResult;
userData.resultFig = userData.crtPackage.processes_{procID}.resultDisplay;


set(handles.figure1, 'UserData', userData);


% --- Executes on button press in checkbox_10.
function checkbox_10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_10
userfcn_checkAllMovies(10, get(hObject,'value'), handles);
userfcn_lampSwitch(10, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set_10.
function pushbutton_set_10_Callback(hObject, eventdata, handles)
% The process setting panel this button triggers is defined by 'procID', 
% who is the index of corresponding process in current package's process list
userData = get(handles.figure1, 'UserData');
procID = 10;
userData.setFig(procID) = photobleachCorrectionProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_show_10.
function pushbutton_show_10_Callback(hObject, eventdata, handles)
procID = 10;
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'resultFig') && ishandle(userData.resultFig)
    delete(userData.resultFig)
end
%     userData.resultFig = userData.crtPackage.processes_{procID}.showResult;
    userData.resultFig = userData.crtPackage.processes_{procID}.resultDisplay;


set(handles.figure1, 'UserData', userData);


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
status = web('http://lccb.hms.harvard.edu/', '-browser');
if status
    switch status
        case 1
            msg = 'System default web browser is not found.';
        case 2
            msg = 'System default web browser is found but could not be launched.';
        otherwise
            msg = 'Fail to open browser for unknown reason.';
    end
    warndlg(msg,'Fail to open browser','modal');
end

% --------------------------------------------------------------------
function menu_file_open_Callback(hObject, eventdata, handles)
% Call back function of 'New' in menu bar
userData = get(handles.figure1,'Userdata');
if isfield(userData,'MD')
    for i = 1: length(userData.MD)
        userData.MD(i).saveMovieData
    end
end
movieSelectorGUI(handles.packageName);
delete(handles.figure1)


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
if isfield(userData, 'MD')
    MD = userData.MD;
else
    delete(handles.figure1);
    return;
end

user_response = questdlg('Do you want to save the current progress?', ...
    'Biosensors Package Control Panel');
switch lower(user_response)
    case 'yes'
        for i = 1: length(userData.MD)
            userData.MD(i).saveMovieData
        end
        delete(handles.figure1);
    case 'no'
        delete(handles.figure1);
    case 'cancel'
end


% --- Executes on button press in checkbox_11.
function checkbox_11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_11
userfcn_checkAllMovies(11, get(hObject,'value'), handles);
userfcn_lampSwitch(11, get(hObject,'value'), handles);

% --- Executes on button press in pushbutton_set_11.
function pushbutton_set_11_Callback(hObject, eventdata, handles)
% The process setting panel this button triggers is defined by 'procID', 
% who is the index of corresponding process in current package's process list
userData = get(handles.figure1, 'UserData');
procID = 11;
userData.setFig(procID) = outputRatioProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_show_11.
function pushbutton_show_11_Callback(hObject, eventdata, handles)

procID = 11;
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'resultFig') && ishandle(userData.resultFig)
    delete(userData.resultFig)
end
%     userData.resultFig = userData.crtPackage.processes_{procID}.showResult;
    userData.resultFig = userData.crtPackage.processes_{procID}.resultDisplay;


set(handles.figure1, 'UserData', userData);


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

for i = 1: length(userData.MD)
    userData.MD(i).saveMovieData
end

set(handles.text_body3, 'Visible', 'on')
pause(1)
set(handles.text_body3, 'Visible', 'off')


% --------------------------------------------------------------------
function menu_file_save_Callback(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');
userData.MD(userData.id).saveMovieData

set(handles.text_body3, 'Visible', 'on')
pause(1)
set(handles.text_body3, 'Visible', 'off')


% --------------------------------------------------------------------
function menu_about_update_Callback(hObject, eventdata, handles)
status = web('http://lccb.hms.harvard.edu/software.html', '-browser');
if status
    switch status
        case 1
            msg = 'System default web browser is not found.';
        case 2
            msg = 'System default web browser is found but could not be launched.';
        otherwise
            msg = 'Fail to open browser for unknown reason.';
    end
    warndlg(msg,'Fail to open browser','modal');
end


% --- Executes on button press in pushbutton_left.
function pushbutton_left_Callback(hObject, eventdata, handles)
% userData.id
% userData.crtPackage
%
userData = get(handles.figure1, 'UserData');
l = length(userData.MD);

userData.statusM(userData.id).Checked = userfcn_saveCheckbox(handles);

userData.id = userData.id - 1;

if userData.id < 1
   userData.id = l;
end

userData.crtPackage = userData.package(userData.id);
set(handles.figure1, 'UserData', userData)

% Set up movie explorer
set(handles.popupmenu_movie, 'Value', userData.id)

% Set up GUI
if userData.statusM(userData.id).Visited
   userfcn_updateGUI(handles, 'refresh') 
else
   userfcn_updateGUI(handles, 'initialize') 
end
    

% --- Executes on button press in pushbutton_right.
function pushbutton_right_Callback(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');
l = length(userData.MD);

userData.statusM(userData.id).Checked = userfcn_saveCheckbox(handles);

userData.id = userData.id + 1;

if userData.id > l
   userData.id = mod(userData.id, l);
end

userData.crtPackage = userData.package(userData.id);
set(handles.figure1, 'UserData', userData)

% Set up movie explorer
set(handles.popupmenu_movie, 'Value', userData.id)

% Set up GUI
if userData.statusM(userData.id).Visited
   userfcn_updateGUI(handles, 'refresh') 
else
   userfcn_updateGUI(handles, 'initialize') 
end


% --- Executes on selection change in popupmenu_movie.
function popupmenu_movie_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

if get(hObject, 'Value') == userData.id
   return 
end

l = length(userData.MD);
userData.statusM(userData.id).Checked = userfcn_saveCheckbox(handles);

userData.id = get(hObject, 'Value');
userData.crtPackage = userData.package(userData.id);
set(handles.figure1, 'UserData', userData)

% Set up GUI
if userData.statusM(userData.id).Visited
   userfcn_updateGUI(handles, 'refresh') 
else
   userfcn_updateGUI(handles, 'initialize') 
end

% --- Executes during object creation, after setting all properties.
function popupmenu_movie_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
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
if strcmp(eventdata.Key, 'leftarrow')
    pushbutton_left_Callback(handles.pushbutton_left, [], handles);
end
if strcmp(eventdata.Key, 'rightarrow')
    pushbutton_right_Callback(handles.pushbutton_right, [], handles);
end



function checkbox_runall_Callback(hObject, eventdata, handles)



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
