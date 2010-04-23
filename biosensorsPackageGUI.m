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

% Last Modified by GUIDE v2.5 23-Apr-2010 10:42:29

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
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to biosensorsPackageGUI (see VARARGIN)

% Choose default command line output for biosensorsPackageGUI
handles.output = hObject;
% For test reason, define a dependency matrix in here
handles.dependM = [0 0 0 0;
                   1 0 0 0;
                   0 1 0 0;
                   0 1 0 0];
% Initial set up
userfcn_enable(find (any(handles.dependM,2)), 'off',handles);
% Load icon images from dialogicons.mat
load lccbGuiIcons.mat
supermap(1,:) = get(hObject,'color');
set(hObject,'colormap',supermap);

axes(handles.axes_help);
Img = image(questIconData); 
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off','YDir','reverse');

for i = 1:4
    eval (['axes(handles.axes_help' num2str(i) ')']);
    Img = image(questIconData); 
    set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
        'visible','off','YDir','reverse');
end
 eval (['axes(handles.axes_icon' num2str(1) ')']);
 Img = image(passIconData);
 set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off','YDir','reverse');
 eval (['axes(handles.axes_icon' num2str(2) ')']);
 Img = image(warnIconData);
 set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off','YDir','reverse');
 eval (['axes(handles.axes_icon' num2str(3) ')']);
 Img = image(errorIconData);
 set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off','YDir','reverse');

% Update handles structure
guidata(hObject, handles);

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
% hObject    handle to pushbutton_set1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
% hObject    handle to pushbutton_set2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
% hObject    handle to pushbutton_set3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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

% Hint: get(hObject,'Value') returns toggle state of checkbox_all
switch get(hObject,'value')
    case 1
        userfcn_enable(1:length(handles.dependM(:,1)),'on',handles,true);
    case 0
        userfcn_enable(find(any(handles.dependM,2)),'off',handles);
        for i = 1:length(handles.dependM(:,1))
            eval( ['set(handles.checkbox_', ...
                                          num2str(i),', ''value'', 0)'] );
        end
end


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_status.
function pushbutton_status_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_4.
function checkbox_4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_4
userfcn_lampSwitch(4, get(hObject,'value'), handles);


% --- Executes on button press in pushbutton_set4.
function pushbutton_set4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_set4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_show4.
function pushbutton_show4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function userfcn_enable (index, onoff, handles, enable)
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
%       onoff - 'on' or 'off'
%       handles - handles of control panel
%       select - (Optional) true or false. It provides a option to select/unselect 
%       the checkboxs that have been enabled/disabled.
% 
if nargin < 4
    enable = false;
end

for i = 1: length(index)
    eval (['set(handles.checkbox_', num2str(index(i)),...
                                        ',''enable'',''',onoff,''')']);
    eval (['set(handles.pushbutton_set', num2str(index(i)),...
                                        ',''enable'',''',onoff,''')']);
    eval (['set(handles.pushbutton_show', num2str(index(i)),...
                                        ',''enable'',''',onoff,''')']);
end
if enable
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
M = handles.dependM;

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
            for i =1:length(subindex)
                % Turn off and uncheck the follower checkboxes
                userfcn_enable(subindex(i),'off',handles,true);
                
%                 childrenI = find( M(:,subindex(i)) );
                userfcn_lampSwitch(subindex(i),0,handles);
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
