function varargout = checkIdlistGui(varargin)
% CHECKIDLISTGUI M-file for checkIdlistGui.fig
%      CHECKIDLISTGUI, by itself, creates a new CHECKIDLISTGUI or raises the existing
%      singleton*.
%
%      H = CHECKIDLISTGUI returns the handle to a new CHECKIDLISTGUI or the handle to
%      the existing singleton*.
%
%      CHECKIDLISTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHECKIDLISTGUI.M with the given input arguments.
%
%      CHECKIDLISTGUI('Property','Value',...) creates a new CHECKIDLISTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before checkIdlistGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to checkIdlistGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help checkIdlistGui

% Last Modified by GUIDE v2.5 07-Jan-2008 11:29:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @checkIdlistGui_OpeningFcn, ...
                   'gui_OutputFcn',  @checkIdlistGui_OutputFcn, ...
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


% --- Executes just before checkIdlistGui is made visible.
function checkIdlistGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to checkIdlistGui (see VARARGIN)

% preset output
handles.outA2a = 0;
handles.outCell = -1;

% set title
if length(varargin)>1
    set(handles.ci_figure,'Name',varargin{2});
end

% store info-list
handles.pulldownList = varargin{1};

% set pulldown-text
pdHandles = findall(handles.ci_figure,'Style','popupmenu');
set(pdHandles,'String',handles.pulldownList(:,1));
handles.pdHandles = pdHandles(end:-1:1); % reverse list

% collect edit, text handles
etHandles = cell(4,2);
etHandles{1,1} = [handles.edit1,handles.edit2];
etHandles{1,2} = [handles.text1,handles.text2];
etHandles{2,1} = [handles.edit3,handles.edit4];
etHandles{2,2} = [handles.text3,handles.text4];
etHandles{3,1} = [handles.edit5,handles.edit6];
etHandles{3,2} = [handles.text5,handles.text6];
etHandles{4,1} = [handles.edit7,handles.edit8];
etHandles{4,2} = [handles.text7,handles.text8];

handles.etHandles = etHandles;

% hide options
set([handles.etHandles{:}],'Visible','off');


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes checkIdlistGui wait for user response (see UIRESUME)
 uiwait(handles.ci_figure);


% --- Outputs from this function are returned to the command line.
function varargout = checkIdlistGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.outCell;
varargout{2} = handles.outA2a;
delete(handles.ci_figure);


% --- Executes on button press in ci_okButton.
function ci_okButton_Callback(hObject, eventdata, handles)
% hObject    handle to ci_okButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% read current selections
selection = get(handles.pdHandles,'Value');
selection = [selection{:}];
% loop selection to write.
outCell = cell(4,3);
for i = 1:4
    % has something been selected?
    if selection(i) > 1
        
        % write selection into out cell
        outCell(i,1) = handles.pulldownList(selection(i),2);
        
        % check for input
        switch handles.pulldownList{selection(i),3}
            case 0
                % no additional input needed
            case {1}
                % one additional input needed, numeric
                opt = eval(sprintf('[%s];',get(handles.etHandles{i,1}(1),'String')));
                if isempty(opt) || any(isnan(opt))
                    h = errordlg('Please provide correct, nonempty input');
                    uiwait(h);
                    return
                else
                    outCell{i,2} = opt;
                end
                
            case 2
                % one additional input needed, string
                opt = get(handles.etHandles{i,1}(1),'String');
                if isempty(opt) 
                    h = errordlg('Please provide correct, nonempty input');
                    uiwait(h);
                    return
                else
                    outCell{i,2} = opt;
                end
                
            case 11
                % two additional inputs, numeric
                
                 opt = eval(sprintf('[%s];',get(handles.etHandles{i,1}(1),'String')));
                if isempty(opt) || any(isnan(opt))
                    h = errordlg('Please provide correct, nonempty input');
                    uiwait(h);
                    return
                else
                    outCell{i,2} = opt;
                end
                
                 opt = eval(sprintf('[%s];',get(handles.etHandles{i,1}(2),'String')));
                if isempty(opt) || any(isnan(opt))
                    h = errordlg('Please provide correct, nonempty input');
                    uiwait(h);
                    return
                else
                    outCell{i,3} = opt;
                end
                
            case 12
                % numeric, string
                    opt = eval(sprintf('[%s];',get(handles.etHandles{i,1}(1),'String')));
                if isempty(opt) || any(isnan(opt))
                    h = errordlg('Please provide correct, nonempty input');
                    uiwait(h);
                    return
                else
                    outCell{i,2} = opt;
                end
                
                 opt = get(handles.etHandles{i,1}(2),'String');
                if isempty(opt) 
                    h = errordlg('Please provide correct, nonempty input');
                    uiwait(h);
                    return
                else
                    outCell{i,3} = opt;
                end
                
            case 21
                % string, numeric
                 opt = get(handles.etHandles{i,1}(1),'String');
                if isempty(opt) 
                    h = errordlg('Please provide correct, nonempty input');
                    uiwait(h);
                    return
                else
                    outCell{i,2} = opt;
                end
                
                 opt = eval(sprintf('[%s];',get(handles.etHandles{i,1}(2),'String')));
                if isempty(opt) || any(isnan(opt))
                    h = errordlg('Please provide correct, nonempty input');
                    uiwait(h);
                    return
                else
                    outCell{i,3} = opt;
                end
                
            case 22
                % string, string
                    opt = get(handles.etHandles{i,1}(1),'String');
                if isempty(opt) 
                    h = errordlg('Please provide correct, nonempty input');
                    uiwait(h);
                    return
                else
                    outCell{i,2} = opt;
                end
                
                  opt = get(handles.etHandles{i,1}(2),'String');
                if isempty(opt) 
                    h = errordlg('Please provide correct, nonempty input');
                    uiwait(h);
                    return
                else
                    outCell{i,3} = opt;
                end
                
        end % switch
    end % if selected
end % loop

% remove non-selections
outCell(selection==1,:) = [];

% save
handles.outCell = outCell;

guidata(hObject,handles)
uiresume(handles.ci_figure)
    

% --- Executes on button press in ci_a2aButton.
function ci_a2aButton_Callback(hObject, eventdata, handles)
% hObject    handle to ci_a2aButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% remember apply to all
handles.outA2a = 1;
% simply call ok button
ci_okButton_Callback(hObject,eventdata,handles);



% --- Executes on button press in ci_cancelButton.
function ci_cancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to ci_cancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ask for abort
a = questdlg('Do you really want to abort?','Warning','Yes','No','No');
if strcmp(a,'Yes')
    % kill figure
    uiresume(handles.ci_figure);
end

% --- Executes on selection change in popupmenu.
function popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get number of popupmenu
popTag = get(hObject,'Tag');
popNumber = str2double(popTag(end));

% get selection
selection = get(hObject,'Value');

% look up selection
selectionCell = handles.pulldownList(selection,:);
etH = handles.etHandles(popNumber,:);

% check for visibility of options and insert defaults
switch selectionCell{3}
    case 0
        set([etH{:}],'Visible','off');
    case {1,2}
        % write default
        set(etH{1}(1),'String',num2str(selectionCell{4}),'Visible','on')
        % label
        set(etH{2}(1),'String',selectionCell{5},'Visible','on')
    case {11,12,21,22}
        % write defaults
        set(etH{1}(1),'String',num2str(selectionCell{4}),'Visible','on')
        set(etH{1}(2),'String',num2str(selectionCell{6}),'Visible','on')
        % labels
        set(etH{2}(1),'String',selectionCell{5},'Visible','on')
        set(etH{2}(2),'String',selectionCell{7},'Visible','on')
end



% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


