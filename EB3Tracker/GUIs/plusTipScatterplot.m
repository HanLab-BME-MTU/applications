function varargout = plusTipScatterplot(varargin)
% PLUSTIPSCATTERPLOT M-file for plusTipScatterplot.fig
%      PLUSTIPSCATTERPLOT, by itself, creates a new PLUSTIPSCATTERPLOT or raises the existing
%      singleton*.
%
%      H = PLUSTIPSCATTERPLOT returns the handle to a new PLUSTIPSCATTERPLOT or the handle to
%      the existing singleton*.
%
%      PLUSTIPSCATTERPLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLUSTIPSCATTERPLOT.M with the given input arguments.
%
%      PLUSTIPSCATTERPLOT('Property','Value',...) creates a new PLUSTIPSCATTERPLOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plusTipScatterplot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plusTipScatterplot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plusTipScatterplot

% Last Modified by GUIDE v2.5 20-May-2009 13:13:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plusTipScatterplot_OpeningFcn, ...
                   'gui_OutputFcn',  @plusTipScatterplot_OutputFcn, ...
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


% --- Executes just before plusTipScatterplot is made visible.
function plusTipScatterplot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plusTipScatterplot (see VARARGIN)

% Choose default command line output for plusTipScatterplot
handles.output = hObject;

handles.getStr=0;

handles.xaxis='growthSpeedMean';
handles.xlabel='Growth Speed (microns/min)';
handles.yaxis='pauseSpeedMean';
handles.ylabel='Pause Speed (microns/min)';

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes plusTipScatterplot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plusTipScatterplot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in selectGroupsPush.
function selectGroupsPush_Callback(hObject, eventdata, handles)
% hObject    handle to selectGroupsPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.movDataSet]=plusTipPickGroups;
guidata(hObject, handles);


% --- Executes on selection change in xaxisDrop.
function xaxisDrop_Callback(hObject, eventdata, handles)
% hObject    handle to xaxisDrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns xaxisDrop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xaxisDrop
val = get(hObject,'Value');
switch val
    case 1 % growth speed
        handles.xaxis='growthSpeedMean';
        handles.xlabel='Mean Growth Speed (microns/min)';
    case 2
        handles.xaxis='pauseSpeedMean';
        handles.xlabel='Mean Pause Speed (microns/min)';
    case 3
        handles.xaxis='shrinkSpeedMean';
        handles.xlabel='Mean Shrinkage Speed (microns/min)';
    case 4
        handles.xaxis='Ppause';
        handles.xlabel='Probability of pause';
    case 5
        handles.xaxis='Pcat';
        handles.xlabel='Probability of catastrophe';
end
guidata(hObject, handles);
        


% --- Executes during object creation, after setting all properties.
function xaxisDrop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xaxisDrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in makePlotButton.
function makePlotButton_Callback(hObject, eventdata, handles)
% hObject    handle to makePlotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%subset = ismember(movDataSet.groupName,unique(movDataSet.groupName)');
scattergroup = handles.movDataSet.groupName;
colorMap='rbgbcm'; %varycolor(length(unique(scattergroup)));
colorMap=colorMap(1:length(unique(scattergroup)));
figure
gscatter(handles.movDataSet.(handles.xaxis),handles.movDataSet.(handles.yaxis),scattergroup,colorMap)
xlabel(handles.xlabel)
ylabel(handles.ylabel)
legend('location','best')



% --- Executes on selection change in yaxisDrop.
function yaxisDrop_Callback(hObject, eventdata, handles)
% hObject    handle to yaxisDrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns yaxisDrop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from yaxisDrop
val = get(hObject,'Value');
switch val
    case 1 % growth speed
        handles.yaxis='growthSpeedMean';
        handles.ylabel='Mean Growth Speed (microns/min)';
    case 2
        handles.yaxis='pauseSpeedMean';
        handles.ylabel='Mean Pause Speed (microns/min)';
    case 3
        handles.yaxis='shrinkSpeedMean';
        handles.ylabel='Mean Shrinkage Speed (microns/min)';
    case 4
        handles.yaxis='Ppause';
        handles.ylabel='Probability of pause';
    case 5
        handles.yaxis='Pcat';
        handles.ylabel='Probability of catastrophe';
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function yaxisDrop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yaxisDrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in getQueryStr_Check.
function getQueryStr_Check_Callback(hObject, eventdata, handles)
% hObject    handle to getQueryStr_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of getQueryStr_Check
handles.getStr=get(hObject,'Value');
guidata(hObject, handles);

