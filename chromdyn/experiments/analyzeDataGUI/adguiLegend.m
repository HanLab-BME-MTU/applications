function varargout = adguiLegend(varargin)
% ADGUILEGEND M-file for adguiLegend.fig
%      ADGUILEGEND, by itself, creates a new ADGUILEGEND or raises the existing
%      singleton*.
%
%      H = ADGUILEGEND returns the handle to a new ADGUILEGEND or the handle to
%      the existing singleton*.
%
%      ADGUILEGEND('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADGUILEGEND.M with the given input arguments.
%
%      ADGUILEGEND('Property','Value',...) creates a new ADGUILEGEND or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before adguiLegend_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to adguiLegend_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help adguiLegend

% Last Modified by GUIDE v2.5 02-May-2003 16:32:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @adguiLegend_OpeningFcn, ...
                   'gui_OutputFcn',  @adguiLegend_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before adguiLegend is made visible.
function adguiLegend_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to adguiLegend (see VARARGIN)

% Choose default command line output for adguiLegend
handles.output = hObject;

%make it a white figure
set(hObject,'Color',[1,1,1]);

%-----write pulldown-menu

%get adguiH
adguiH = findall(0,'Tag','adgui');
adguiHandles = guidata(adguiH);

%get plotH
plotFigureH = adguiHandles.plotFigureH;

%check whether there are any figures at all
if ~any(ishandle(plotFigureH))
    close(hObject)
    errordlg('You need at least one valid figure to display legend!');
    return
end

numberTitleStr = {'Figure No. '};
numFig = length(plotFigureH);

%create PD-String
pdString = cellstr([char(numberTitleStr(ones(numFig,1))),num2str([1:numFig]')]);

%set string
set(handles.adguiLeg_figure_PD,'String',pdString);

%find active figure
cFigH = gcf;
activeFigure = find(cFigH == plotFigureH);

if activeFigure %if there is a valid figure active: select for legend
    set(handles.adguiLeg_figure_PD,'Value',activeFigure);
else
    set(handles.adguiLeg_figure_PD,'Value',1);
end

%remember selection in PD-menu
handles.lastSelection = 1;

%store handles of axes and legend text
handles.axH = handles.adguiLeg_ax1;
handles.txtH = handles.adguiLeg_txt1;

handles.pdH = [handles.adguiLeg_figure_PD;handles.adguiLeg_pd_txt];

% Update handles structure
guidata(hObject, handles);

%draw legend items
adguiLeg_redrawLegend;


% --- Outputs from this function are returned to the command line.
function varargout = adguiLegend_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function adguiLeg_figure_PD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adguiLeg_figure_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in adguiLeg_figure_PD.
function adguiLeg_figure_PD_Callback(hObject, eventdata, handles)
% hObject    handle to adguiLeg_figure_PD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get selected entry and last selected entry in PD-String
selectedFigure = get(hObject,'Value');
lastSelection = handles.lastSelection;

%----------check that it is a valid figure handle

%get corresponding figure handles
adguiH = findall(0,'Tag','adgui');
adguiHandles = guidata(adguiH);

%get plotH
plotFigureH = adguiHandles.plotFigureH;

%prepare resetting figH
numberTitleStr = {'Figure No. '};
numFig = length(plotFigureH);

%if for some reason there are more legend entries than figures: reset pdString
if selectedFigure > numFig
    %create new PD-String
    pdString = cellstr([char(numberTitleStr(ones(numFig,1))),num2str([plotFigureH]')]);
    selectedFigure = numFig;
    if selectedFigure == 0
        close(handles.adguiLegend);
        errordlg('no valid figure - can''t display legend');
    else
        set(hObject,'String',pdString);
    end
end

%check whether there are any figures at all. If not, delete selected entry. If
%this was the only entry, delete legend (is this still necessary?)
while ~ishandle(plotFigureH(selectedFigure))
    set(hObject,'Value',1);
    pdString = get(hObject,'String');
    pdString(selectedFigure) = [];
    plotFigureH(selectedFigure) = [];
    
    %check whether this would result in an empty string
    if isempty(pdString)
        close(handles.adguiLegend);
        h = errordlg('no valid figure - can''t display legend');
        uiwait(h);
    else
        set(hObject,'String',pdString);
    end
    
    %reset selection (if last selection was)
    if handles.lastSelection < selectedFigure;
        selectedFigure = lastSelection;
        set(hObject,'Value',lastSelection);
    else
        selectedFigure = lastSelection-1;
        set(hObject,'Value',lastSelection-1);
    end
    if selectedFigure == 0
        selectedFigure = 1;
    end
    lastSelection = selectedFigure;
end

%remove double entries
cFigH = plotFigureH(selectedFigure);
plotFigureH = unique(plotFigureH);
selectedFigure = find(plotFigureH == cFigH);
set(hObject,'Value',selectedFigure);
numFig = length(plotFigureH);
pdString = cellstr([char(numberTitleStr(ones(numFig,1))),num2str([plotFigureH]')]);
set(hObject,'String',pdString);

%show figure
figure(plotFigureH(selectedFigure));

%store data
adguiHandles.plotFigureH = plotFigureH;
guidata(adguiH,adguiHandles);

handles.lastSelection = selectedFigure;
guidata(hObject,handles);

%update legend
adguiLeg_redrawLegend;





