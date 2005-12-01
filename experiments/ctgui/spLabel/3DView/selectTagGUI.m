function varargout = selectTagGUI(varargin)
% SELECTTAGGUI M-file for selectTagGUI.fig
%      SELECTTAGGUI, by itself, creates a new SELECTTAGGUI or raises the existing
%      singleton*.
%
%      H = SELECTTAGGUI returns the handle to a new SELECTTAGGUI or the handle to
%      the existing singleton*.
%
%      SELECTTAGGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECTTAGGUI.M with the given input arguments.
%
%      SELECTTAGGUI('Property','Value',...) creates a new SELECTTAGGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before selectTagGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to selectTagGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help selectTagGUI

% Last Modified by GUIDE v2.5 19-Mar-2003 10:53:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @selectTagGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @selectTagGUI_OutputFcn, ...
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


% --- Executes just before selectTagGUI is made visible.
function selectTagGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to selectTagGUI (see VARARGIN)

% Choose default command line output for selectTagGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%------------------initialize-----------------
%get view3D-data
labelPanelH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
view3DH = GetUserData(labelPanelH,'view3DGUIH');
view3D_handles = guidata(view3DH);
nTags = view3D_handles.nTags;
cMap = view3D_handles.cMap;
cMapFact = view3D_handles.cMapFact;
labelList = view3D_handles.labelList;

setTagGUIH = hObject;

%------------------draw full GUI-----------------
%adjust GUI size
if ~strcmp(get(setTagGUIH,'Units'),'pixels')
    set(setTagGUIH,'Units','pixels');
end

pos = get(setTagGUIH,'Position');
pos(4) = 50+nTags*50;
set(setTagGUIH,'Position',pos);

%draw buttons and text
for i = 1:nTags
    buttonColor = cMap(2^(i-1)*cMapFact,:);
    buttonH(i) = uicontrol('Style','togglebutton','BackgroundColor',buttonColor,...
        'Tag',['Button_',num2str(i)],...
        'Position',[20,60+(i-1)*50,30,30],...
        'Callback','selectTagGUI(''select_TB_CB'',gcbo,[],guidata(gcbo))');
    txtH(i) = uicontrol('Style','text',...
        'FontSize',12,'HorizontalAlignment','left',...
        'Position',[90,70+(i-1)*50,110,20],...
        'String',char(labelList(i)));
end

handles.buttonH = buttonH;
handles.txtH = txtH;

guidata(hObject,handles);


% --- Outputs from this function are returned to the command line.
function varargout = selectTagGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in selectTag_OK_PB.
function selectTag_OK_PB_Callback(hObject, eventdata, handles)
% hObject    handle to selectTag_OK_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

buttonValueC = get(handles.buttonH,'Value');
buttonValue = cat(1,buttonValueC{1:end});

switch handles.num2select
    case 1 %one tag has to be selected
        %test if enough buttons
        if sum(buttonValue)<1
            h = warndlg('Select one tag!','Input not acceptable!');
            uiwait(h)
            return %end evaluation here
        end
        view3DH = findall(0,'Tag','view3DGUI');
        view3D_handles = guidata(view3DH);
        view3D_handles.centerTagNumber = find(buttonValue);
        guidata(view3DH,view3D_handles);
    case 2 %two tags have to be selected
        %test if enough buttons
        if sum(buttonValue)<2
            h = warndlg('Select two tags!','Input not acceptable!');
            uiwait(h)
            return %end evaluation here
        end
        view3DH = findall(0,'Tag','view3DGUI');
        view3D_handles = guidata(view3DH);
        view3D_handles.axisTagNumber = find(buttonValue==1);
        view3D_handles.nonAxisTagNumber = find(buttonValue==0);
        guidata(view3DH,view3D_handles);
end

close(handles.STGUI);


% --- Executes on button press in selectTag_cancel_PB.
function selectTag_cancel_PB_Callback(hObject, eventdata, handles)
% hObject    handle to selectTag_cancel_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

%confirm user input
really = questdlg('Do you really want to quit without saving changes?','Quit?','yes','no','yes');
if strcmp(really,'no')
    return %end evaluation here
end

labelPanelH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
view3DH = GetUserData(labelPanelH,'view3DGUIH');

switch handles.num2select
    
    case 1 %one tag has to be selected
        view3D_handles = guidata(view3DH);
        view3D_handles.centerTagNumber = 0;
        guidata(view3DH,view3D_handles);
    case 2 %two tags have to be selected
        view3D_handles = guidata(view3DH);
        view3D_handles.axisTagNumber = 0;
        guidata(view3DH,view3D_handles);
end

close(handles.STGUI);


function select_TB_CB(hObject, eventdata, handles)

%test if being selected, if not, do nothing
if get(hObject,'Value')==0
    return %do nothing
end

%get button number
myButtonStr = get(hObject,'Tag');
myButtonNum = str2num(myButtonStr(8:end));

%switch how many to select
switch handles.num2select
    case 1 %one tag has to be selected
        set(handles.buttonH,'Value',0);
        set(hObject,'Value',1);
    case 2
        buttonValueC = get(handles.buttonH,'Value');
        buttonValue = cat(1,buttonValueC{1:end});
        
        if sum(buttonValue)>2
            %too many buttons are selected: select only current
            set(handles.buttonH,'Value',0);
            set(hObject,'Value',1);
        end
end

%just to make sure
guidata(handles.STGUI,handles);