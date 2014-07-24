function varargout = inputGUI(varargin)
% INPUTGUI M-file for inputGUI.fig
%
% [strList,useAnd]=inputGUI to get cell array of strings corresponding to
% user input. useAnd will be 1 if 'and' 0 if 'or'.
%
%
%      INPUTGUI, by itself, creates a new INPUTGUI or raises the existing
%      singleton*.
%
%      H = INPUTGUI returns the handle to a new INPUTGUI or the handle to
%      the existing singleton*.
%
%      INPUTGUI('CALLBACK',hObject,eventData,h,...) calls the local
%      function named CALLBACK in INPUTGUI.M with the given input arguments.
%
%      INPUTGUI('Property','Value',...) creates a new INPUTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before inputGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to inputGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help inputGUI

% Last Modified by GUIDE v2.5 01-May-2009 13:00:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @inputGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @inputGUI_OutputFcn, ...
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


% --- Executes just before inputGUI is made visible.
function inputGUI_OpeningFcn(hObject, eventdata, h, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
% varargin   command line arguments to inputGUI (see VARARGIN)

% Choose default command line output for inputGUI

h.query=[];
h.inputQuery=[];

h.output = hObject;

% Update h structure
guidata(hObject, h);

% UIWAIT makes inputGUI wait for user response (see UIRESUME)
uiwait(h.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = inputGUI_OutputFcn(hObject, eventdata, h) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Get default command line output from h structure
if nargin>=1
    if ~isempty(h);
        varargout{1} = h.inputQuery;
        close(h.figure1)
    else
        varargout{1} = [];
    end
end



function addQuery_pushbutton_Callback(hObject, eventdata, h)
%gets input file(s) from user
input_file=h.query;
 
%gets the current data file names inside the listbox
inputQuery = get(h.inputFiles_listbox,'String');
inputQuery{end+1} = input_file;

 
%updates the gui to display all filenames in the listbox
set(h.inputFiles_listbox,'String',inputQuery);
 
%make sure first file is always selected so it doesn't go out of range
%the GUI will break if this value is out of range
set(h.inputFiles_listbox,'Value',1);
 
% Update h structure
guidata(hObject, h);


function deleteQuery_pushbutton_Callback(hObject, eventdata, h)
%get the current list of file names from the listbox
inputQuery = get(h.inputFiles_listbox,'String');
 
%get the values for the selected file names
option = get(h.inputFiles_listbox,'Value');
 
%is there is nothing to delete, nothing happens
if (isempty(option) == 1 || option(1) == 0 || isempty(inputQuery))
    return
end
 
%erases the contents of highlighted item in data array
inputQuery(option) = [];
 
%updates the gui, erasing the selected item from the listbox
set(h.inputFiles_listbox,'String',inputQuery);
 
%moves the highlighted item to an appropiate value or else will get error
if option(end) > length(inputQuery)
    set(h.inputFiles_listbox,'Value',length(inputQuery));
end
 
% Update h structure
guidata(hObject, h);


function inputFiles_listbox_Callback(hObject, eventdata, h)
%nothing needed here.  we are just using the listbox to display data!


function closeGui_pushbutton_Callback(hObject, eventdata, h)
%get the current list of file names from the listbox
inputQuery = get(h.inputFiles_listbox,'String');
if iscell(inputQuery)
h.inputQuery=inputQuery(~cellfun(@isempty, inputQuery),1);
else
    h.inputQuery=inputQuery;
end
guidata(hObject, h);
uiresume(h.figure1);


function inputFiles_listbox_CreateFcn(hObject, eventdata, h)

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function strEditBox_Callback(hObject, eventdata, h)
% hObject    handle to strEditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of strEditBox as text
%        str2double(get(hObject,'String')) returns contents of strEditBox as a double
h.query=get(hObject,'String');
guidata(hObject, h);


% --- Executes during object creation, after setting all properties.
function strEditBox_CreateFcn(hObject, eventdata, h)
% hObject    handle to strEditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




