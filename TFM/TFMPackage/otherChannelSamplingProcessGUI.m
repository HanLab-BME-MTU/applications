function varargout = otherChannelSamplingProcessGUI(varargin)
% otherChannelSamplingProcessGUI
% GUI for OtherChannelSamplingProcess (TFMPackage optional step).
%
% Fields supported:
%   - ChannelIndex (scalar)
%   - Measure: 'mean' or 'median'
%   - ComputeDFF0 (true/false)
%   - BaselineFrames (vector or range string)
%
% This GUI is created programmatically (no .fig needed) but follows
% GUIDE/gui_mainfcn conventions so it works with packageGUI "Set" button.
%
% Sangyoon Han / patched by ChatGPT

% ===== Begin initialization code (GUIDE-like) =====
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @otherChannelSamplingProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @otherChannelSamplingProcessGUI_OutputFcn, ...
                   'gui_LayoutFcn',  @otherChannelSamplingProcessGUI_LayoutFcn, ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% ===== End initialization code =====


% --- Executes just before otherChannelSamplingProcessGUI is made visible.
function otherChannelSamplingProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

% This handles embedding into packageGUI when called as:
% otherChannelSamplingProcessGUI('mainFig',handles.figure1,procID)
processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:});

% Load current process params
userData  = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

% ---- Defaults if missing ----
if ~isfield(funParams,'ChannelIndex') || isempty(funParams.ChannelIndex)
    funParams.ChannelIndex = 1;
end
if ~isfield(funParams,'Measure') || isempty(funParams.Measure)
    funParams.Measure = 'mean';
end
if ~isfield(funParams,'ComputeDFF0') || isempty(funParams.ComputeDFF0)
    funParams.ComputeDFF0 = false;
end
if ~isfield(funParams,'BaselineFrames') || isempty(funParams.BaselineFrames)
    funParams.BaselineFrames = 1;
end

% Populate UI
set(handles.edit_channelIndex, 'String', num2str(funParams.ChannelIndex));

m = lower(string(funParams.Measure));
if m == "median"
    set(handles.popup_measure, 'Value', 2);
else
    set(handles.popup_measure, 'Value', 1);
end

set(handles.checkbox_dff0, 'Value', logical(funParams.ComputeDFF0));

set(handles.edit_baselineFrames, 'String', framesToString(funParams.BaselineFrames));

% Enable/disable baseline entry based on dF/F0
toggleBaselineUI(handles);

% Update handles
handles.output = hObject;
set(hObject,'UserData',userData);
guidata(hObject, handles);


function varargout = otherChannelSamplingProcessGUI_OutputFcn(~, ~, handles)
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(~, ~, handles)

userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

% ---- ChannelIndex ----
chStr = strtrim(get(handles.edit_channelIndex,'String'));
chVal = str2double(chStr);
if ~isfinite(chVal) || chVal < 1 || mod(chVal,1) ~= 0
    errordlg('ChannelIndex must be a positive integer (e.g., 1).','Invalid input');
    return;
end
funParams.ChannelIndex = chVal;

% ---- Measure ----
measList = get(handles.popup_measure,'String');
measVal  = get(handles.popup_measure,'Value');
if iscell(measList)
    funParams.Measure = lower(strtrim(measList{measVal}));
else
    % older MATLAB might return char array
    funParams.Measure = lower(strtrim(measList(measVal,:)));
end
if ~ismember(funParams.Measure, {'mean','median'})
    funParams.Measure = 'mean';
end

% ---- dF/F0 ----
funParams.ComputeDFF0 = logical(get(handles.checkbox_dff0,'Value'));

% ---- Baseline frames ----
bfStr = strtrim(get(handles.edit_baselineFrames,'String'));
if funParams.ComputeDFF0
    try
        funParams.BaselineFrames = parseFrameSpec(bfStr);
        if isempty(funParams.BaselineFrames)
            error('Empty baseline frames.');
        end
    catch ME
        errordlg(sprintf('BaselineFrames is invalid.\n\nExamples:\n  1:10\n  1,2,3\n  5-12\n\nError: %s', ME.message), ...
            'Invalid BaselineFrames');
        return;
    end
else
    % Keep whatever, but normalize to something sane
    if isempty(bfStr)
        funParams.BaselineFrames = 1;
    else
        try
            funParams.BaselineFrames = parseFrameSpec(bfStr);
            if isempty(funParams.BaselineFrames), funParams.BaselineFrames = 1; end
        catch
            funParams.BaselineFrames = 1;
        end
    end
end

% ---- Apply ----
try
    userData.crtProc.setPara(funParams);
    userData.crtProc.sanityCheck();
catch ME
    errordlg(sprintf('Failed to apply parameters:\n\n%s', ME.message), 'Process parameter error');
    return;
end

set(handles.figure1,'UserData',userData);
guidata(handles.figure1,handles);

delete(handles.figure1);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(~, ~, handles)
delete(handles.figure1);


% --- Executes on checkbox toggle
function checkbox_dff0_Callback(~, ~, handles)
toggleBaselineUI(handles);


function toggleBaselineUI(handles)
useDFF0 = logical(get(handles.checkbox_dff0,'Value'));
if useDFF0
    set(handles.edit_baselineFrames,'Enable','on');
    set(handles.text_baselineFrames,'Enable','on');
else
    set(handles.edit_baselineFrames,'Enable','off');
    set(handles.text_baselineFrames,'Enable','off');
end


% ===== Layout function (programmatic; no .fig needed) =====
function h1 = otherChannelSamplingProcessGUI_LayoutFcn(varargin)

% Create figure
h1 = figure('Units','pixels', ...
    'Position',[200 200 420 240], ...
    'MenuBar','none', ...
    'Name','Other Channel Sampling', ...
    'NumberTitle','off', ...
    'Color',get(0,'defaultUicontrolBackgroundColor'), ...
    'Tag','figure1', ...
    'Visible','off', ...
    'HandleVisibility','callback', ...
    'IntegerHandle','off');

% GUIDEOptions to avoid syscolorfig errors in some MATLAB/GUI stacks
try
    setappdata(h1,'GUIDEOptions',struct( ...
        'active_h',         [], ...
        'taginfo',          struct(), ...
        'override',         0, ...
        'release',          13, ...
        'resize',           'none', ...
        'accessibility',    'on', ...
        'syscolorfig',      1)); %#ok<STRNU>
catch
end

% ---- Controls ----
uicontrol('Parent',h1,'Style','text', ...
    'Units','pixels','Position',[20 190 120 18], ...
    'String','ChannelIndex','HorizontalAlignment','left');

uicontrol('Parent',h1,'Style','edit', ...
    'Units','pixels','Position',[150 188 80 24], ...
    'String','1','BackgroundColor','white', ...
    'Tag','edit_channelIndex');

uicontrol('Parent',h1,'Style','text', ...
    'Units','pixels','Position',[20 155 120 18], ...
    'String','Measure','HorizontalAlignment','left');

uicontrol('Parent',h1,'Style','popupmenu', ...
    'Units','pixels','Position',[150 152 120 26], ...
    'String',{'mean','median'}, ...
    'BackgroundColor','white', ...
    'Tag','popup_measure');

uicontrol('Parent',h1,'Style','checkbox', ...
    'Units','pixels','Position',[20 118 220 20], ...
    'String','Compute dF/F0', ...
    'Value',0, ...
    'Tag','checkbox_dff0', ...
    'Callback',@checkbox_dff0_Callback);

uicontrol('Parent',h1,'Style','text', ...
    'Units','pixels','Position',[20 85 120 18], ...
    'String','BaselineFrames','HorizontalAlignment','left', ...
    'Tag','text_baselineFrames');

uicontrol('Parent',h1,'Style','edit', ...
    'Units','pixels','Position',[150 82 180 24], ...
    'String','1:10', ...
    'BackgroundColor','white', ...
    'Tag','edit_baselineFrames');

uicontrol('Parent',h1,'Style','pushbutton', ...
    'Units','pixels','Position',[210 20 90 30], ...
    'String','Done', ...
    'Tag','pushbutton_done', ...
    'Callback',@pushbutton_done_Callback);

uicontrol('Parent',h1,'Style','pushbutton', ...
    'Units','pixels','Position',[310 20 90 30], ...
    'String','Cancel', ...
    'Tag','pushbutton_cancel', ...
    'Callback',@pushbutton_cancel_Callback);

% Create handles structure like GUIDE
handles = guihandles(h1);
guidata(h1, handles);

set(h1,'Visible','on');


% ===== Helpers =====
function s = framesToString(frames)
if isnumeric(frames)
    frames = frames(:)';
    if numel(frames) >= 2 && all(diff(frames)==1)
        s = sprintf('%d:%d', frames(1), frames(end));
    else
        s = strtrim(sprintf('%d,', frames));
        if ~isempty(s), s(end) = []; end
    end
else
    s = char(string(frames));
end

function v = parseFrameSpec(spec)
% Accept:
%   "1:10"
%   "5-12"
%   "1,2,3"
%   "1 2 3"
spec = strtrim(spec);
if isempty(spec)
    v = [];
    return;
end

% normalize hyphen range
spec = regexprep(spec,'(\d)\s*-\s*(\d)','$1:$2');

% allow spaces as commas
spec = regexprep(spec,'\s+',',');

% safety: only allow digits, commas, colons
if ~isempty(regexp(spec,'[^0-9,:]','once'))
    error('BaselineFrames contains invalid characters.');
end

% Evaluate safely
parts = split(string(spec), ",");
v = [];
for i=1:numel(parts)
    p = strtrim(parts(i));
    if p == "", continue; end
    if contains(p, ":")
        rr = split(p, ":");
        if numel(rr) ~= 2, error('Invalid range: %s', p); end
        a = str2double(rr(1)); b = str2double(rr(2));
        if ~isfinite(a) || ~isfinite(b) || a<1 || b<1 || mod(a,1)~=0 || mod(b,1)~=0
            error('Invalid range endpoints: %s', p);
        end
        if b < a
            v = [v, a:-1:b]; %#ok<AGROW>
        else
            v = [v, a:b]; %#ok<AGROW>
        end
    else
        a = str2double(p);
        if ~isfinite(a) || a<1 || mod(a,1)~=0
            error('Invalid frame index: %s', p);
        end
        v = [v, a]; %#ok<AGROW>
    end
end
v = unique(v(:)','stable');