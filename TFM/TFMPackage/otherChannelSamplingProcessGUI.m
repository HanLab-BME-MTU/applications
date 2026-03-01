function varargout = otherChannelSamplingProcessGUI(varargin)
% otherChannelSamplingProcessGUI
% GUI for OtherChannelSamplingProcess
%
% This GUI is implemented programmatically (no .fig) to avoid GUIDE/FIG
% compatibility issues across MATLAB versions.
%
% Expected invocation pattern from packageGUI:
%   otherChannelSamplingProcessGUI('mainFig', hMain, 'procID', procID)
%
% Sangyoon Han lab, 2026

% --- Standard GUIDE-style state struct so packageGUI can call it ---

gui_Singleton = 1;
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
end

%% ===================== LAYOUT =====================
function hFig = otherChannelSamplingProcessGUI_LayoutFcn()

hFig = figure('Units','pixels', ...
    'Position',[300 200 620 360], ...
    'MenuBar','none', ...
    'ToolBar','none', ...
    'NumberTitle','off', ...
    'Name','Other channel sampling', ...
    'Color','w', ...
    'Tag','figure1', ...
    'Visible','off');

% Title
uicontrol('Parent',hFig,'Style','text','String','Other channel sampling (optional Step 6)', ...
    'Units','pixels','Position',[20 325 580 22], 'HorizontalAlignment','left', ...
    'FontWeight','bold','BackgroundColor','w','Tag','text_processName');

% Output dir
uicontrol('Parent',hFig,'Style','text','String','Output directory:', ...
    'Units','pixels','Position',[20 285 120 18],'HorizontalAlignment','left','BackgroundColor','w');

uicontrol('Parent',hFig,'Style','edit','String','', ...
    'Units','pixels','Position',[150 280 360 26],'HorizontalAlignment','left', ...
    'BackgroundColor','white','Tag','edit_outputDir');

uicontrol('Parent',hFig,'Style','pushbutton','String','Browse...', ...
    'Units','pixels','Position',[520 280 80 26],'Tag','pushbutton_browse', ...
    'Callback',@pushbutton_browse_Callback);

% Channel index
uicontrol('Parent',hFig,'Style','text','String','Other channel index:', ...
    'Units','pixels','Position',[20 240 140 18],'HorizontalAlignment','left','BackgroundColor','w');

uicontrol('Parent',hFig,'Style','edit','String','1', ...
    'Units','pixels','Position',[170 235 60 26],'BackgroundColor','white','Tag','edit_channelIndex');

% Mask source
uicontrol('Parent',hFig,'Style','text','String','Mask source:', ...
    'Units','pixels','Position',[20 200 140 18],'HorizontalAlignment','left','BackgroundColor','w');

uicontrol('Parent',hFig,'Style','popupmenu','String',{'MaskRefinementProcess','MaskIntersectionProcess','None (whole FOV)'}, ...
    'Units','pixels','Position',[170 195 220 26],'BackgroundColor','white','Tag','popup_maskSource');

uicontrol('Parent',hFig,'Style','text','String','Mask channel index:', ...
    'Units','pixels','Position',[410 200 120 18],'HorizontalAlignment','left','BackgroundColor','w');

uicontrol('Parent',hFig,'Style','edit','String','1', ...
    'Units','pixels','Position',[535 195 60 26],'BackgroundColor','white','Tag','edit_maskChan');

% Statistics
uicontrol('Parent',hFig,'Style','text','String','Statistic:', ...
    'Units','pixels','Position',[20 160 140 18],'HorizontalAlignment','left','BackgroundColor','w');

uicontrol('Parent',hFig,'Style','popupmenu','String',{'mean','median','sum'}, ...
    'Units','pixels','Position',[170 155 120 26],'BackgroundColor','white','Tag','popup_stat');

% dF/F0
uicontrol('Parent',hFig,'Style','checkbox','String','Compute dF/F0 (per ROI)', ...
    'Units','pixels','Position',[20 118 220 22],'BackgroundColor','w','Tag','checkbox_doDFF0');

uicontrol('Parent',hFig,'Style','text','String','Baseline frames (e.g., 1:10):', ...
    'Units','pixels','Position',[260 120 200 18],'HorizontalAlignment','left','BackgroundColor','w');

uicontrol('Parent',hFig,'Style','edit','String','1:10', ...
    'Units','pixels','Position',[465 115 130 26],'BackgroundColor','white','Tag','edit_baselineFrames');

% Buttons
uicontrol('Parent',hFig,'Style','pushbutton','String','Done', ...
    'Units','pixels','Position',[420 20 80 30],'Tag','pushbutton_done', ...
    'Callback',@pushbutton_done_Callback);

uicontrol('Parent',hFig,'Style','pushbutton','String','Cancel', ...
    'Units','pixels','Position',[515 20 80 30],'Tag','pushbutton_cancel', ...
    'Callback',@pushbutton_cancel_Callback);

% Label expected by processGUI helpers
uicontrol('Parent',hFig,'Style','text','String','', ...
    'Units','pixels','Position',[20 15 360 18],'HorizontalAlignment','left','BackgroundColor','w','Tag','text_copyright');

end

%% ===================== OPENING =====================
function otherChannelSamplingProcessGUI_OpeningFcn(hObject, ~, handles, varargin)

% Wire up the standard MovieManagement userData/crtProc machinery.
processGUI_OpeningFcn(hObject, [], handles, varargin{:}, ...
    'initChannel', 1, ...
    'initProcess', @OtherChannelSamplingProcess);

handles = guidata(hObject);
proc = handles.userData.crtProc;
fp = proc.funParams_;

% Populate fields
setSafe(handles,'edit_outputDir', fp.OutputDirectory);
setSafe(handles,'edit_channelIndex', num2str(fp.ChannelIndex));

% Mask selection
maskChoices = getPopupChoices(handles,'popup_maskSource');
if ~isempty(fp.MaskProcessName)
    idx = find(strcmpi(maskChoices, fp.MaskProcessName), 1);
    if isempty(idx), idx = 1; end
else
    idx = 3; % None
end
set(handles.popup_maskSource,'Value',idx);

setSafe(handles,'edit_maskChan', num2str(fp.MaskChannelIndex));

% Statistic
statChoices = getPopupChoices(handles,'popup_stat');
idx = find(strcmpi(statChoices, fp.Statistic), 1);
if isempty(idx), idx = 1; end
set(handles.popup_stat,'Value',idx);

% dF/F0
set(handles.checkbox_doDFF0,'Value', logical(fp.DoDFF0));

% Baseline frames
bfStr = '1:10';
if isnumeric(fp.BaselineFrames) && ~isempty(fp.BaselineFrames)
    bf = fp.BaselineFrames(:)';
    if numel(bf)>=2 && all(diff(bf)==1)
        bfStr = sprintf('%d:%d', bf(1), bf(end));
    else
        bfStr = mat2str(bf);
    end
end
setSafe(handles,'edit_baselineFrames', bfStr);

set(hObject,'Visible','on');

guidata(hObject, handles);
end

%% ===================== OUTPUT =====================
function varargout = otherChannelSamplingProcessGUI_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
end

%% ===================== CALLBACKS =====================
function pushbutton_browse_Callback(hObject, ~)
handles = guidata(hObject);
startDir = get(handles.edit_outputDir,'String');
if isempty(startDir) || exist(startDir,'dir')~=7
    startDir = pwd;
end
p = uigetdir(startDir, 'Select output directory');
if isequal(p,0), return; end
set(handles.edit_outputDir,'String',p);
end

function pushbutton_done_Callback(hObject, ~)
handles = guidata(hObject);
proc = handles.userData.crtProc;
fp = proc.funParams_;

% Output directory
outDir = strtrim(get(handles.edit_outputDir,'String'));
if isempty(outDir)
    outDir = fp.OutputDirectory;
end
fp.OutputDirectory = outDir;

% Channel index
ch = str2double(get(handles.edit_channelIndex,'String'));
if ~isfinite(ch) || ch < 1
    errordlg('Other channel index must be a positive integer.','Invalid input');
    return;
end
fp.ChannelIndex = round(ch);

% Mask source
maskChoices = getPopupChoices(handles,'popup_maskSource');
maskSel = maskChoices{get(handles.popup_maskSource,'Value')};
if contains(lower(maskSel),'none')
    fp.MaskProcessName = '';
else
    fp.MaskProcessName = maskSel;
end

% Mask channel
mc = str2double(get(handles.edit_maskChan,'String'));
if ~isfinite(mc) || mc < 1
    mc = fp.MaskChannelIndex;
end
fp.MaskChannelIndex = round(mc);

% Statistic
statChoices = getPopupChoices(handles,'popup_stat');
fp.Statistic = statChoices{get(handles.popup_stat,'Value')};

% dF/F0 + baseline frames
fp.DoDFF0 = logical(get(handles.checkbox_doDFF0,'Value'));

bfStr = strtrim(get(handles.edit_baselineFrames,'String'));
if isempty(bfStr)
    bf = fp.BaselineFrames;
else
    try
        bf = eval(bfStr); %#ok<EVLDIR>
    catch
        errordlg('Baseline frames must be a valid MATLAB expression (e.g., 1:10).','Invalid input');
        return;
    end
end
if ~isnumeric(bf) || isempty(bf)
    errordlg('Baseline frames must evaluate to a numeric vector.','Invalid input');
    return;
end
fp.BaselineFrames = unique(round(bf(:)'));

% Apply (writes back into MD + refresh main GUI)
proc.setPara(fp);
processGUI_ApplyFcn(hObject, [], handles, fp);

delete(handles.figure1);
end

function pushbutton_cancel_Callback(hObject, ~)
handles = guidata(hObject);
delete(handles.figure1);
end

%% ===================== small helpers =====================
function setSafe(handles, field, val)
if isfield(handles, field)
    if isempty(val), val = ''; end
    set(handles.(field), 'String', char(string(val)));
end
end

function choices = getPopupChoices(handles, field)
raw = get(handles.(field), 'String');
if ischar(raw)
    choices = cellstr(raw);
else
    choices = raw;
end
end
