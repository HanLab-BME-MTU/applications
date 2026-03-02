function varargout = otherChannelSamplingProcessGUI(varargin)
% otherChannelSamplingProcessGUI M-file for otherChannelSamplingProcessGUI.fig
%
% Settings GUI for OtherChannelSamplingProcess.
% This GUI follows the same structure as subResolutionProcessGUI.m:
%   - Uses processGUI_OpeningFcn with 'initChannel'
%   - Uses standard channel selection listboxes from the common processGUI template
%   - Adds analysis-specific parameter controls dynamically inside the Parameters panel.
%
% 2026 - Sangyoon Han lab

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @otherChannelSamplingProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @otherChannelSamplingProcessGUI_OutputFcn, ...
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
end


% --- Executes just before otherChannelSamplingProcessGUI is made visible.
function otherChannelSamplingProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

% Initialize the standard process GUI (channel listboxes, etc.)
processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:}, 'initChannel', 1);

% Refresh handles after processGUI_OpeningFcn
handles = guidata(hObject);

% Parameter Setup
userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end
funParams = userData.crtProc.funParams_;

% Title
safeSetString(handles, 'text_processName', 'OtherChannelSamplingProcess settings');
set(handles.figure1, 'Name', 'Setting - Other Channel Sampling');

% --------------------
% Build parameter panel
% --------------------
handles = ocsp_buildParamControls(handles);
handles = ocsp_initParamControls(handles, funParams, userData);

% Choose default command line output
handles.output = hObject;

set(hObject, 'UserData', userData);
uicontrol(handles.pushbutton_done);
guidata(hObject, handles);
end


% --- Outputs from this function are returned to the command line.
function varargout = otherChannelSamplingProcessGUI_OutputFcn(hObject, eventdata, handles)
% Return figure handle like most process GUIs
varargout{1} = handles.figure1;
end


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% Close without applying
try
    delete(handles.figure1);
catch
    try, close(handles.figure1); end %#ok<TRYNC>
end
end


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% Apply button
userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end

% Must select at least one sampling channel
if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Input Channels''.', 'Setting Error', 'modal');
    return;
end
channelIndex = get(handles.listbox_selectedChannels, 'Userdata');
if isempty(channelIndex)
    % Fallback: if UserData not set, use listbox value
    channelIndex = get(handles.listbox_selectedChannels, 'Value');
end

funParams = userData.crtProc.funParams_;
funParams.ChannelIndex = channelIndex;

% -------- Mask source --------
if isfield(handles, 'ocsp_popup_maskProcess')
    ud = get(handles.ocsp_popup_maskProcess, 'UserData');
    v  = get(handles.ocsp_popup_maskProcess, 'Value');
    if iscell(ud) && v>=1 && v<=numel(ud)
        funParams.MaskProcessName = ud{v};
    end
end

if isfield(handles, 'ocsp_edit_maskChan')
    maskChan = str2double(get(handles.ocsp_edit_maskChan,'String'));
    if isnan(maskChan) || maskChan < 1
        errordlg('Mask channel index must be a positive integer.', 'Setting Error', 'modal');
        return;
    end
    funParams.MaskChannelIndex = round(maskChan);
end

% -------- Options --------
funParams.UseStageDriftCorrection = logical(get(handles.ocsp_check_useSDC,'Value'));
funParams.SavePerFrameTifPreview  = logical(get(handles.ocsp_check_preview,'Value'));
funParams.UseLabeling            = logical(get(handles.ocsp_check_useLabeling,'Value'));
funParams.ComputeDFF0            = logical(get(handles.ocsp_check_computeDFF0,'Value'));

minArea = str2double(get(handles.ocsp_edit_minArea,'String'));
if isnan(minArea) || minArea < 0
    errordlg('Min object area (pixels) must be a nonnegative number.', 'Setting Error', 'modal');
    return;
end
funParams.MinAreaPix = round(minArea);

% Baseline frames
bfStr = strtrim(get(handles.ocsp_edit_baselineFrames,'String'));
if isempty(bfStr)
    bf = 1;
else
    bf = str2num(bfStr); %#ok<ST2NM>
end
if isempty(bf) || ~isnumeric(bf)
    errordlg('Baseline frames must be numeric (e.g., 1:10 or [1 2 3]).', 'Setting Error', 'modal');
    return;
end
funParams.BaselineFrames = bf(:)';

bfN = str2double(get(handles.ocsp_edit_baselineFallbackN,'String'));
if isnan(bfN) || bfN < 0
    errordlg('Baseline fallback N must be >= 0.', 'Setting Error', 'modal');
    return;
end
funParams.BaselineFallbackN = round(bfN);

% -------- Tracking --------
funParams.TrackCells = logical(get(handles.ocsp_check_trackCells,'Value'));
funParams.TrackByOverlap = logical(get(handles.ocsp_check_trackByOverlap,'Value'));

minIoU = str2double(get(handles.ocsp_edit_minIoU,'String'));
if isnan(minIoU) || minIoU < 0 || minIoU > 1
    errordlg('Min IoU must be between 0 and 1.', 'Setting Error', 'modal');
    return;
end
funParams.MinIoU = minIoU;

maxDist = str2double(get(handles.ocsp_edit_maxCentroidDist,'String'));
if isnan(maxDist) || maxDist <= 0
    errordlg('Max centroid distance must be a positive number.', 'Setting Error', 'modal');
    return;
end
funParams.MaxCentroidDist = maxDist;

% -------- Instance split --------
ud = get(handles.ocsp_popup_instanceMethod,'UserData');
v  = get(handles.ocsp_popup_instanceMethod,'Value');
if iscell(ud) && v>=1 && v<=numel(ud)
    funParams.InstanceSegMethod = ud{v};
end

sig = str2double(get(handles.ocsp_edit_smoothSigma,'String'));
if isnan(sig) || sig < 0
    errordlg('SmoothSigma must be >= 0.', 'Setting Error', 'modal');
    return;
end
funParams.SmoothSigma = sig;

hval = str2double(get(handles.ocsp_edit_minSeedH,'String'));
if isnan(hval) || hval < 0
    errordlg('MinSeedH must be >= 0.', 'Setting Error', 'modal');
    return;
end
funParams.MinSeedH = hval;

funParams.BorderClear = logical(get(handles.ocsp_check_borderClear,'Value'));

% -------- Progress bar --------
funParams.UseProgressBar = logical(get(handles.ocsp_check_useProgressBar,'Value'));
ud2 = get(handles.ocsp_popup_progressPos,'UserData');
v2  = get(handles.ocsp_popup_progressPos,'Value');
if isnumeric(ud2) && v2>=1 && v2<=numel(ud2)
    funParams.ProgressBarPosition = ud2(v2);
end

% Sanity check
try
    userData.crtProc.sanityCheck;
catch ME
    errordlg([ME.message ' Please double check your data.'], 'Setting Error', 'modal');
    return;
end

% Apply
processGUI_ApplyFcn(hObject, eventdata, handles, funParams);
end


% =========================
% Build + init param controls
% =========================

function handles = ocsp_buildParamControls(handles)
% Create UI controls inside the Parameters panel without modifying the .fig

% Find parameter panel
if isfield(handles,'uipanel_1')
    parent = handles.uipanel_1;
else
    parent = handles.figure1;
end

set(parent,'Units','pixels');
pPos = get(parent,'Position');

% Clear any previously created OCSP controls (if GUI reopened)
ch = allchild(parent);
for i=1:numel(ch)
    tg = '';
    try tg = get(ch(i),'Tag'); end %#ok<TRYNC>
    if startsWith(tg,'ocsp_')
        delete(ch(i));
    end
end

% Font sizing (Linux tends to scale GUIDE fonts huge)
baseFont = 10;
if isunix && ~ismac
    baseFont = 9;
end

% Layout grid
margin = 10;
colW = (pPos(3) - 3*margin) / 2;
rowH = 18;
lineGap = 4;

xL = margin;
xR = margin*2 + colW;
yTop = pPos(4) - margin;

% Section titles
handles.ocsp_text_tracking = uicontrol(parent,'Style','text','Tag','ocsp_text_tracking',...
    'String','Tracking','Units','pixels','HorizontalAlignment','left',...
    'FontWeight','bold','FontSize',baseFont+2,'Position',[xL yTop-20 colW 18]);

handles.ocsp_text_instance = uicontrol(parent,'Style','text','Tag','ocsp_text_instance',...
    'String','Instance split','Units','pixels','HorizontalAlignment','left',...
    'FontWeight','bold','FontSize',baseFont+2,'Position',[xR yTop-20 colW 18]);

% --- Left column: Tracking + dF/F0 + Mask/Options
y = yTop - 42;

handles.ocsp_check_trackCells = uicontrol(parent,'Style','checkbox','Tag','ocsp_check_trackCells',...
    'String','Track cells across frames','Units','pixels','FontSize',baseFont,...
    'Position',[xL y colW 18]);

y = y - (rowH+lineGap);
handles.ocsp_check_trackByOverlap = uicontrol(parent,'Style','checkbox','Tag','ocsp_check_trackByOverlap',...
    'String','Link by overlap (IoU)','Units','pixels','FontSize',baseFont,...
    'Position',[xL y colW 18]);

y = y - (rowH+lineGap);
handles.ocsp_text_minIoU = uicontrol(parent,'Style','text','Tag','ocsp_text_minIoU',...
    'String','Min IoU','Units','pixels','HorizontalAlignment','left','FontSize',baseFont,...
    'Position',[xL y 80 18]);
handles.ocsp_edit_minIoU = uicontrol(parent,'Style','edit','Tag','ocsp_edit_minIoU',...
    'String','0.10','Units','pixels','FontSize',baseFont,...
    'Position',[xL+90 y 80 20]);

y = y - (rowH+lineGap);
handles.ocsp_text_maxDist = uicontrol(parent,'Style','text','Tag','ocsp_text_maxDist',...
    'String','Max centroid dist (px)','Units','pixels','HorizontalAlignment','left','FontSize',baseFont,...
    'Position',[xL y 140 18]);
handles.ocsp_edit_maxCentroidDist = uicontrol(parent,'Style','edit','Tag','ocsp_edit_maxCentroidDist',...
    'String','40','Units','pixels','FontSize',baseFont,...
    'Position',[xL+150 y 60 20]);

% dF/F0 section
y = y - 28;
handles.ocsp_text_dff0 = uicontrol(parent,'Style','text','Tag','ocsp_text_dff0',...
    'String','dF/F0','Units','pixels','HorizontalAlignment','left',...
    'FontWeight','bold','FontSize',baseFont+2,'Position',[xL y colW 18]);

y = y - (rowH+lineGap);
handles.ocsp_check_computeDFF0 = uicontrol(parent,'Style','checkbox','Tag','ocsp_check_computeDFF0',...
    'String','Compute dF/F0 (whole + per-track)','Units','pixels','FontSize',baseFont,...
    'Position',[xL y colW 18]);

y = y - (rowH+lineGap);
handles.ocsp_text_baselineFrames = uicontrol(parent,'Style','text','Tag','ocsp_text_baselineFrames',...
    'String','Baseline frames','Units','pixels','HorizontalAlignment','left','FontSize',baseFont,...
    'Position',[xL y 110 18]);
handles.ocsp_edit_baselineFrames = uicontrol(parent,'Style','edit','Tag','ocsp_edit_baselineFrames',...
    'String','1:10','Units','pixels','FontSize',baseFont,...
    'Position',[xL+120 y 120 20]);

y = y - (rowH+lineGap);
handles.ocsp_text_baselineFallbackN = uicontrol(parent,'Style','text','Tag','ocsp_text_baselineFallbackN',...
    'String','Fallback N','Units','pixels','HorizontalAlignment','left','FontSize',baseFont,...
    'Position',[xL y 110 18]);
handles.ocsp_edit_baselineFallbackN = uicontrol(parent,'Style','edit','Tag','ocsp_edit_baselineFallbackN',...
    'String','5','Units','pixels','FontSize',baseFont,...
    'Position',[xL+120 y 60 20]);

% Mask & Options section
y = y - 32;
handles.ocsp_text_mask = uicontrol(parent,'Style','text','Tag','ocsp_text_mask',...
    'String','Mask / Options','Units','pixels','HorizontalAlignment','left',...
    'FontWeight','bold','FontSize',baseFont+2,'Position',[xL y colW 18]);

y = y - (rowH+lineGap);
handles.ocsp_text_maskProc = uicontrol(parent,'Style','text','Tag','ocsp_text_maskProc',...
    'String','Mask process','Units','pixels','HorizontalAlignment','left','FontSize',baseFont,...
    'Position',[xL y 110 18]);
handles.ocsp_popup_maskProcess = uicontrol(parent,'Style','popupmenu','Tag','ocsp_popup_maskProcess',...
    'String',{'(auto-detect)'},'Units','pixels','FontSize',baseFont,...
    'Position',[xL+120 y 160 22]);

y = y - (rowH+lineGap);
handles.ocsp_text_maskChan = uicontrol(parent,'Style','text','Tag','ocsp_text_maskChan',...
    'String','Mask channel','Units','pixels','HorizontalAlignment','left','FontSize',baseFont,...
    'Position',[xL y 110 18]);
handles.ocsp_edit_maskChan = uicontrol(parent,'Style','edit','Tag','ocsp_edit_maskChan',...
    'String','1','Units','pixels','FontSize',baseFont,...
    'Position',[xL+120 y 60 20]);

y = y - (rowH+lineGap);
handles.ocsp_check_useSDC = uicontrol(parent,'Style','checkbox','Tag','ocsp_check_useSDC',...
    'String','Warp masks using stage drift correction','Units','pixels','FontSize',baseFont,...
    'Position',[xL y colW 18]);

y = y - (rowH+lineGap);
handles.ocsp_check_preview = uicontrol(parent,'Style','checkbox','Tag','ocsp_check_preview',...
    'String','Save per-frame preview overlay','Units','pixels','FontSize',baseFont,...
    'Position',[xL y colW 18]);

y = y - (rowH+lineGap);
handles.ocsp_check_useLabeling = uicontrol(parent,'Style','checkbox','Tag','ocsp_check_useLabeling',...
    'String','Compute per-object stats (labeling)','Units','pixels','FontSize',baseFont,...
    'Position',[xL y colW 18], 'Callback', @(h,e) ocsp_useLabeling_Callback(h,e,guidata(h))); 

y = y - (rowH+lineGap);
handles.ocsp_text_minArea = uicontrol(parent,'Style','text','Tag','ocsp_text_minArea',...
    'String','Min area (px)','Units','pixels','HorizontalAlignment','left','FontSize',baseFont,...
    'Position',[xL y 110 18]);
handles.ocsp_edit_minArea = uicontrol(parent,'Style','edit','Tag','ocsp_edit_minArea',...
    'String','200','Units','pixels','FontSize',baseFont,...
    'Position',[xL+120 y 80 20]);

% --- Right column: Instance split + progress bar
yR = yTop - 42;
handles.ocsp_text_method = uicontrol(parent,'Style','text','Tag','ocsp_text_method',...
    'String','Method','Units','pixels','HorizontalAlignment','left','FontSize',baseFont,...
    'Position',[xR yR 60 18]);
handles.ocsp_popup_instanceMethod = uicontrol(parent,'Style','popupmenu','Tag','ocsp_popup_instanceMethod',...
    'String',{'cc (connected components)','watershed'},'Units','pixels','FontSize',baseFont,...
    'Position',[xR+70 yR 190 22], 'Callback', @(h,e) ocsp_instanceMethod_Callback(h,e,guidata(h)));
set(handles.ocsp_popup_instanceMethod,'UserData',{'cc','watershed'});

yR = yR - (rowH+lineGap);
handles.ocsp_text_sigma = uicontrol(parent,'Style','text','Tag','ocsp_text_sigma',...
    'String','Smooth sigma','Units','pixels','HorizontalAlignment','left','FontSize',baseFont,...
    'Position',[xR yR 90 18]);
handles.ocsp_edit_smoothSigma = uicontrol(parent,'Style','edit','Tag','ocsp_edit_smoothSigma',...
    'String','1.5','Units','pixels','FontSize',baseFont,...
    'Position',[xR+100 yR 60 20]);

yR = yR - (rowH+lineGap);
handles.ocsp_text_seedh = uicontrol(parent,'Style','text','Tag','ocsp_text_seedh',...
    'String','Seed h','Units','pixels','HorizontalAlignment','left','FontSize',baseFont,...
    'Position',[xR yR 90 18]);
handles.ocsp_edit_minSeedH = uicontrol(parent,'Style','edit','Tag','ocsp_edit_minSeedH',...
    'String','2.0','Units','pixels','FontSize',baseFont,...
    'Position',[xR+100 yR 60 20]);

yR = yR - (rowH+lineGap);
handles.ocsp_check_borderClear = uicontrol(parent,'Style','checkbox','Tag','ocsp_check_borderClear',...
    'String','Clear border objects','Units','pixels','FontSize',baseFont,...
    'Position',[xR yR colW 18]);

% Progress bar section
yR = yR - 32;
handles.ocsp_text_pb = uicontrol(parent,'Style','text','Tag','ocsp_text_pb',...
    'String','Progress bar','Units','pixels','HorizontalAlignment','left',...
    'FontWeight','bold','FontSize',baseFont+2,'Position',[xR yR colW 18]);

yR = yR - (rowH+lineGap);
handles.ocsp_check_useProgressBar = uicontrol(parent,'Style','checkbox','Tag','ocsp_check_useProgressBar',...
    'String','Show progress bar','Units','pixels','FontSize',baseFont,...
    'Position',[xR yR colW 18]);

yR = yR - (rowH+lineGap);
handles.ocsp_text_pbpos = uicontrol(parent,'Style','text','Tag','ocsp_text_pbpos',...
    'String','Position','Units','pixels','HorizontalAlignment','left','FontSize',baseFont,...
    'Position',[xR yR 70 18]);
handles.ocsp_popup_progressPos = uicontrol(parent,'Style','popupmenu','Tag','ocsp_popup_progressPos',...
    'String',{'0 Center','1 Upper right','2 Upper left','3 Lower left','4 Lower right','5 Random'},...
    'Units','pixels','FontSize',baseFont,'Position',[xR+80 yR 160 22]);
set(handles.ocsp_popup_progressPos,'UserData',[0 1 2 3 4 5]);

% Initial enable/disable for watershed controls
ocsp_instanceMethod_Callback(handles.ocsp_popup_instanceMethod, [], handles);
end


function handles = ocsp_initParamControls(handles, funParams, userData)
% Fill controls from funParams

% Mask process list
MD = [];
try
    MD = userData.MD;
catch
end
if isempty(MD)
    try
        MD = userData.crtProc.getOwner();
    catch
        MD = [];
    end
end
[maskString, maskData] = ocsp_buildMaskProcessList(MD);
set(handles.ocsp_popup_maskProcess, 'String', maskString, 'UserData', maskData);

mp = getFieldOr(funParams,'MaskProcessName','');
idx = find(strcmp(mp, maskData), 1, 'first');
if isempty(idx), idx = 1; end
set(handles.ocsp_popup_maskProcess,'Value',idx);

set(handles.ocsp_edit_maskChan,'String', num2str(getFieldOr(funParams,'MaskChannelIndex',1)));

safeSetValueDirect(handles.ocsp_check_useSDC, getFieldOr(funParams,'UseStageDriftCorrection',false));
safeSetValueDirect(handles.ocsp_check_preview, getFieldOr(funParams,'SavePerFrameTifPreview',false));

safeSetValueDirect(handles.ocsp_check_useLabeling, getFieldOr(funParams,'UseLabeling',true));
set(handles.ocsp_edit_minArea,'String', num2str(getFieldOr(funParams,'MinAreaPix',200)));

safeSetValueDirect(handles.ocsp_check_trackCells, getFieldOr(funParams,'TrackCells',true));
safeSetValueDirect(handles.ocsp_check_trackByOverlap, getFieldOr(funParams,'TrackByOverlap',true));
set(handles.ocsp_edit_minIoU,'String', num2str(getFieldOr(funParams,'MinIoU',0.10)));
set(handles.ocsp_edit_maxCentroidDist,'String', num2str(getFieldOr(funParams,'MaxCentroidDist',40)));

safeSetValueDirect(handles.ocsp_check_computeDFF0, getFieldOr(funParams,'ComputeDFF0',true));
set(handles.ocsp_edit_baselineFallbackN,'String', num2str(getFieldOr(funParams,'BaselineFallbackN',5)));

bf = getFieldOr(funParams,'BaselineFrames',1:10);
if isnumeric(bf)
    if isscalar(bf)
        bfStr = num2str(bf);
    else
        bfStr = mat2str(bf);
    end
else
    bfStr = char(bf);
end
set(handles.ocsp_edit_baselineFrames,'String', bfStr);

% Instance method
method = getFieldOr(funParams,'InstanceSegMethod','cc');
ud = get(handles.ocsp_popup_instanceMethod,'UserData');
vidx = find(strcmp(method, ud), 1, 'first');
if isempty(vidx), vidx = 1; end
set(handles.ocsp_popup_instanceMethod,'Value',vidx);

set(handles.ocsp_edit_smoothSigma,'String', num2str(getFieldOr(funParams,'SmoothSigma',1.5)));
set(handles.ocsp_edit_minSeedH,'String', num2str(getFieldOr(funParams,'MinSeedH',2.0)));
safeSetValueDirect(handles.ocsp_check_borderClear, getFieldOr(funParams,'BorderClear',false));

% Progress bar
safeSetValueDirect(handles.ocsp_check_useProgressBar, getFieldOr(funParams,'UseProgressBar',true));
ud2 = get(handles.ocsp_popup_progressPos,'UserData');
pp  = getFieldOr(funParams,'ProgressBarPosition',4);
vidx = find(ud2==pp, 1, 'first');
if isempty(vidx), vidx = find(ud2==4,1,'first'); end
set(handles.ocsp_popup_progressPos,'Value',vidx);

% Enable/disable based on labeling
ocsp_useLabeling_Callback(handles.ocsp_check_useLabeling, [], handles);
ocsp_instanceMethod_Callback(handles.ocsp_popup_instanceMethod, [], handles);
end


% =========================
% Callbacks for dynamic controls
% =========================

function ocsp_useLabeling_Callback(hObject, eventdata, handles)
useLab = logical(get(hObject,'Value'));
state = onOff(useLab);
% tracking related
set(handles.ocsp_check_trackCells,'Enable',state);
set(handles.ocsp_check_trackByOverlap,'Enable',state);
set(handles.ocsp_edit_minIoU,'Enable',state);
set(handles.ocsp_edit_maxCentroidDist,'Enable',state);
% min area still relevant when labeling
set(handles.ocsp_edit_minArea,'Enable',state);
end


function ocsp_instanceMethod_Callback(hObject, eventdata, handles)
ud = get(hObject,'UserData');
val = get(hObject,'Value');
method = ud{val};
ws = strcmpi(method,'watershed');
state = onOff(ws);
set(handles.ocsp_edit_smoothSigma,'Enable',state);
set(handles.ocsp_edit_minSeedH,'Enable',state);
set(handles.ocsp_check_borderClear,'Enable',state);
end


% =========================
% Helpers
% =========================

function [maskString, maskData] = ocsp_buildMaskProcessList(MD)
% Candidate processes that likely produce masks.
maskData = {'', 'MaskIntersectionProcess','MaskRefinementProcess','ThresholdProcess'};
maskString = {'(auto-detect)', 'MaskIntersectionProcess', 'MaskRefinementProcess', 'ThresholdProcess'};

try
    if ~isempty(MD) && isprop(MD,'processes_') && ~isempty(MD.processes_)
        extra = {};
        for i = 1:numel(MD.processes_)
            pr = MD.processes_{i};
            if isempty(pr), continue; end
            cn = class(pr);
            if contains(cn,'Mask','IgnoreCase',true) || contains(cn,'Threshold','IgnoreCase',true) || contains(cn,'Segmentation','IgnoreCase',true)
                extra{end+1} = cn; %#ok<AGROW>
            end
        end
        extra = unique(extra, 'stable');
        extra = setdiff(extra, maskData, 'stable');
        if ~isempty(extra)
            maskData = [maskData, extra];
            maskString = [maskString, extra];
        end
    end
catch
end
end


function v = getFieldOr(s, f, defaultVal)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f))
    v = s.(f);
else
    v = defaultVal;
end
end

function safeSetString(handles, tag, str)
if isfield(handles, tag)
    try, set(handles.(tag), 'String', str); end %#ok<TRYNC>
end
end

function safeSetValueDirect(h, val)
try
    set(h,'Value',double(val));
catch
end
end

function s = onOff(tf)
if tf, s = 'on'; else, s = 'off'; end
end
