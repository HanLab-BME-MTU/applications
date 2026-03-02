function varargout = otherChannelSamplingProcessGUI(varargin)
% otherChannelSamplingProcessGUI M-file for otherChannelSamplingProcessGUI.fig
%
% This GUI is intentionally based on forceFieldCalculationProcessGUI's FIG
% layout, but the controls are repurposed for OtherChannelSamplingProcess.
%
% Mapping (FIG tag -> funParams field):
%   edit_basisClassTblPath      -> OutputDirectory
%   popupmenu_method            -> ChannelIndex
%   popupmenu_solMethodBEM      -> MaskProcessName
%   edit_meshPtsFwdSol          -> MaskChannelIndex
%   checkbox_lastToFirst        -> UseStageDriftCorrection
%   checkbox_everyframe         -> SavePerFrameTifPreview
%   useLcurve                   -> UseLabeling
%   edit_PoissonRatio           -> MinAreaPix
%   setROIfromForcemap          -> ComputeDFF0
%   edit_LcurveFactor           -> BaselineFrames (string, parsed by str2num)
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

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:});

% Ensure handles contains all GUI component fields (GUIDE sometimes passes empty handles).
handles = guihandles(hObject);


processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:});

% Refresh handles after processGUI_OpeningFcn may have updated guidata
handles = guidata(hObject);


% Repurpose labels / hide unused widgets (we keep the FIG but retitle controls)
set(handles.figure1, 'Name', 'Setting - Other Channel Sampling');

% Update static text labels
safeSetString(handles, 'text_processName', 'OtherChannelSamplingProcess settings');
safeSetString(handles, 'text_method', 'Sampling channel');
safeSetString(handles, 'text49', 'Mask process');
safeSetString(handles, 'text_meshPtsFwdSol', 'Mask channel index');
safeSetString(handles, 'text_basisClassTblPath', 'Output directory');
safeSetString(handles, 'pushbutton_basisClassTblPath', 'Browse');
safeSetString(handles, 'text_PoissonRatio', 'Min object area (pixels)');
safeSetString(handles, 'text_LcurveFactor', 'Baseline frames (e.g., 1:10)');

safeSetString(handles, 'checkbox_lastToFirst', 'Use stage drift correction (warp masks)');
safeSetString(handles, 'checkbox_everyframe', 'Save per-frame preview');
safeSetString(handles, 'useLcurve', 'Compute per-object stats (labeling)');
safeSetString(handles, 'setROIfromForcemap', 'Compute dF/F0');

% Hide unrelated controls/panels from the inherited FIG
safeSetVisible(handles, 'uipanel_BEM', 'off');
safeSetVisible(handles, 'groupCornerOptimal', 'off');
safeSetVisible(handles, 'text_YoungModulus', 'off');
safeSetVisible(handles, 'edit_YoungModulus', 'off');
safeSetVisible(handles, 'text_thickness', 'off');
safeSetVisible(handles, 'edit_thickness', 'off');
safeSetVisible(handles, 'text_regParam', 'off');
safeSetVisible(handles, 'edit_regParam', 'off');
safeSetVisible(handles, 'text_tolx', 'off');
safeSetVisible(handles, 'edit_tolx', 'off');
safeSetVisible(handles, 'text52', 'off');
safeSetVisible(handles, 'text54', 'off');
safeSetVisible(handles, 'lcorner', 'off');
safeSetVisible(handles, 'optimal', 'off');

% Load current funParams
userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end
funParams = userData.crtProc.funParams_;

% Build channel popup
try
    MD = userData.MD;
catch
    MD = [];
end
if isempty(MD)
    % processGUI_OpeningFcn usually sets userData.MD; but keep safe fallback
    try
        MD = userData.crtProc.getOwner();
    catch
        MD = [];
    end
end

% [chString, chData] = buildChannelList(MD);
% set(handles.popupmenu_method, 'String', chString, 'UserData', chData);
% 
% if isfield(funParams,'ChannelIndex') && ~isempty(funParams.ChannelIndex) && funParams.ChannelIndex >= 1 && funParams.ChannelIndex <= numel(chData)
%     set(handles.popupmenu_method, 'Value', funParams.ChannelIndex);
% else
%     set(handles.popupmenu_method, 'Value', 1);
% end


% Checkboxes
safeSetValue(handles, 'checkbox_lastToFirst', getFieldOr(funParams,'UseStageDriftCorrection',false));
safeSetValue(handles, 'checkbox_everyframe', getFieldOr(funParams,'SavePerFrameTifPreview',false));
safeSetValue(handles, 'useLcurve', getFieldOr(funParams,'UseLabeling',false));
safeSetValue(handles, 'setROIfromForcemap', getFieldOr(funParams,'ComputeDFF0',false));

% % Min area
% if isfield(funParams,'MinAreaPix')
%     set(handles.edit_PoissonRatio, 'String', num2str(funParams.MinAreaPix));
% else
%     set(handles.edit_PoissonRatio, 'String', '50');
% end

% Baseline frames
if isfield(funParams,'BaselineFrames')
    bf = funParams.BaselineFrames;
    if isnumeric(bf)
        if isscalar(bf)
            bfStr = num2str(bf);
        else
            bfStr = mat2str(bf);
        end
    else
        bfStr = bf;
    end
    % set(handles.edit_LcurveFactor, 'String', bfStr);
else
    % set(handles.edit_LcurveFactor, 'String', '1');
end


% Create/initialize additional parameter controls (dynamic; no FIG edit needed)
handles = ensureOtherChannelParamControls(handles, funParams);
guidata(hObject, handles);
% (handles updated)

% Output directory
% if isfield(funParams, 'OutputDirectory')
%     set(handles.output, 'String', funParams.OutputDirectory);
% end

guidata(hObject, handles);
end

% --- Outputs from this function are returned to the command line.
function varargout = otherChannelSamplingProcessGUI_OutputFcn(~, ~, handles)
    varargout{1} = handles.figure1;

end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(~, eventdata, handles)
if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end
end


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end

funParams = userData.crtProc.funParams_;

% % Output directory
% outDir = get(handles.edit_basisClassTblPath, 'String');
% if isempty(outDir)
%     errordlg('Please enter a valid output directory.','Setting Error','modal');
%     return;
% end
% funParams.OutputDirectory = outDir;

% ChannelIndex
% props = get(handles.popupmenu_method, {'UserData','Value'});
% chData = props{1};
% chVal  = props{2};
% if isempty(chData) || chVal < 1 || chVal > numel(chData)
%     funParams.ChannelIndex = 1;
% else
%     funParams.ChannelIndex = chVal; % index into channel list
% end

% Mask process name
% props = get(handles.popupmenu_solMethodBEM, {'UserData','Value'});
% maskData = props{1};
% maskVal  = props{2};
% if isempty(maskData)
%     funParams.MaskProcessName = '';
% else
%     funParams.MaskProcessName = maskData{maskVal};
% end

% Mask channel index
% maskChanStr = get(handles.edit_meshPtsFwdSol, 'String');
% maskChan = str2double(maskChanStr);
% if isnan(maskChan) || maskChan < 1
%     errordlg('Mask channel index must be a positive integer.','Setting Error','modal');
%     return;
% end
% funParams.MaskChannelIndex = round(maskChan);

% Options
funParams.UseStageDriftCorrection = logical(get(handles.checkbox_lastToFirst, 'Value'));
funParams.SavePerFrameTifPreview  = logical(get(handles.checkbox_everyframe, 'Value'));
funParams.UseLabeling            = logical(get(handles.useLcurve, 'Value'));
funParams.ComputeDFF0            = logical(get(handles.setROIfromForcemap, 'Value'));

% MinAreaPix
minAreaStr = get(handles.edit_PoissonRatio, 'String');
minArea = str2double(minAreaStr);
if isnan(minArea) || minArea < 0
    errordlg('Min object area (pixels) must be a nonnegative number.','Setting Error','modal');
    return;
end
funParams.MinAreaPix = round(minArea);

% BaselineFrames
bfStr = get(handles.edit_LcurveFactor, 'String');
bf = str2num(bfStr); %#ok<ST2NM> allow MATLAB-style vectors like 1:10
if isempty(bf) || ~isnumeric(bf)
    errordlg('Baseline frames must be a numeric vector (e.g., 1:10 or [1 2 3]).','Setting Error','modal');
    return;
end
funParams.BaselineFrames = bf(:)';

% -------- Additional controls (tracking / instance split / progressbar / dF/F0 per track) --------
if isfield(handles,'checkbox_trackCells')
    funParams.TrackCells = logical(get(handles.checkbox_trackCells,'Value'));
end
if isfield(handles,'checkbox_trackByOverlap')
    funParams.TrackByOverlap = logical(get(handles.checkbox_trackByOverlap,'Value'));
end
if isfield(handles,'edit_minIoU')
    v = str2double(get(handles.edit_minIoU,'String'));
    if isnan(v) || v < 0 || v > 1
        errordlg('Min IoU must be between 0 and 1.','Setting Error','modal'); return;
    end
    funParams.MinIoU = v;
end
if isfield(handles,'edit_maxCentroidDist')
    v = str2double(get(handles.edit_maxCentroidDist,'String'));
    if isnan(v) || v <= 0
        errordlg('Max centroid distance must be a positive number.','Setting Error','modal'); return;
    end
    funParams.MaxCentroidDist = v;
end
if isfield(handles,'edit_baselineFallbackN')
    v = str2double(get(handles.edit_baselineFallbackN,'String'));
    if isnan(v) || v < 0
        errordlg('BaselineFallbackN must be >= 0.','Setting Error','modal'); return;
    end
    funParams.BaselineFallbackN = round(v);
end
if isfield(handles,'popupmenu_instanceSegMethod')
    ud = get(handles.popupmenu_instanceSegMethod,'UserData');
    val = get(handles.popupmenu_instanceSegMethod,'Value');
    if iscell(ud) && val>=1 && val<=numel(ud)
        funParams.InstanceSegMethod = ud{val};
    end
end
if isfield(handles,'edit_smoothSigma')
    v = str2double(get(handles.edit_smoothSigma,'String'));
    if isnan(v) || v < 0
        errordlg('SmoothSigma must be >= 0.','Setting Error','modal'); return;
    end
    funParams.SmoothSigma = v;
end
if isfield(handles,'edit_minSeedH')
    v = str2double(get(handles.edit_minSeedH,'String'));
    if isnan(v) || v < 0
        errordlg('MinSeedH must be >= 0.','Setting Error','modal'); return;
    end
    funParams.MinSeedH = v;
end
if isfield(handles,'checkbox_borderClear')
    funParams.BorderClear = logical(get(handles.checkbox_borderClear,'Value'));
end
if isfield(handles,'checkbox_useProgressBar')
    funParams.UseProgressBar = logical(get(handles.checkbox_useProgressBar,'Value'));
end
if isfield(handles,'popupmenu_progressBarPosition')
    ud = get(handles.popupmenu_progressBarPosition,'UserData');
    val = get(handles.popupmenu_progressBarPosition,'Value');
    if isnumeric(ud) && val>=1 && val<=numel(ud)
        funParams.ProgressBarPosition = ud(val);
    end
end

% Sanity check (data only)
try
    userData.crtProc.sanityCheck;
catch ME
    errordlg([ME.message ' Please double check your data.'], 'Setting Error', 'modal');
    return;
end

% Apply
processGUI_ApplyFcn(hObject, eventdata, handles, funParams);
end

% --- Executes on button press in pushbutton_basisClassTblPath.
function pushbutton_basisClassTblPath_Callback(hObject, eventdata, handles)
% Browse output directory
startDir = get(handles.edit_basisClassTblPath,'String');
if isempty(startDir) || ~exist(startDir,'dir')
    startDir = pwd;
end
selDir = uigetdir(startDir, 'Select output directory');
if isequal(selDir,0), return; end
set(handles.edit_basisClassTblPath, 'String', selDir);
end

% ---------------- Helper functions ----------------
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



function safeSetVisible(handles, tag, vis)
if isfield(handles, tag)
    try, set(handles.(tag), 'Visible', vis); end %#ok<TRYNC>
end
end

function safeSetValue(handles, tag, val)
if isfield(handles, tag)
    try, set(handles.(tag), 'Value', double(val)); end %#ok<TRYNC>
end
end

function [chString, chData] = buildChannelList(MD)
chString = {'Channel 1'};
chData = {1};
try
    if ~isempty(MD) && isprop(MD,'channels_') && ~isempty(MD.channels_)
        nCh = numel(MD.channels_);
        chString = cell(nCh,1);
        chData = cell(nCh,1);
        for i = 1:nCh
            nm = '';
            try, nm = MD.channels_(i).name_; end %#ok<TRYNC>
            if isempty(nm)
                chString{i} = sprintf('Channel %d', i);
            else
                chString{i} = sprintf('Channel %d: %s', i, nm);
            end
            chData{i} = i;
        end
    end
catch
end
end

function [maskString, maskData] = buildMaskProcessList(MD)
% Return list of candidate processes that likely produce masks.
maskData = {'', 'MaskIntersectionProcess','MaskRefinementProcess','ThresholdProcess'};
maskString = {'(auto-detect)', 'MaskIntersectionProcess', 'MaskRefinementProcess', 'ThresholdProcess'};

try
    if ~isempty(MD) && isprop(MD,'processes_') && ~isempty(MD.processes_)
        extra = {};
        for i = 1:numel(MD.processes_)
            p = MD.processes_{i};
            if isempty(p), continue; end
            cn = class(p);
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
end % buildMaskProcessList

% =========================
% Dynamic parameter controls
% =========================
function handles = ensureOtherChannelParamControls(handles, funParams)
% Create parameter controls inside the Parameters panel at runtime (so the .fig doesn't need editing).
% If controls already exist (because you added them in GUIDE), this function will only initialize values.

% Find parent panel for parameters
if isfield(handles,'uipanel_1')
    parent = handles.uipanel_1;
else
    parent = handles.figure1;
end

% Use pixel units for predictable layout
try, set(parent,'Units','pixels'); end
posP = get(parent,'Position'); %#ok<NASGU>

% Helper to create if missing
    function h = getOrCreate(tag, createFcn)
        if isfield(handles, tag) && ishghandle(handles.(tag))
            h = handles.(tag);
        else
            h = createFcn();
            handles.(tag) = h;
        end
    end

y0 = 135; dy = 22;
xL = 15; wL = 160;
xC = 185; wC = 80;
xR = 285; wR = 205;

% --- Tracking header ---
getOrCreate('text_tracking_header', @() uicontrol('Parent',parent,'Style','text','Tag','text_tracking_header', ...
    'String','Tracking','Units','pixels','HorizontalAlignment','left','FontWeight','bold', ...
    'Position',[xL y0 wL 18]));

% TrackCells checkbox
getOrCreate('checkbox_trackCells', @() uicontrol('Parent',parent,'Style','checkbox','Tag','checkbox_trackCells', ...
    'String','Track cells across frames','Units','pixels','Position',[xL y0-dy wR 18]));

% TrackByOverlap checkbox
getOrCreate('checkbox_trackByOverlap', @() uicontrol('Parent',parent,'Style','checkbox','Tag','checkbox_trackByOverlap', ...
    'String','Link by overlap (IoU)','Units','pixels','Position',[xL y0-2*dy wR 18]));

% MinIoU label + edit
getOrCreate('text_minIoU', @() uicontrol('Parent',parent,'Style','text','Tag','text_minIoU', ...
    'String','Min IoU','Units','pixels','HorizontalAlignment','left','Position',[xL y0-3*dy wL 18]));
getOrCreate('edit_minIoU', @() uicontrol('Parent',parent,'Style','edit','Tag','edit_minIoU', ...
    'Units','pixels','Position',[xC y0-3*dy wC 20]));

% MaxCentroidDist label + edit
getOrCreate('text_maxCentroidDist', @() uicontrol('Parent',parent,'Style','text','Tag','text_maxCentroidDist', ...
    'String','Max centroid dist (px)','Units','pixels','HorizontalAlignment','left','Position',[xL y0-4*dy wL 18]));
getOrCreate('edit_maxCentroidDist', @() uicontrol('Parent',parent,'Style','edit','Tag','edit_maxCentroidDist', ...
    'Units','pixels','Position',[xC y0-4*dy wC 20]));

% --- dF/F0 header ---
y1 = y0-5*dy-5;
getOrCreate('text_dff_header', @() uicontrol('Parent',parent,'Style','text','Tag','text_dff_header', ...
    'String','dF/F0','Units','pixels','HorizontalAlignment','left','FontWeight','bold', ...
    'Position',[xL y1 wL 18]));
getOrCreate('text_baselineFallbackN', @() uicontrol('Parent',parent,'Style','text','Tag','text_baselineFallbackN', ...
    'String','Baseline fallback N','Units','pixels','HorizontalAlignment','left','Position',[xL y1-dy wL 18]));
getOrCreate('edit_baselineFallbackN', @() uicontrol('Parent',parent,'Style','edit','Tag','edit_baselineFallbackN', ...
    'Units','pixels','Position',[xC y1-dy wC 20]));

% --- Instance split header ---
y2 = y1-2*dy-8;
getOrCreate('text_instance_header', @() uicontrol('Parent',parent,'Style','text','Tag','text_instance_header', ...
    'String','Instance split (from mask)','Units','pixels','HorizontalAlignment','left','FontWeight','bold', ...
    'Position',[xR y0 wR 18]));
getOrCreate('text_instanceMethod', @() uicontrol('Parent',parent,'Style','text','Tag','text_instanceMethod', ...
    'String','Method','Units','pixels','HorizontalAlignment','left','Position',[xR y0-dy wL 18]));
getOrCreate('popupmenu_instanceSegMethod', @() uicontrol('Parent',parent,'Style','popupmenu','Tag','popupmenu_instanceSegMethod', ...
    'Units','pixels','Position',[xR+70 y0-dy 140 20]));
getOrCreate('text_smoothSigma', @() uicontrol('Parent',parent,'Style','text','Tag','text_smoothSigma', ...
    'String','Smooth sigma','Units','pixels','HorizontalAlignment','left','Position',[xR y0-2*dy wL 18]));
getOrCreate('edit_smoothSigma', @() uicontrol('Parent',parent,'Style','edit','Tag','edit_smoothSigma', ...
    'Units','pixels','Position',[xR+70 y0-2*dy 60 20]));
getOrCreate('text_minSeedH', @() uicontrol('Parent',parent,'Style','text','Tag','text_minSeedH', ...
    'String','Seed h','Units','pixels','HorizontalAlignment','left','Position',[xR+140 y0-2*dy 60 18]));
getOrCreate('edit_minSeedH', @() uicontrol('Parent',parent,'Style','edit','Tag','edit_minSeedH', ...
    'Units','pixels','Position',[xR+190 y0-2*dy 60 20]));
getOrCreate('checkbox_borderClear', @() uicontrol('Parent',parent,'Style','checkbox','Tag','checkbox_borderClear', ...
    'String','Clear border objects','Units','pixels','Position',[xR y0-3*dy wR 18]));

% --- Progress bar controls ---
getOrCreate('text_pb_header', @() uicontrol('Parent',parent,'Style','text','Tag','text_pb_header', ...
    'String','Progress','Units','pixels','HorizontalAlignment','left','FontWeight','bold', ...
    'Position',[xR y2 wR 18]));
getOrCreate('checkbox_useProgressBar', @() uicontrol('Parent',parent,'Style','checkbox','Tag','checkbox_useProgressBar', ...
    'String','Show progress bar','Units','pixels','Position',[xR y2-dy wR 18]));
getOrCreate('text_pb_pos', @() uicontrol('Parent',parent,'Style','text','Tag','text_pb_pos', ...
    'String','Position','Units','pixels','HorizontalAlignment','left','Position',[xR y2-2*dy 60 18]));
getOrCreate('popupmenu_progressBarPosition', @() uicontrol('Parent',parent,'Style','popupmenu','Tag','popupmenu_progressBarPosition', ...
    'Units','pixels','Position',[xR+70 y2-2*dy 180 20]));

% Initialize values from funParams
safeSetValue(handles,'checkbox_trackCells', getFieldOr(funParams,'TrackCells',true));
safeSetValue(handles,'checkbox_trackByOverlap', getFieldOr(funParams,'TrackByOverlap',true));
safeSetString(handles,'edit_minIoU', num2str(getFieldOr(funParams,'MinIoU',0.10)));
safeSetString(handles,'edit_maxCentroidDist', num2str(getFieldOr(funParams,'MaxCentroidDist',40)));
safeSetString(handles,'edit_baselineFallbackN', num2str(getFieldOr(funParams,'BaselineFallbackN',5)));

% Instance method popup
set(handles.popupmenu_instanceSegMethod,'String',{'cc (connected components)','watershed'},'UserData',{'cc','watershed'});
m = getFieldOr(funParams,'InstanceSegMethod','cc');
ud = get(handles.popupmenu_instanceSegMethod,'UserData');
v = find(strcmpi(m, ud), 1, 'first'); if isempty(v), v = 1; end
set(handles.popupmenu_instanceSegMethod,'Value',v);
safeSetString(handles,'edit_smoothSigma', num2str(getFieldOr(funParams,'SmoothSigma',1.5)));
safeSetString(handles,'edit_minSeedH', num2str(getFieldOr(funParams,'MinSeedH',2.0)));
safeSetValue(handles,'checkbox_borderClear', getFieldOr(funParams,'BorderClear',false));

% progress bar popup
set(handles.popupmenu_progressBarPosition,'String',{'0 Center','1 Upper right','2 Upper left','3 Lower left','4 Lower right','5 Random'},'UserData',[0 1 2 3 4 5]);
pos = getFieldOr(funParams,'ProgressBarPosition',4);
ud2 = get(handles.popupmenu_progressBarPosition,'UserData');
vidx = find(ud2==pos,1,'first'); if isempty(vidx), vidx = find(ud2==4,1,'first'); end
set(handles.popupmenu_progressBarPosition,'Value',vidx);
safeSetValue(handles,'checkbox_useProgressBar', getFieldOr(funParams,'UseProgressBar',true));

% Wire callbacks for enable/disable behavior
set(handles.popupmenu_instanceSegMethod,'Callback',@(h,e) ocsp_instanceSegMethod_Callback(h,e,guidata(h)));

% Apply initial enabling
ocsp_instanceSegMethod_Callback(handles.popupmenu_instanceSegMethod, [], handles);
end

function ocsp_useLabeling_Callback(hObject, ~, handles)
tf = logical(get(hObject,'Value'));
safeSetEnable(handles,'checkbox_trackCells', tf);
safeSetEnable(handles,'checkbox_trackByOverlap', tf);
safeSetEnable(handles,'edit_minIoU', tf);
safeSetEnable(handles,'edit_maxCentroidDist', tf);
safeSetEnable(handles,'edit_baselineFallbackN', tf);
end

function ocsp_instanceSegMethod_Callback(hObject, ~, handles)
ud = get(hObject,'UserData');
v = get(hObject,'Value');
method = 'cc';
if iscell(ud) && v>=1 && v<=numel(ud), method = ud{v}; end
isWS = strcmpi(method,'watershed');
safeSetEnable(handles,'edit_smoothSigma', isWS);
safeSetEnable(handles,'edit_minSeedH', isWS);
safeSetEnable(handles,'checkbox_borderClear', isWS);
end

function safeSetEnable(handles, tag, state)
try
    h = handles.(tag);
    if state, set(h,'Enable','on'); else, set(h,'Enable','off'); end
catch
end
end




