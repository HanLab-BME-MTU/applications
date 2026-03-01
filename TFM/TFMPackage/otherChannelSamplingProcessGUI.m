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


% --- Executes just before otherChannelSamplingProcessGUI is made visible.
function otherChannelSamplingProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:});

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

% Output directory
if isfield(funParams, 'OutputDirectory')
    set(handles.edit_basisClassTblPath, 'String', funParams.OutputDirectory);
end

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

[chString, chData] = buildChannelList(MD);
set(handles.popupmenu_method, 'String', chString, 'UserData', chData);

if isfield(funParams,'ChannelIndex') && ~isempty(funParams.ChannelIndex) && funParams.ChannelIndex >= 1 && funParams.ChannelIndex <= numel(chData)
    set(handles.popupmenu_method, 'Value', funParams.ChannelIndex);
else
    set(handles.popupmenu_method, 'Value', 1);
end

% Build mask-process popup
[maskString, maskData] = buildMaskProcessList(MD);
set(handles.popupmenu_solMethodBEM, 'String', maskString, 'UserData', maskData);

% Set mask selection
if isfield(funParams,'MaskProcessName') && ~isempty(funParams.MaskProcessName)
    v = find(strcmp(funParams.MaskProcessName, maskData), 1, 'first');
    if ~isempty(v), set(handles.popupmenu_solMethodBEM, 'Value', v); end
else
    set(handles.popupmenu_solMethodBEM, 'Value', 1);
end

% Mask channel index
if isfield(funParams,'MaskChannelIndex')
    set(handles.edit_meshPtsFwdSol, 'String', num2str(funParams.MaskChannelIndex));
else
    set(handles.edit_meshPtsFwdSol, 'String', '1');
end

% Checkboxes
safeSetValue(handles, 'checkbox_lastToFirst', getFieldOr(funParams,'UseStageDriftCorrection',false));
safeSetValue(handles, 'checkbox_everyframe', getFieldOr(funParams,'SavePerFrameTifPreview',false));
safeSetValue(handles, 'useLcurve', getFieldOr(funParams,'UseLabeling',false));
safeSetValue(handles, 'setROIfromForcemap', getFieldOr(funParams,'ComputeDFF0',false));

% Min area
if isfield(funParams,'MinAreaPix')
    set(handles.edit_PoissonRatio, 'String', num2str(funParams.MinAreaPix));
else
    set(handles.edit_PoissonRatio, 'String', '50');
end

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
    set(handles.edit_LcurveFactor, 'String', bfStr);
else
    set(handles.edit_LcurveFactor, 'String', '1');
end

guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = otherChannelSamplingProcessGUI_OutputFcn(~, ~, handles)
varargout{1} = handles.output;


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(~, eventdata, handles)
if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end

funParams = userData.crtProc.funParams_;

% Output directory
outDir = get(handles.edit_basisClassTblPath, 'String');
if isempty(outDir)
    errordlg('Please enter a valid output directory.','Setting Error','modal');
    return;
end
funParams.OutputDirectory = outDir;

% ChannelIndex
props = get(handles.popupmenu_method, {'UserData','Value'});
chData = props{1};
chVal  = props{2};
if isempty(chData) || chVal < 1 || chVal > numel(chData)
    funParams.ChannelIndex = 1;
else
    funParams.ChannelIndex = chVal; % index into channel list
end

% Mask process name
props = get(handles.popupmenu_solMethodBEM, {'UserData','Value'});
maskData = props{1};
maskVal  = props{2};
if isempty(maskData)
    funParams.MaskProcessName = '';
else
    funParams.MaskProcessName = maskData{maskVal};
end

% Mask channel index
maskChanStr = get(handles.edit_meshPtsFwdSol, 'String');
maskChan = str2double(maskChanStr);
if isnan(maskChan) || maskChan < 1
    errordlg('Mask channel index must be a positive integer.','Setting Error','modal');
    return;
end
funParams.MaskChannelIndex = round(maskChan);

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

% Sanity check (data only)
try
    userData.crtProc.sanityCheck;
catch ME
    errordlg([ME.message ' Please double check your data.'], 'Setting Error', 'modal');
    return;
end

% Apply
processGUI_ApplyFcn(hObject, eventdata, handles, funParams);


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


% ---------------- Helper functions ----------------
function v = getFieldOr(s, f, defaultVal)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f))
    v = s.(f);
else
    v = defaultVal;
end

function safeSetString(handles, tag, str)
if isfield(handles, tag)
    try, set(handles.(tag), 'String', str); end %#ok<TRYNC>
end

function safeSetVisible(handles, tag, vis)
if isfield(handles, tag)
    try, set(handles.(tag), 'Visible', vis); end %#ok<TRYNC>
end

function safeSetValue(handles, tag, val)
if isfield(handles, tag)
    try, set(handles.(tag), 'Value', double(val)); end %#ok<TRYNC>
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
