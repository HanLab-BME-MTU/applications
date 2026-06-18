function varargout = filoGUI_helper(action, varargin)
%FILOGUI_HELPER  Shared utilities for FilopodiaForcePackage process GUIs.
%
% Usage (called from each process GUI):
%   ud = filoGUI_helper('init', mainFig, procID, paramDefs)
%       Parse caller args, create uifigure, populate param rows.
%       Returns userData struct stored in figure.
%   filoGUI_helper('apply', fig)
%       Read current values from fig, call proc.setPara, MD.save, refresh.
%   filoGUI_helper('cancel', fig)
%       Close without saving.
%
% paramDefs is a struct array with fields:
%   .name    : funParams field name
%   .label   : human-readable label
%   .type    : 'edit' (numeric) | 'editstr' (string) | 'checkbox' (logical)
%   .tooltip : tooltip string
% Sangyoon J. Han / 2026

switch action
    case 'init',  varargout{1} = initGUI(varargin{:});
    case 'apply', applyGUI(varargin{1});
    case 'cancel', delete(varargin{1});
end
end

% ===================================================================
function ud = initGUI(mainFig, procID, procClass, paramDefs)
% parse mainFig
ud_main = get(mainFig, 'UserData');
if isfield(ud_main,'MD'), MD = ud_main.MD(ud_main.id); else, MD = []; end
pkg = ud_main.crtPackage;
proc = pkg.getProcess(procID);
if isempty(proc)
    % create with defaults
    constr = pkg.getDefaultProcessConstructors(procID);
    proc = constr{1}(MD);
    MD.addProcess(proc);
end
funParams = proc.funParams_;

% figure
ROW_H = 32; PAD = 10; LBL_W = 260; EDT_W = 160; BTN_H = 36;
nP = numel(paramDefs);
FIG_H = PAD + nP*ROW_H + PAD + BTN_H + PAD;
FIG_W = PAD + LBL_W + PAD + EDT_W + PAD;
fig = uifigure('Name', [proc.getName() ' Settings'], ...
    'Position', [200 200 FIG_W FIG_H], 'Resize','off');

ud.mainFig   = mainFig;
ud.procID    = procID;
ud.MD        = MD;
ud.pkg       = pkg;
ud.proc      = proc;
ud.funParams = funParams;
ud.paramDefs = paramDefs;
ud.controls  = struct();

% param rows (bottom-up layout)
for k = 1:nP
    d = paramDefs(k);
    y = FIG_H - PAD - k*ROW_H + 4;
    uilabel(fig,'Text',[d.label ':'], ...
        'Position',[PAD, y, LBL_W, ROW_H-4], ...
        'HorizontalAlignment','right','FontSize',12,'Tooltip',d.tooltip);
    curVal = gf(funParams, d.name, []);
    switch d.type
        case 'checkbox'
            c = uicheckbox(fig,'Text','','Value',logical(curVal), ...
                'Position',[PAD+LBL_W+PAD, y+4, 30, ROW_H-8], ...
                'Tooltip',d.tooltip);
        case 'editstr'
            c = uieditfield(fig,'text','Value',num2str(curVal), ...
                'Position',[PAD+LBL_W+PAD, y, EDT_W, ROW_H-4], ...
                'Tooltip',d.tooltip);
        otherwise  % 'edit' numeric
            if ischar(curVal) || isstring(curVal)
                str = char(curVal);
            elseif numel(curVal)>1
                str = mat2str(curVal);
            else
                str = num2str(curVal);
            end
            c = uieditfield(fig,'text','Value',str, ...
                'Position',[PAD+LBL_W+PAD, y, EDT_W, ROW_H-4], ...
                'Tooltip',d.tooltip);
    end
    ud.controls.(d.name) = c;
end

% Done / Cancel buttons
y = PAD;
uibutton(fig,'Text','Done','Position',[FIG_W-2*(BTN_H*2+PAD), y, BTN_H*2, BTN_H], ...
    'ButtonPushedFcn', @(~,~) filoGUI_helper('apply', fig));
uibutton(fig,'Text','Cancel','Position',[FIG_W-(BTN_H*2+PAD), y, BTN_H*2, BTN_H], ...
    'ButtonPushedFcn', @(~,~) filoGUI_helper('cancel', fig));

fig.UserData = ud;
end

% ===================================================================
function applyGUI(fig)
ud = fig.UserData;
fp = ud.funParams;
pd = ud.paramDefs;
for k = 1:numel(pd)
    d = pd(k); c = ud.controls.(d.name);
    switch d.type
        case 'checkbox', fp.(d.name) = c.Value;
        case 'editstr',  fp.(d.name) = c.Value;
        otherwise
            str = strtrim(c.Value);
            v = str2num(str); %#ok<ST2NM>
            if isempty(v), v = str; end
            fp.(d.name) = v;
    end
end
try
    ud.proc.setPara(fp);
    ud.MD.save();
    % refresh packageGUI
    handles_main = guidata(ud.mainFig);
    packageGUI('refreshPackage_Callback', handles_main.figure1, [], handles_main);
catch ME
    uialert(fig, ME.message, 'Error saving params');
    return;
end
delete(fig);
end

function v = gf(s,n,d), if isfield(s,n)&&~isempty(s.(n)), v=s.(n); else, v=d; end, end
