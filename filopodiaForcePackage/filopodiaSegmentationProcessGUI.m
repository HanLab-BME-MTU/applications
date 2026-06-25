function varargout = filopodiaSegmentationProcessGUI(varargin)
%FILOPODIASEGMENTATIONPROCESSGUI  Settings window for Process 1 (body
%segmentation) of the FilopodiaForcePackage. Provides channel selection and
%an interactive preview: pick a frame, adjust the body threshold (auto
%method + scale, or a manual intensity level via the slider) and morphology,
%and see the resulting cell-body outline overlaid on the talin image live.
%
% Called by packageGUI:  filopodiaSegmentationProcessGUI('mainFig',fig,procID)
% Settings are written to the process for THIS MovieData only.
% Sangyoon J. Han / 2026

if numel(varargin) < 3
    error('Call from packageGUI only: crtProcGUI(''mainFig'',mainFig,procID)');
end
mainFig = varargin{2};
procID  = varargin{3};

ud_main = get(mainFig,'UserData');
pkg = ud_main.crtPackage;
MD  = pkg.getOwner();
proc = pkg.getProcess(procID);
if isempty(proc)
    constrs = pkg.getDefaultProcessConstructors(procID);
    proc = constrs{1}(MD, MD.outputDirectory_);
    MD.addProcess(proc);
end
fp = proc.funParams_;

W = 980; H = 620;
fig = uifigure('Name',[proc.getName() ' Settings'],'Position',[120 120 W H],'Resize','on');
g = uigridlayout(fig,[1 2]);
g.ColumnWidth = {360,'1x'}; g.RowHeight = {'1x'};

left = uipanel(g,'Title','Settings');
sp = uigridlayout(left,[23 2]);
sp.ColumnWidth = {150,'1x'}; sp.RowHeight = repmat({26},1,23);
sp.Padding=[10 10 10 10]; sp.RowSpacing=4;

ud = struct();
ud.mainFig=mainFig; ud.procID=procID; ud.MD=MD; ud.pkg=pkg; ud.proc=proc; ud.fp=fp;
ud.previewFig = -1;

r = 0;
r=r+1; lbl(sp,r,'Channel');
chanNames = MD.getChannelPaths();
shortNames = cellfun(@(s) shortenPath(s), chanNames, 'unif',false);
ci = fp.ChannelIndex; if isempty(ci)||ci<1, ci=1; end
ud.ddChannel = uidropdown(sp,'Items',shortNames,'ItemsData',1:numel(chanNames),'Value',min(ci,numel(chanNames)));
ud.ddChannel.Layout.Row=r; ud.ddChannel.Layout.Column=2;

r=r+1; lbl(sp,r,'Threshold method');
ud.ddMethod = uidropdown(sp,'Items',{'rosin','otsu','manual'},'Value',methodFromParams(fp));
ud.ddMethod.Layout.Row=r; ud.ddMethod.Layout.Column=2;

r=r+1; lbl(sp,r,'Threshold scale (x)');
ud.efScale = uieditfield(sp,'numeric','Value',gf(fp,'ThreshScale',1),'Limits',[0 Inf],'Tooltip','Multiplier on auto level; <1 recovers dim talin band');
ud.efScale.Layout.Row=r; ud.efScale.Layout.Column=2;

r=r+1; lbl(sp,r,'Manual level');
ud.slLevel = uislider(sp,'Limits',[0 1],'Value',0);
ud.slLevel.Layout.Row=r; ud.slLevel.Layout.Column=2;
r=r+1; lbl(sp,r,'  fine adjust');
fineP = uipanel(sp,'BorderType','none'); fineP.Layout.Row=r; fineP.Layout.Column=2;
fg = uigridlayout(fineP,[1 3]); fg.Padding=[0 0 0 0]; fg.ColumnWidth={36,'1x',36}; fg.ColumnSpacing=4;
ud.btnMinus = uibutton(fg,'Text',char(9664),'ButtonPushedFcn',@(~,~)nudgeLevel(fig,-1));
ud.btnMinus.Layout.Row=1; ud.btnMinus.Layout.Column=1;
ud.efLevel = uieditfield(fg,'numeric','Value',0,'Tooltip','Used only when method = manual');
ud.efLevel.Layout.Row=1; ud.efLevel.Layout.Column=2;
ud.btnPlus = uibutton(fg,'Text',char(9654),'ButtonPushedFcn',@(~,~)nudgeLevel(fig,+1));
ud.btnPlus.Layout.Row=1; ud.btnPlus.Layout.Column=3;

r=r+1; lbl(sp,r,'Display brightness');
ud.slBright = uislider(sp,'Limits',[0.05 1],'Value',1, ...
    'Tooltip','Display only: lower = brighter dim features (does not affect segmentation)');
ud.slBright.Layout.Row=r; ud.slBright.Layout.Column=2;

r=r+1; lbl(sp,r,'Body blur sigma (px)');
ud.efBlur = uieditfield(sp,'numeric','Value',gf(fp,'GaussianBlurSigma',2),'Limits',[0 Inf]);
ud.efBlur.Layout.Row=r; ud.efBlur.Layout.Column=2;

r=r+1; lbl(sp,r,'Body open radius (px)');
ud.efOpen = uieditfield(sp,'numeric','Value',gf(fp,'BodyOpenRadius',8),'Limits',[0 Inf]);
ud.efOpen.Layout.Row=r; ud.efOpen.Layout.Column=2;

r=r+1; lbl(sp,r,'Body close radius (px)');
ud.efClose = uieditfield(sp,'numeric','Value',gf(fp,'BodyClosingRadius',8),'Limits',[0 Inf]);
ud.efClose.Layout.Row=r; ud.efClose.Layout.Column=2;

r=r+1; lbl(sp,r,'Body min area (px)');
ud.efMinArea = uieditfield(sp,'numeric','Value',gf(fp,'BodyMinArea',500),'Limits',[0 Inf]);
ud.efMinArea.Layout.Row=r; ud.efMinArea.Layout.Column=2;

r=r+1; lbl(sp,r,'Steerable order');
ud.efOrder = uieditfield(sp,'numeric','Value',gf(fp,'SteerableOrder',4),'Limits',[1 Inf]);
ud.efOrder.Layout.Row=r; ud.efOrder.Layout.Column=2;

r=r+1; lbl(sp,r,'Steerable sigma(s)');
ud.efSigma = uieditfield(sp,'text','Value',mat2str(gf(fp,'SigmaArray',[1 2])));
ud.efSigma.Layout.Row=r; ud.efSigma.Layout.Column=2;

r=r+1; lbl(sp,r,'Preview frame');
ud.efFrame = uieditfield(sp,'numeric','Value',max(1,round(MD.nFrames_/2)),'Limits',[1 MD.nFrames_],'RoundFractionalValues',true);
ud.efFrame.Layout.Row=r; ud.efFrame.Layout.Column=2;

r=r+1;
ud.btnPreview = uibutton(sp,'Text','Update preview','ButtonPushedFcn',@(~,~)doPreview(fig));
ud.btnPreview.Layout.Row=r; ud.btnPreview.Layout.Column=[1 2];

r=r+1;
ud.cbApplyAll = uicheckbox(sp,'Text','Apply settings to all movies','Value',false);
ud.cbApplyAll.Layout.Row=r; ud.cbApplyAll.Layout.Column=[1 2];

r=r+1;
ud.btnDone = uibutton(sp,'Text','Done (save to this movie)','ButtonPushedFcn',@(~,~)doApply(fig));
ud.btnDone.Layout.Row=r; ud.btnDone.Layout.Column=1;
ud.btnCancel = uibutton(sp,'Text','Cancel','ButtonPushedFcn',@(~,~)delete(fig));
ud.btnCancel.Layout.Row=r; ud.btnCancel.Layout.Column=2;

right = uipanel(g,'Title','Preview (cyan = body outline)');
rg = uigridlayout(right,[1 1]); rg.Padding=[6 6 6 6];
ud.ax = uiaxes(rg);
title(ud.ax,'Click "Update preview"'); ud.ax.XTick=[]; ud.ax.YTick=[];

fig.UserData = ud;
ud.slLevel.ValueChangedFcn = @(s,~) onLevelSlider(fig);
ud.efLevel.ValueChangedFcn = @(s,~) onLevelEdit(fig);
ud.slBright.ValueChangedFcn = @(s,~) doPreview(fig);
ud.ddMethod.ValueChangedFcn = @(s,~) doPreview(fig);
ud.ddChannel.ValueChangedFcn = @(s,~) onChannelChange(fig);
ud.efFrame.ValueChangedFcn = @(s,~) onChannelChange(fig);
fig.UserData = ud;

onChannelChange(fig);

if nargout>0, varargout{1}=fig; end
end

% ===================================================================
function onChannelChange(fig)
ud = fig.UserData;
fr = round(ud.efFrame.Value);
img = double(ud.MD.channels_(ud.ddChannel.Value).loadImage(fr));
mx = max(img(:)); mn = min(img(:));
if mx<=mn, mx=mn+1; end
ud.slLevel.Limits = [mn mx];

% restore a previously saved manual level the first time we open; otherwise
% only re-center if the current value falls outside the new channel's range.
savedML = [];
if isfield(ud.fp,'ManualThreshold') && ~isempty(ud.fp.ManualThreshold)
    savedML = ud.fp.ManualThreshold;
end
if ~isfield(ud,'levelInit') || ~ud.levelInit
    if ~isempty(savedML)
        v = min(max(savedML, mn), mx);
    else
        v = mn + 0.5*(mx-mn);
    end
    ud.slLevel.Value = v; ud.efLevel.Value = v;
    ud.levelInit = true;
elseif ud.slLevel.Value<mn || ud.slLevel.Value>mx
    v = mn + 0.5*(mx-mn);
    ud.slLevel.Value = v; ud.efLevel.Value = v;
end
fig.UserData = ud;
doPreview(fig);
end

function onLevelSlider(fig)
ud = fig.UserData;
ud.efLevel.Value = ud.slLevel.Value;
ud.ddMethod.Value = 'manual';
fig.UserData = ud;
doPreview(fig);
end

function onLevelEdit(fig)
ud = fig.UserData;
v = min(max(ud.efLevel.Value, ud.slLevel.Limits(1)), ud.slLevel.Limits(2));
ud.slLevel.Value = v;
ud.ddMethod.Value = 'manual';
fig.UserData = ud;
doPreview(fig);
end

function nudgeLevel(fig, dir)
% fine adjust: step = 0.5% of the slider range
ud = fig.UserData;
lim = ud.slLevel.Limits;
step = 0.005*(lim(2)-lim(1));
v = min(max(ud.slLevel.Value + dir*step, lim(1)), lim(2));
ud.slLevel.Value = v;
ud.efLevel.Value = v;
ud.ddMethod.Value = 'manual';
fig.UserData = ud;
doPreview(fig);
end

% ===================================================================
function doPreview(fig)
ud = fig.UserData;
fr = min(max(round(ud.efFrame.Value),1),ud.MD.nFrames_);
img = double(ud.MD.channels_(ud.ddChannel.Value).loadImage(fr));
p = collectParams(ud);
manual = [];
if strcmp(ud.ddMethod.Value,'manual'), manual = ud.efLevel.Value; end
try
    [bodyMask, level] = segmentFilopodiaBody(img, p, manual);
catch ME
    title(ud.ax, ['preview error: ' ME.message],'Interpreter','none');
    return;
end
imshow(img,[],'Parent',ud.ax); hold(ud.ax,'on');
% display-only contrast: cap the upper display limit so dim filopodia brighten
mn = min(img(:)); mx = max(img(:));
b = 1; try, b = ud.slBright.Value; catch, end
if mx>mn
    hi = mn + b*(mx-mn);
    if hi<=mn, hi = mn + 1e-6; end
    ud.ax.CLim = [mn hi];
end
B = bwboundaries(bodyMask);
for k=1:numel(B)
    plot(ud.ax,B{k}(:,2),B{k}(:,1),'c-','LineWidth',1.5);
end
hold(ud.ax,'off');
title(ud.ax, sprintf('frame %d   level=%.4g   body area=%d px', fr, level, nnz(bodyMask)));
ud.ax.XTick=[]; ud.ax.YTick=[];
end

% ===================================================================
function p = collectParams(ud)
p = ud.fp;
p.ChannelIndex      = ud.ddChannel.Value;
m = ud.ddMethod.Value;
if strcmp(m,'manual'), p.BodyThreshold = 'rosin'; else, p.BodyThreshold = m; end
p.ThreshScale       = ud.efScale.Value;
p.GaussianBlurSigma = ud.efBlur.Value;
p.BodyOpenRadius    = ud.efOpen.Value;
p.BodyClosingRadius = ud.efClose.Value;
p.BodyMinArea       = ud.efMinArea.Value;
p.SteerableOrder    = ud.efOrder.Value;
sig = str2num(ud.efSigma.Value); %#ok<ST2NM>
if ~isempty(sig), p.SigmaArray = sig; end
end

% ===================================================================
function doApply(fig)
ud = fig.UserData;
% Start from the process's CURRENT funParams_ (has all required hidden
% fields like OutputDirectory), then overwrite only what the GUI controls.
p = ud.proc.funParams_;
p.ChannelIndex      = ud.ddChannel.Value;
m = ud.ddMethod.Value;
if strcmp(m,'manual'), p.BodyThreshold = 'rosin'; else, p.BodyThreshold = m; end
p.ThreshScale       = ud.efScale.Value;
p.GaussianBlurSigma = ud.efBlur.Value;
p.BodyOpenRadius    = ud.efOpen.Value;
p.BodyClosingRadius = ud.efClose.Value;
p.BodyMinArea       = ud.efMinArea.Value;
p.SteerableOrder    = ud.efOrder.Value;
sig = str2num(ud.efSigma.Value); %#ok<ST2NM>
if ~isempty(sig), p.SigmaArray = sig; end
if strcmp(m,'manual')
    p.ManualThreshold = ud.efLevel.Value;
else
    p.ManualThreshold = [];
end

ok = false;
% setPara works on this install; use it as the primary path so the process
% records its parameters and dirty/updated state correctly.
try
    ud.proc.setPara(p);
    ok = true;
catch ME1
    try
        ud.proc.funParams_ = p;   % fallback
        ok = true;
    catch ME2
        if isvalid(fig)
            uialert(fig, sprintf('%s\n(fallback: %s)', ME1.message, ME2.message), ...
                'Error saving parameters');
        end
        return;
    end
end
if ~ok, return; end

try
    ud.MD.save();
catch ME
    if isvalid(fig), uialert(fig, ME.message, 'Saved params but MD.save failed'); end
end

% Tick the process checkbox / enable the next process in packageGUI the way
% the standard lccb settings windows do. processGUI_ApplyFcn reads its state
% from handles.figure1's UserData (crtProc/crtPackage/procID/MD/mainFig), so
% assemble that structure and call it.
    % Never let 'apply to all' propagate OutputDirectory: each movie must
    % keep its own output path. Strip it before passing to processGUI_ApplyFcn.
    safeParams = ud.proc.funParams_;
    if isfield(safeParams,'OutputDirectory'), safeParams = rmfield(safeParams,'OutputDirectory'); end
appliedOK = false;
try
    ud2.crtPackage  = ud.pkg;
    ud2.crtProc     = ud.proc;
    ud2.procID      = ud.procID;
    ud2.MD          = ud.MD;
    ud2.mainFig     = ud.mainFig;
    ud2.handles_main = guidata(ud.mainFig);   % main packageGUI handles (required by ApplyFcn)
    % Only create a new process on other movies if they don't already have one
    % of this class. Passing an empty procConstr prevents processGUI_ApplyFcn
    % from adding a duplicate when apply-to-all runs on a movie that already
    % has this process registered.
    ud2.procConstr = @(o,od) safeGetOrCreateProc(o, od, class(ud.proc));
    set(fig, 'UserData', mergeStruct(get(fig,'UserData'), ud2));
    handles.figure1 = fig;
    if isfield(ud,'cbApplyAll') && isvalid(ud.cbApplyAll), handles.checkbox_applytoall = ud.cbApplyAll; end
    processGUI_ApplyFcn(fig, [], handles, safeParams);
    appliedOK = true;
catch ME
    fprintf(2,'[filopodiaSegmentationProcessGUI] processGUI_ApplyFcn failed: %s\n', ME.message);
end
if ~appliedOK
    try
        h = guidata(ud.mainFig);
        packageGUI('refreshPackage_Callback', h.figure1, [], h);
    catch ME
        fprintf(2,'[filopodiaSegmentationProcessGUI] refresh fallback failed: %s\n', ME.message);
    end
end
if ishandle(ud.previewFig), delete(ud.previewFig); end
if isvalid(fig), delete(fig); end
end

% ===================================================================
function s = mergeStruct(s, add)
fn = fieldnames(add);
for i=1:numel(fn), s.(fn{i}) = add.(fn{i}); end
end

% ===================================================================
function lbl(parent,row,txt)
L = uilabel(parent,'Text',[txt ':'],'HorizontalAlignment','right');
L.Layout.Row=row; L.Layout.Column=1;
end
function v = gf(s,n,d), if isfield(s,n)&&~isempty(s.(n)), v=s.(n); else, v=d; end, end
function m = methodFromParams(fp)
if isfield(fp,'ManualThreshold') && ~isempty(fp.ManualThreshold)
    m = 'manual';
else
    bt = lower(num2str(gf(fp,'BodyThreshold','rosin')));
    if any(strcmp(bt,{'rosin','otsu'})), m = bt; else, m = 'manual'; end
end
end
function s = shortenPath(p)
parts = strsplit(p, filesep);
if numel(parts)>=2, s = fullfile(parts{end-1},parts{end}); else, s = p; end
end
