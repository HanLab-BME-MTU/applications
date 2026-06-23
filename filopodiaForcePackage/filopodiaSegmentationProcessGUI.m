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
sp = uigridlayout(left,[20 2]);
sp.ColumnWidth = {150,'1x'}; sp.RowHeight = repmat({26},1,20);
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
r=r+1; lbl(sp,r,'  (manual value)');
ud.efLevel = uieditfield(sp,'numeric','Value',0,'Tooltip','Used only when method = manual');
ud.efLevel.Layout.Row=r; ud.efLevel.Layout.Column=2;

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
if ud.slLevel.Value<mn || ud.slLevel.Value>mx
    ud.slLevel.Value = mn + 0.5*(mx-mn);
    ud.efLevel.Value = ud.slLevel.Value;
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
p = collectParams(ud);
if strcmp(ud.ddMethod.Value,'manual')
    p.ManualThreshold = ud.efLevel.Value;
else
    p.ManualThreshold = [];
end
try
    ud.proc.setPara(p);
    ud.MD.save();
    h = guidata(ud.mainFig);
    packageGUI('refreshPackage_Callback', h.figure1, [], h);
catch ME
    uialert(fig, ME.message, 'Error saving parameters');
    return;
end
if ishandle(ud.previewFig), delete(ud.previewFig); end
delete(fig);
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
