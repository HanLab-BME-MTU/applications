function varargout = filopodiaSamplingProcessGUI(varargin)
%FILOPODIASAMPLINGPROCESSGUI  Settings window for Process 5 (force / intensity
%sampling) of the FilopodiaForcePackage. Channel + sampling parameters, plus
%a preview that shows, for a chosen filopodium and frame, the shaft profile
%(traction and talin vs arc length tip->base) from the SAVED P5 result. P5
%samples the TFM force field along each shaft, so the preview reflects the
%last run; change a parameter and re-run P5 to update it.
%
% Called by packageGUI:  filopodiaSamplingProcessGUI('mainFig',fig,procID)
% Sangyoon J. Han / 2026

if numel(varargin) < 3
    error('Call from packageGUI only: crtProcGUI(''mainFig'',mainFig,procID)');
end
mainFig = varargin{2}; procID = varargin{3};
ud_main = get(mainFig,'UserData');
pkg = ud_main.crtPackage; MD = pkg.getOwner();
proc = pkg.getProcess(procID);
if isempty(proc)
    constrs = pkg.getDefaultProcessConstructors(procID);
    proc = constrs{1}(MD, MD.outputDirectory_); MD.addProcess(proc);
end
fp = proc.funParams_;

W = 1000; H = 640;
fig = uifigure('Name',[proc.getName() ' Settings'],'Position',[120 120 W H],'Resize','on');
g = uigridlayout(fig,[1 2]); g.ColumnWidth = {360,'1x'}; g.RowHeight = {'1x'};

left = uipanel(g,'Title','Settings');
sp = uigridlayout(left,[17 2]);
sp.ColumnWidth = {170,'1x'}; sp.RowHeight = repmat({26},1,17);
sp.Padding=[10 10 10 10]; sp.RowSpacing=4;

ud = struct();
ud.mainFig=mainFig; ud.procID=procID; ud.MD=MD; ud.pkg=pkg; ud.proc=proc; ud.fp=fp;
ud.previewFig=-1;

r=0;
r=r+1; lbl(sp,r,'Channel');
chanNames = MD.getChannelPaths(); shortNames = cellfun(@(s) shortenPath(s), chanNames,'unif',false);
ci=fp.ChannelIndex; if isempty(ci)||ci<1, ci=1; end
ud.ddChannel = uidropdown(sp,'Items',shortNames,'ItemsData',1:numel(chanNames),'Value',min(ci,numel(chanNames)));
ud.ddChannel.Layout.Row=r; ud.ddChannel.Layout.Column=2;

r=r+1; lbl(sp,r,'Force package');
ud.efForcePkg = uieditfield(sp,'text','Value',gf(fp,'ForcePackageName','TFMPackage'));
ud.efForcePkg.Layout.Row=r; ud.efForcePkg.Layout.Column=2;

r=r+1; lbl(sp,r,'Force process');
ud.efForceProc = uieditfield(sp,'text','Value',gf(fp,'ForceProcessName','ForceFieldCalculationProcess'));
ud.efForceProc.Layout.Row=r; ud.efForceProc.Layout.Column=2;

r=r+1; lbl(sp,r,'Shaft sample step (px)');
ud.efStep = uieditfield(sp,'numeric','Value',gf(fp,'ShaftSampleStep',3),'Limits',[0.5 Inf]);
ud.efStep.Layout.Row=r; ud.efStep.Layout.Column=2;

r=r+1; lbl(sp,r,'Sample radius (px)');
ud.efRad = uieditfield(sp,'numeric','Value',gf(fp,'SampleRadius',1),'Limits',[0 Inf]);
ud.efRad.Layout.Row=r; ud.efRad.Layout.Column=2;

r=r+1; lbl(sp,r,'Norm profile points');
ud.efNorm = uieditfield(sp,'numeric','Value',gf(fp,'NormProfileN',50),'Limits',[5 Inf]);
ud.efNorm.Layout.Row=r; ud.efNorm.Layout.Column=2;

r=r+1; lbl(sp,r,'Sample statistic');
ud.ddStat = uidropdown(sp,'Items',{'mean','median','max'},'Value',gf(fp,'SampleStat','mean'));
ud.ddStat.Layout.Row=r; ud.ddStat.Layout.Column=2;

r=r+1; lbl(sp,r,'Preview filopodium #');
ud.efFilo = uieditfield(sp,'numeric','Value',1,'Limits',[1 Inf],'RoundFractionalValues',true);
ud.efFilo.Layout.Row=r; ud.efFilo.Layout.Column=2;

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

right = uipanel(g,'Title','Shaft profile (traction & talin vs arc length, tip\rightarrowbase)');
rg = uigridlayout(right,[1 1]); rg.Padding=[6 6 6 6];
ud.ax = uiaxes(rg);
title(ud.ax,'Click "Update preview" (needs a saved P5 run)');

fig.UserData = ud;
ud.ddChannel.ValueChangedFcn = @(s,~) doPreview(fig);
ud.efFilo.ValueChangedFcn = @(s,~) doPreview(fig);
ud.efFrame.ValueChangedFcn = @(s,~) doPreview(fig);
fig.UserData = ud;

doPreview(fig);
if nargout>0, varargout{1}=fig; end
end

% ===================================================================
function doPreview(fig)
ud = fig.UserData;
iChan = ud.ddChannel.Value;
fNo = round(ud.efFilo.Value);
fr  = round(ud.efFrame.Value);

[fs, ok] = loadSamples(ud.MD, iChan);
cla(ud.ax,'reset');
if ~ok || isempty(fs)
    title(ud.ax,'No saved P5 result yet: set parameters, Done, then run P5.'); return;
end
nFilo = numel(fs);
fNo = min(max(fNo,1),nFilo);
S = fs(fNo);
% pick the requested frame's shaft profile, else the middle valid one
sp = [];
if isfield(S,'shaftProfile') && ~isempty(S.shaftProfile)
    j = find(S.frames==fr,1);
    if isempty(j)
        valid = find(~cellfun(@isempty,S.shaftProfile),1,'first');
        if ~isempty(valid), j = valid; end
    end
    if ~isempty(j) && j<=numel(S.shaftProfile), sp = S.shaftProfile{j}; end
end
if isempty(sp) || isempty(sp.s_nm)
    title(ud.ax, sprintf('Filopodium %d: no shaft profile at frame %d', fNo, fr)); return;
end

yyaxis(ud.ax,'left');
plot(ud.ax, sp.s_nm, sp.force, '-o','LineWidth',1.5); ylabel(ud.ax,'traction (Pa)');
yyaxis(ud.ax,'right');
plot(ud.ax, sp.s_nm, sp.talin, '-s','LineWidth',1.5); ylabel(ud.ax,'talin (a.u.)');
xlabel(ud.ax,'arc length tip\rightarrowbase (nm)');
title(ud.ax, sprintf('Filopodium %d / %d   frame %d   (trackId %d)', fNo, nFilo, fr, S.tipTrackId));
grid(ud.ax,'on');
end

% ===================================================================
function [fs, ok] = loadSamples(MD, iChan)
fs = []; ok = false;
try
    ip = MD.getProcessIndex('FilopodiaSamplingProcess',1,0);
    if isempty(ip), return; end
    pr = MD.processes_{ip};
    f = pr.outFilePaths_{1,iChan};
    if isempty(f)||exist(f,'file')~=2, f = fullfile(pr.funParams_.OutputDirectory,'filoSamples.mat'); end
    if exist(f,'file')~=2, return; end
    S = load(f,'filoSamples'); fs = S.filoSamples; ok = true;
catch
end
end

% ===================================================================
function p = collectParams(ud)
p = ud.proc.funParams_;
p.ChannelIndex     = ud.ddChannel.Value;
p.ForcePackageName = ud.efForcePkg.Value;
p.ForceProcessName = ud.efForceProc.Value;
p.ShaftSampleStep  = ud.efStep.Value;
p.SampleRadius     = ud.efRad.Value;
p.NormProfileN     = ud.efNorm.Value;
p.SampleStat       = ud.ddStat.Value;
end

% ===================================================================
function doApply(fig)
ud = fig.UserData;
p = collectParams(ud);
ok=false;
try, ud.proc.setPara(p); ok=true;
catch ME1
    try, ud.proc.funParams_=p; ok=true;
    catch ME2
        if isvalid(fig), uialert(fig, sprintf('%s\n(fallback: %s)',ME1.message,ME2.message),'Error saving parameters'); end
        return;
    end
end
if ~ok, return; end
try, ud.MD.save(); catch ME
    if isvalid(fig), uialert(fig, ME.message, 'Saved params but MD.save failed'); end
end
    % Never let 'apply to all' propagate OutputDirectory: each movie must
    % keep its own output path. Strip it before passing to processGUI_ApplyFcn.
    safeParams = ud.proc.funParams_;
    if isfield(safeParams,'OutputDirectory'), safeParams = rmfield(safeParams,'OutputDirectory'); end
appliedOK=false;
try
    ud2.crtPackage=ud.pkg; ud2.crtProc=ud.proc; ud2.procID=ud.procID;
    ud2.MD=ud.MD; ud2.mainFig=ud.mainFig; ud2.handles_main=guidata(ud.mainFig);
    % Only create a new process on other movies if they don't already have one
    % of this class. Passing an empty procConstr prevents processGUI_ApplyFcn
    % from adding a duplicate when apply-to-all runs on a movie that already
    % has this process registered.
    % Do not pass procConstr: prevents ApplyFcn from creating duplicate processes.
    % Params are applied to existing processes only (parseProcessParams).
    set(fig,'UserData',mergeStruct(get(fig,'UserData'),ud2));
    handles.figure1 = fig;
    if isfield(ud,'cbApplyAll') && isvalid(ud.cbApplyAll), handles.checkbox_applytoall = ud.cbApplyAll; end
    processGUI_ApplyFcn(fig, [], handles, safeParams);
    appliedOK=true;
catch ME
    fprintf(2,'[filopodiaSamplingProcessGUI] processGUI_ApplyFcn failed: %s\n', ME.message);
end
if ~appliedOK
    try, h=guidata(ud.mainFig); packageGUI('refreshPackage_Callback', h.figure1, [], h);
    catch ME, fprintf(2,'[filopodiaSamplingProcessGUI] refresh fallback failed: %s\n', ME.message); end
end
if ishandle(ud.previewFig), delete(ud.previewFig); end
if isvalid(fig), delete(fig); end
end

% ===================================================================
function lbl(parent,row,txt)
L=uilabel(parent,'Text',[txt ':'],'HorizontalAlignment','right'); L.Layout.Row=row; L.Layout.Column=1;
end
function v = gf(s,n,d), if isfield(s,n)&&~isempty(s.(n)), v=s.(n); else, v=d; end, end
function s = mergeStruct(s, add), fn=fieldnames(add); for i=1:numel(fn), s.(fn{i})=add.(fn{i}); end, end
function s = shortenPath(p)
parts=strsplit(p,filesep); if numel(parts)>=2, s=fullfile(parts{end-1},parts{end}); else, s=p; end
end
