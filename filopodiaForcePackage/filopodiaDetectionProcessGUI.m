function varargout = filopodiaDetectionProcessGUI(varargin)
%FILOPODIADETECTIONPROCESSGUI  Settings window for Process 2 (tip/base
%detection) of the FilopodiaForcePackage. Channel selection + interactive
%preview: adjust PSF sigma / alpha / distance gates and see detected tips
%(yellow) and bases (green) overlaid on the talin frame live.
%
% Called by packageGUI:  filopodiaDetectionProcessGUI('mainFig',fig,procID)
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

W = 980; H = 640;
fig = uifigure('Name',[proc.getName() ' Settings'],'Position',[120 120 W H],'Resize','on');
g = uigridlayout(fig,[1 2]);
g.ColumnWidth = {360,'1x'}; g.RowHeight = {'1x'};

left = uipanel(g,'Title','Settings');
sp = uigridlayout(left,[20 2]);
sp.ColumnWidth = {160,'1x'}; sp.RowHeight = repmat({26},1,20);
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

r=r+1; lbl(sp,r,'Detect mode');
ud.ddMode = uidropdown(sp,'Items',{'auto','all','tip'},'Value',gf(fp,'DetectMode','auto'));
ud.ddMode.Layout.Row=r; ud.ddMode.Layout.Column=2;

r=r+1; lbl(sp,r,'PSF sigma (px)');
ud.efSigma = uieditfield(sp,'numeric','Value',gf(fp,'PSFsigma',2.1),'Limits',[0.2 Inf], ...
    'Tooltip','Tip detection scale; lower = more sensitive to small tips');
ud.efSigma.Layout.Row=r; ud.efSigma.Layout.Column=2;

r=r+1; lbl(sp,r,'Alpha');
ud.efAlpha = uieditfield(sp,'numeric','Value',gf(fp,'Alpha',0.05),'Limits',[1e-6 1], ...
    'Tooltip','pointSourceDetection significance; higher = more (looser) detections');
ud.efAlpha.Layout.Row=r; ud.efAlpha.Layout.Column=2;

r=r+1; lbl(sp,r,'Tip max dist (px)');
ud.efTipMax = uieditfield(sp,'numeric','Value',gf(fp,'TipMaxDistFromBody',130),'Limits',[0 Inf], ...
    'Tooltip','Max distance outside body for a punctum to count as a tip');
ud.efTipMax.Layout.Row=r; ud.efTipMax.Layout.Column=2;

r=r+1; lbl(sp,r,'Base search band (px)');
ud.efBaseBand = uieditfield(sp,'numeric','Value',gf(fp,'BaseSearchBand',5),'Limits',[0 Inf]);
ud.efBaseBand.Layout.Row=r; ud.efBaseBand.Layout.Column=2;

r=r+1; lbl(sp,r,'Base inside band (px)');
ud.efInBand = uieditfield(sp,'numeric','Value',gf(fp,'BaseInsideBand',4),'Limits',[0 Inf]);
ud.efInBand.Layout.Row=r; ud.efInBand.Layout.Column=2;

r=r+1; lbl(sp,r,'Max tip-base dist (px)');
ud.efMaxTB = uieditfield(sp,'numeric','Value',gf(fp,'MaxTipBaseDist',160),'Limits',[0 Inf], ...
    'Tooltip','Max shaft trace length (batch run only; not shown in preview)');
ud.efMaxTB.Layout.Row=r; ud.efMaxTB.Layout.Column=2;

r=r+1; lbl(sp,r,'Min filo length (px)');
ud.efMinLen = uieditfield(sp,'numeric','Value',gf(fp,'MinFiloLength',5),'Limits',[0 Inf]);
ud.efMinLen.Layout.Row=r; ud.efMinLen.Layout.Column=2;

r=r+1; lbl(sp,r,'Display brightness');
ud.slBright = uislider(sp,'Limits',[0.05 1],'Value',1, ...
    'Tooltip','Display only: lower = brighter dim features (does not affect detection)');
ud.slBright.Layout.Row=r; ud.slBright.Layout.Column=2;

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

right = uipanel(g,'Title','Preview (yellow = tips, green = bases, cyan = body)');
rg = uigridlayout(right,[1 1]); rg.Padding=[6 6 6 6];
ud.ax = uiaxes(rg);
title(ud.ax,'Click "Update preview"'); ud.ax.XTick=[]; ud.ax.YTick=[];

fig.UserData = ud;
ud.slBright.ValueChangedFcn = @(s,~) doPreview(fig);
ud.ddChannel.ValueChangedFcn = @(s,~) doPreview(fig);
ud.efFrame.ValueChangedFcn = @(s,~) doPreview(fig);
fig.UserData = ud;

doPreview(fig);
if nargout>0, varargout{1}=fig; end
end

% ===================================================================
function doPreview(fig)
ud = fig.UserData;
fr = min(max(round(ud.efFrame.Value),1),ud.MD.nFrames_);
img = double(ud.MD.channels_(ud.ddChannel.Value).loadImage(fr));
p = collectParams(ud);

bodyMask = loadBodyMask(ud.MD, fr, ud.ddChannel.Value);

try
    [tipXY, baseXY] = detectFilopodiaPointsPreview(img, bodyMask, p);
catch ME
    title(ud.ax, ['preview error: ' ME.message],'Interpreter','none'); return;
end

imshow(img,[],'Parent',ud.ax); hold(ud.ax,'on');
mn = min(img(:)); mx = max(img(:)); b = 1; try, b = ud.slBright.Value; catch, end
if mx>mn
    hi = mn + b*(mx-mn); if hi<=mn, hi=mn+1e-6; end
    ud.ax.CLim = [mn hi];
end
if ~isempty(bodyMask) && any(bodyMask(:))
    B = bwboundaries(bodyMask);
    for k=1:numel(B), plot(ud.ax,B{k}(:,2),B{k}(:,1),'c-','LineWidth',1.2); end
end
if ~isempty(baseXY), plot(ud.ax,baseXY(:,1),baseXY(:,2),'g.','MarkerSize',10); end
if ~isempty(tipXY),  plot(ud.ax,tipXY(:,1), tipXY(:,2), 'y.','MarkerSize',12); end
hold(ud.ax,'off');
title(ud.ax, sprintf('frame %d   tips=%d   bases=%d', fr, size(tipXY,1), size(baseXY,1)));
ud.ax.XTick=[]; ud.ax.YTick=[];
end

% ===================================================================
function bodyMask = loadBodyMask(MD, fr, iChan)
% prefer the saved P1 body mask; fall back to computing it on the fly
bodyMask = [];
try
    ip = MD.getProcessIndex('FilopodiaSegmentationProcess',1,0);
    if ~isempty(ip)
        pr = MD.processes_{ip};
        od = pr.funParams_.OutputDirectory;
        fn = fullfile(od, sprintf('filoSeg_frame_%04d.mat', fr));
        if exist(fn,'file')==2
            S = load(fn,'bodyMask'); bodyMask = S.bodyMask; return;
        end
    end
catch
end
try
    ip = MD.getProcessIndex('FilopodiaSegmentationProcess',1,0);
    if ~isempty(ip)
        img = double(MD.channels_(iChan).loadImage(fr));
        bodyMask = segmentFilopodiaBody(img, MD.processes_{ip}.funParams_);
    end
catch
end
end

% ===================================================================
function p = collectParams(ud)
p = ud.fp;
p.ChannelIndex       = ud.ddChannel.Value;
p.DetectMode         = ud.ddMode.Value;
p.PSFsigma           = ud.efSigma.Value;
p.Alpha              = ud.efAlpha.Value;
p.TipMaxDistFromBody = ud.efTipMax.Value;
p.BaseSearchBand     = ud.efBaseBand.Value;
p.BaseInsideBand     = ud.efInBand.Value;
p.MaxTipBaseDist     = ud.efMaxTB.Value;
p.MinFiloLength      = ud.efMinLen.Value;
end

% ===================================================================
function doApply(fig)
ud = fig.UserData;
p = ud.proc.funParams_;
p.ChannelIndex       = ud.ddChannel.Value;
p.DetectMode         = ud.ddMode.Value;
p.PSFsigma           = ud.efSigma.Value;
p.Alpha              = ud.efAlpha.Value;
p.TipMaxDistFromBody = ud.efTipMax.Value;
p.BaseSearchBand     = ud.efBaseBand.Value;
p.BaseInsideBand     = ud.efInBand.Value;
p.MaxTipBaseDist     = ud.efMaxTB.Value;
p.MinFiloLength      = ud.efMinLen.Value;

ok = false;
try, ud.proc.setPara(p); ok = true;
catch ME1
    try, ud.proc.funParams_ = p; ok = true;
    catch ME2
        if isvalid(fig), uialert(fig, sprintf('%s\n(fallback: %s)',ME1.message,ME2.message),'Error saving parameters'); end
        return;
    end
end
if ~ok, return; end
try, ud.MD.save(); catch ME
    if isvalid(fig), uialert(fig, ME.message, 'Saved params but MD.save failed'); end
end

appliedOK = false;
try
    ud2.crtPackage  = ud.pkg; ud2.crtProc = ud.proc; ud2.procID = ud.procID;
    ud2.MD = ud.MD; ud2.mainFig = ud.mainFig; ud2.handles_main = guidata(ud.mainFig);
    set(fig,'UserData', mergeStruct(get(fig,'UserData'),ud2));
    handles.figure1 = fig;
    processGUI_ApplyFcn(fig, [], handles, ud.proc.funParams_);
    appliedOK = true;
catch ME
    fprintf(2,'[filopodiaDetectionProcessGUI] processGUI_ApplyFcn failed: %s\n', ME.message);
end
if ~appliedOK
    try, h = guidata(ud.mainFig); packageGUI('refreshPackage_Callback', h.figure1, [], h);
    catch ME, fprintf(2,'[filopodiaDetectionProcessGUI] refresh fallback failed: %s\n', ME.message); end
end
if ishandle(ud.previewFig), delete(ud.previewFig); end
if isvalid(fig), delete(fig); end
end

% ===================================================================
function lbl(parent,row,txt)
L = uilabel(parent,'Text',[txt ':'],'HorizontalAlignment','right');
L.Layout.Row=row; L.Layout.Column=1;
end
function v = gf(s,n,d), if isfield(s,n)&&~isempty(s.(n)), v=s.(n); else, v=d; end, end
function s = mergeStruct(s, add)
fn = fieldnames(add); for i=1:numel(fn), s.(fn{i}) = add.(fn{i}); end
end
function s = shortenPath(p)
parts = strsplit(p, filesep);
if numel(parts)>=2, s = fullfile(parts{end-1},parts{end}); else, s = p; end
end