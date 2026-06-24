function varargout = filopodiaClassificationProcessGUI(varargin)
%FILOPODIACLASSIFICATIONPROCESSGUI  Settings window for Process 4 (filopodia
%classification) of the FilopodiaForcePackage. Channel + the main acceptance
%and shaft-tracing parameters, plus a preview that shows accepted filopodia
%(tip->base shaft, colored) versus rejected tip tracks (grey) on the talin
%frame. Classification uses the whole movie, so the preview reflects the
%SAVED result; change a parameter and re-run P4 to update it.
%
% Called by packageGUI:  filopodiaClassificationProcessGUI('mainFig',fig,procID)
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

W = 1000; H = 680;
fig = uifigure('Name',[proc.getName() ' Settings'],'Position',[120 120 W H],'Resize','on');
g = uigridlayout(fig,[1 2]); g.ColumnWidth = {380,'1x'}; g.RowHeight = {'1x'};

left = uipanel(g,'Title','Settings');
sp = uigridlayout(left,[20 2]);
sp.ColumnWidth = {190,'1x'}; sp.RowHeight = repmat({26},1,20);
sp.Padding=[10 10 10 10]; sp.RowSpacing=4;

ud = struct();
ud.mainFig=mainFig; ud.procID=procID; ud.MD=MD; ud.pkg=pkg; ud.proc=proc; ud.fp=fp;
ud.previewFig = -1;

r=0;
r=r+1; lbl(sp,r,'Channel');
chanNames = MD.getChannelPaths(); shortNames = cellfun(@(s) shortenPath(s), chanNames,'unif',false);
ci = fp.ChannelIndex; if isempty(ci)||ci<1, ci=1; end
ud.ddChannel = uidropdown(sp,'Items',shortNames,'ItemsData',1:numel(chanNames),'Value',min(ci,numel(chanNames)));
ud.ddChannel.Layout.Row=r; ud.ddChannel.Layout.Column=2;

r=r+1; lbl(sp,r,'Min tip lifetime (frames)');
ud.efLife = uieditfield(sp,'numeric','Value',gf(fp,'MinTipLifetime',5),'Limits',[1 Inf], ...
    'Tooltip','Tip track must persist at least this many frames');
ud.efLife.Layout.Row=r; ud.efLife.Layout.Column=2;

r=r+1; lbl(sp,r,'Min linear fraction');
ud.efLin = uieditfield(sp,'numeric','Value',gf(fp,'MinLinearFrac',0.85),'Limits',[0 1], ...
    'Tooltip','Fraction of trajectory variance on principal axis to count as linear');
ud.efLin.Layout.Row=r; ud.efLin.Layout.Column=2;

r=r+1; lbl(sp,r,'Min tip dist (px)');
ud.efTipDist = uieditfield(sp,'numeric','Value',gf(fp,'MinTipDist',6),'Limits',[0 Inf], ...
    'Tooltip','Tip must reach at least this far outside the body');
ud.efTipDist.Layout.Row=r; ud.efTipDist.Layout.Column=2;

r=r+1; lbl(sp,r,'Shaft band (px)');
ud.efShaftBand = uieditfield(sp,'numeric','Value',gf(fp,'ShaftBand',4),'Limits',[0 Inf]);
ud.efShaftBand.Layout.Row=r; ud.efShaftBand.Layout.Column=2;

r=r+1; lbl(sp,r,'Max shaft len (px)');
ud.efMaxShaft = uieditfield(sp,'numeric','Value',gf(fp,'MaxShaftLen',160),'Limits',[0 Inf]);
ud.efMaxShaft.Layout.Row=r; ud.efMaxShaft.Layout.Column=2;

r=r+1; lbl(sp,r,'Body max angle (deg)');
ud.efBodyAng = uieditfield(sp,'numeric','Value',gf(fp,'BodyMaxAngle',75),'Limits',[0 180], ...
    'Tooltip','Shaft must point within this angle of body-ward');
ud.efBodyAng.Layout.Row=r; ud.efBodyAng.Layout.Column=2;

r=r+1; lbl(sp,r,'Min reach fraction');
ud.efReach = uieditfield(sp,'numeric','Value',gf(fp,'MinReachFrac',0.5),'Limits',[0 1], ...
    'Tooltip','Accept tip track if it acts as a tip in >= this fraction of frames');
ud.efReach.Layout.Row=r; ud.efReach.Layout.Column=2;

r=r+1; lbl(sp,r,'Min base separation (px)');
ud.efBaseSep = uieditfield(sp,'numeric','Value',gf(fp,'MinBaseSep',8),'Limits',[0 Inf]);
ud.efBaseSep.Layout.Row=r; ud.efBaseSep.Layout.Column=2;

r=r+1; lbl(sp,r,'Length penalty');
ud.efLenPen = uieditfield(sp,'numeric','Value',gf(fp,'LenPenalty',0.6),'Limits',[0 Inf], ...
    'Tooltip','Higher prefers shorter, more radial shafts');
ud.efLenPen.Layout.Row=r; ud.efLenPen.Layout.Column=2;

r=r+1; lbl(sp,r,'Display brightness');
ud.slBright = uislider(sp,'Limits',[0.05 1],'Value',1);
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

right = uipanel(g,'Title','Live preview (colored = traced tip\rightarrowbase; grey = no base found)');
rg = uigridlayout(right,[1 1]); rg.Padding=[6 6 6 6];
ud.ax = uiaxes(rg); title(ud.ax,'Click "Update preview"'); ud.ax.XTick=[]; ud.ax.YTick=[];

fig.UserData = ud;
ud.slBright.ValueChangedFcn = @(s,~) doPreview(fig);
ud.ddChannel.ValueChangedFcn = @(s,~) doPreview(fig);
ud.efFrame.ValueChangedFcn = @(s,~) doPreview(fig);
% shaft-tracing params recompute the preview live (true tuning)
ud.efMaxShaft.ValueChangedFcn = @(s,~) doPreview(fig);
ud.efBodyAng.ValueChangedFcn  = @(s,~) doPreview(fig);
ud.efTipDist.ValueChangedFcn  = @(s,~) doPreview(fig);
ud.efShaftBand.ValueChangedFcn= @(s,~) doPreview(fig);
fig.UserData = ud;

doPreview(fig);
if nargout>0, varargout{1}=fig; end
end

% ===================================================================
function doPreview(fig)
ud = fig.UserData;
fr = min(max(round(ud.efFrame.Value),1),ud.MD.nFrames_);
iChan = ud.ddChannel.Value;
img = double(ud.MD.channels_(iChan).loadImage(fr));
p = collectParams(ud);

imshow(img,[],'Parent',ud.ax); hold(ud.ax,'on');
mn=min(img(:)); mx=max(img(:)); b=1; try, b=ud.slBright.Value; catch, end
if mx>mn, hi=mn+b*(mx-mn); if hi<=mn, hi=mn+1e-6; end, ud.ax.CLim=[mn hi]; end

% body outline for this frame
bodyMask = loadBodyMask(ud.MD, fr, iChan);
if ~isempty(bodyMask)&&any(bodyMask(:))
    B=bwboundaries(bodyMask); for k=1:numel(B), plot(ud.ax,B{k}(:,2),B{k}(:,1),'c-','LineWidth',0.8); end
end

% tips present in this frame, taken from the P3 tracks (live, not from a
% previous P4 run), then shaft-traced with the CURRENT parameters.
tipXY = loadTipsAtFrame(ud.MD, iChan, fr);
if isempty(tipXY)
    hold(ud.ax,'off');
    title(ud.ax, sprintf('frame %d   (no tracks: run P3 first)', fr));
    ud.ax.XTick=[]; ud.ax.YTick=[]; return;
end

filo = classifyFilopodiaPreviewFrame(tipXY, bodyMask, p);
cmap = lines(7); nAcc=0; nRej=0;
for f = 1:numel(filo)
    tp = filo(f).tipPos;
    if filo(f).accepted
        nAcc = nAcc + 1; c = cmap(mod(f-1,7)+1,:);
        bp = filo(f).basePos;
        plot(ud.ax,[bp(1) tp(1)],[bp(2) tp(2)],'-','Color',c,'LineWidth',1.5);
        plot(ud.ax,tp(1),tp(2),'.','Color',c,'MarkerSize',12);
        plot(ud.ax,bp(1),bp(2),'o','Color',c,'MarkerSize',4,'LineWidth',1);
    else
        nRej = nRej + 1;
        plot(ud.ax,tp(1),tp(2),'.','Color',[0.6 0.6 0.6],'MarkerSize',9);
    end
end
hold(ud.ax,'off');
title(ud.ax, sprintf('frame %d   traced=%d   no-base=%d   (live preview)', fr, nAcc, nRej));
ud.ax.XTick=[]; ud.ax.YTick=[];
end

% ===================================================================
function tipXY = loadTipsAtFrame(MD, iChan, fr)
% tip candidate positions in this frame = all tracked positions present at fr
tipXY = zeros(0,2);
try
    ip = MD.getProcessIndex('FilopodiaTrackingProcess',1,0);
    if isempty(ip), return; end
    pr = MD.processes_{ip};
    f = pr.outFilePaths_{1,iChan};
    if isempty(f)||exist(f,'file')~=2, f = fullfile(pr.funParams_.OutputDirectory,'filoTracks.mat'); end
    if exist(f,'file')~=2, return; end
    vars = who('-file', f);
    if ~ismember('adhesionTracks', vars), return; end
    S = load(f,'adhesionTracks'); T = S.adhesionTracks;
    for t = 1:numel(T)
        j = find(T(t).frames==fr,1);
        if ~isempty(j), tipXY(end+1,:) = T(t).pos(j,:); end %#ok<AGROW>
    end
catch
end
end

% ===================================================================
function bodyMask = loadBodyMask(MD, fr, iChan)
bodyMask = [];
try
    ip = MD.getProcessIndex('FilopodiaSegmentationProcess',1,0);
    if isempty(ip), return; end
    pr = MD.processes_{ip};
    fn = fullfile(pr.funParams_.OutputDirectory, sprintf('filoSeg_frame_%04d.mat', fr));
    if exist(fn,'file')==2, S=load(fn,'bodyMask'); bodyMask=S.bodyMask; return; end
    img = double(MD.channels_(iChan).loadImage(fr));
    bodyMask = segmentFilopodiaBody(img, pr.funParams_);
catch
end
end

% ===================================================================
function p = collectParams(ud)
p = ud.proc.funParams_;
p.ChannelIndex   = ud.ddChannel.Value;
p.MinTipLifetime = ud.efLife.Value;
p.MinLinearFrac  = ud.efLin.Value;
p.MinTipDist     = ud.efTipDist.Value;
p.ShaftBand      = ud.efShaftBand.Value;
p.MaxShaftLen    = ud.efMaxShaft.Value;
p.BodyMaxAngle   = ud.efBodyAng.Value;
p.MinReachFrac   = ud.efReach.Value;
p.MinBaseSep     = ud.efBaseSep.Value;
p.LenPenalty     = ud.efLenPen.Value;
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

appliedOK=false;
try
    ud2.crtPackage=ud.pkg; ud2.crtProc=ud.proc; ud2.procID=ud.procID;
    ud2.MD=ud.MD; ud2.mainFig=ud.mainFig; ud2.handles_main=guidata(ud.mainFig);
    set(fig,'UserData',mergeStruct(get(fig,'UserData'),ud2));
    handles.figure1=fig;
    processGUI_ApplyFcn(fig, [], handles, ud.proc.funParams_);
    appliedOK=true;
catch ME
    fprintf(2,'[filopodiaClassificationProcessGUI] processGUI_ApplyFcn failed: %s\n', ME.message);
end
if ~appliedOK
    try, h=guidata(ud.mainFig); packageGUI('refreshPackage_Callback', h.figure1, [], h);
    catch ME, fprintf(2,'[filopodiaClassificationProcessGUI] refresh fallback failed: %s\n', ME.message); end
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
