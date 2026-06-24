function varargout = filopodiaTrackingProcessGUI(varargin)
%FILOPODIATRACKINGPROCESSGUI  Settings window for Process 3 (tip tracking) of
%the FilopodiaForcePackage. Channel selection + parameters, plus a preview
%that shows the existing tracks (one color per track, with a short trailing
%trajectory) overlaid on the talin frame. Tracking links the whole movie, so
%the preview reflects the SAVED result; change a parameter and re-run P3 to
%update it. Before P3 has run, the preview shows the raw P2 detections.
%
% Called by packageGUI:  filopodiaTrackingProcessGUI('mainFig',fig,procID)
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
sp = uigridlayout(left,[18 2]);
sp.ColumnWidth = {160,'1x'}; sp.RowHeight = repmat({26},1,18);
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

r=r+1; lbl(sp,r,'Max link dist (px/frame)');
ud.efLink = uieditfield(sp,'numeric','Value',gf(fp,'MaxLinkDist',8),'Limits',[0 Inf], ...
    'Tooltip','Max distance a tip can move between consecutive frames');
ud.efLink.Layout.Row=r; ud.efLink.Layout.Column=2;

r=r+1; lbl(sp,r,'Link using base');
ud.cbBase = uicheckbox(sp,'Text','','Value',logical(gf(fp,'LinkUseBase',true)));
ud.cbBase.Layout.Row=r; ud.cbBase.Layout.Column=2;

r=r+1; lbl(sp,r,'Max gap frames');
ud.efGap = uieditfield(sp,'numeric','Value',gf(fp,'MaxGapFrames',3),'Limits',[0 Inf], ...
    'Tooltip','Bridge disappearances up to this many frames');
ud.efGap.Layout.Row=r; ud.efGap.Layout.Column=2;

r=r+1; lbl(sp,r,'Min track length (frames)');
ud.efMinTrk = uieditfield(sp,'numeric','Value',gf(fp,'MinTrackLength',3),'Limits',[1 Inf]);
ud.efMinTrk.Layout.Row=r; ud.efMinTrk.Layout.Column=2;

r=r+1; lbl(sp,r,'Vel smooth window');
ud.efVel = uieditfield(sp,'numeric','Value',gf(fp,'VelSmoothWin',3),'Limits',[1 Inf]);
ud.efVel.Layout.Row=r; ud.efVel.Layout.Column=2;

r=r+1; lbl(sp,r,'Trajectory trail (frames)');
ud.efTrail = uieditfield(sp,'numeric','Value',10,'Limits',[0 Inf], ...
    'Tooltip','Preview only: how many past frames of each track to draw');
ud.efTrail.Layout.Row=r; ud.efTrail.Layout.Column=2;

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

right = uipanel(g,'Title','Preview (colored = tracks; white dot = current tip)');
rg = uigridlayout(right,[1 1]); rg.Padding=[6 6 6 6];
ud.ax = uiaxes(rg);
title(ud.ax,'Click "Update preview"'); ud.ax.XTick=[]; ud.ax.YTick=[];

fig.UserData = ud;
ud.slBright.ValueChangedFcn = @(s,~) doPreview(fig);
ud.ddChannel.ValueChangedFcn = @(s,~) doPreview(fig);
ud.efFrame.ValueChangedFcn = @(s,~) doPreview(fig);
ud.efTrail.ValueChangedFcn = @(s,~) doPreview(fig);
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

imshow(img,[],'Parent',ud.ax); hold(ud.ax,'on');
mn=min(img(:)); mx=max(img(:)); b=1; try, b=ud.slBright.Value; catch, end
if mx>mn, hi=mn+b*(mx-mn); if hi<=mn, hi=mn+1e-6; end, ud.ax.CLim=[mn hi]; end

tracks = loadTracks(ud.MD, iChan);
if isempty(tracks)
    % no tracking yet: show raw detections as a hint
    [tipXY, baseXY] = loadDetections(ud.MD, iChan, fr);
    if ~isempty(baseXY), plot(ud.ax,baseXY(:,1),baseXY(:,2),'g.','MarkerSize',8); end
    if ~isempty(tipXY),  plot(ud.ax,tipXY(:,1), tipXY(:,2),'y.','MarkerSize',10); end
    hold(ud.ax,'off');
    title(ud.ax, sprintf('frame %d   (no tracks yet: showing P2 detections)', fr));
    ud.ax.XTick=[]; ud.ax.YTick=[]; return;
end

trail = round(ud.efTrail.Value);
nLive = 0;
cmap = lines(7);
for t = 1:numel(tracks)
    fr_t = tracks(t).frames; pos = tracks(t).pos;
    j = find(fr_t==fr,1);
    if isempty(j), continue; end
    nLive = nLive + 1;
    c = cmap(mod(t-1,7)+1,:);
    j0 = max(1, j-trail);
    plot(ud.ax, pos(j0:j,1), pos(j0:j,2), '-', 'Color',c, 'LineWidth',1.2);
    plot(ud.ax, pos(j,1), pos(j,2), '.', 'Color','w', 'MarkerSize',9);
end
hold(ud.ax,'off');
title(ud.ax, sprintf('frame %d   live tracks=%d / %d total', fr, nLive, numel(tracks)));
ud.ax.XTick=[]; ud.ax.YTick=[];
end

% ===================================================================
function tracks = loadTracks(MD, iChan)
tracks = [];
try
    ip = MD.getProcessIndex('FilopodiaTrackingProcess',1,0);
    if isempty(ip), return; end
    pr = MD.processes_{ip};
    f = pr.outFilePaths_{1,iChan};
    if isempty(f)||exist(f,'file')~=2
        f = fullfile(pr.funParams_.OutputDirectory,'filoTracks.mat');
    end
    if exist(f,'file')~=2, return; end
    S = load(f,'adhesionTracks'); tracks = S.adhesionTracks;
catch
end
end

% ===================================================================
function [tipXY, baseXY] = loadDetections(MD, iChan, fr)
tipXY=zeros(0,2); baseXY=zeros(0,2);
try
    ip = MD.getProcessIndex('FilopodiaDetectionProcess',1,0);
    if isempty(ip), return; end
    pr = MD.processes_{ip};
    f = pr.outFilePaths_{1,iChan};
    if isempty(f)||exist(f,'file')~=2
        f = fullfile(pr.funParams_.OutputDirectory,'filoDetection.mat');
    end
    if exist(f,'file')~=2, return; end
    S = load(f,'adhesionInfo');
    if ~isfield(S,'adhesionInfo') || fr>numel(S.adhesionInfo) || isempty(S.adhesionInfo{fr}), return; end
    a = S.adhesionInfo{fr};
    for k=1:numel(a)
        if isfield(a(k),'dist') && a(k).dist > 0
            tipXY(end+1,:) = a(k).pos; %#ok<AGROW>
        else
            baseXY(end+1,:) = a(k).pos; %#ok<AGROW>
        end
    end
catch
end
end

% ===================================================================
function doApply(fig)
ud = fig.UserData;
p = ud.proc.funParams_;
p.ChannelIndex   = ud.ddChannel.Value;
p.MaxLinkDist    = ud.efLink.Value;
p.LinkUseBase    = ud.cbBase.Value;
p.MaxGapFrames   = ud.efGap.Value;
p.MinTrackLength = ud.efMinTrk.Value;
p.VelSmoothWin   = ud.efVel.Value;

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
    ud2.crtPackage=ud.pkg; ud2.crtProc=ud.proc; ud2.procID=ud.procID;
    ud2.MD=ud.MD; ud2.mainFig=ud.mainFig; ud2.handles_main=guidata(ud.mainFig);
    set(fig,'UserData', mergeStruct(get(fig,'UserData'),ud2));
    handles.figure1 = fig;
    processGUI_ApplyFcn(fig, [], handles, ud.proc.funParams_);
    appliedOK = true;
catch ME
    fprintf(2,'[filopodiaTrackingProcessGUI] processGUI_ApplyFcn failed: %s\n', ME.message);
end
if ~appliedOK
    try, h=guidata(ud.mainFig); packageGUI('refreshPackage_Callback', h.figure1, [], h);
    catch ME, fprintf(2,'[filopodiaTrackingProcessGUI] refresh fallback failed: %s\n', ME.message); end
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
