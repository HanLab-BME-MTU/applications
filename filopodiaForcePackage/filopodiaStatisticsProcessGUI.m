function varargout = filopodiaStatisticsProcessGUI(varargin)
%FILOPODIASTATISTICSPROCESSGUI  Settings window for Process 6 (statistics) of
%the FilopodiaForcePackage. Channel + a few options, plus a preview that
%summarizes the SAVED stats (filopodia count, pooled length / lifetime /
%tip force distributions) so you can sanity-check the result of the last run.
%
% Called by packageGUI:  filopodiaStatisticsProcessGUI('mainFig',fig,procID)
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

W = 1000; H = 600;
fig = uifigure('Name',[proc.getName() ' Settings'],'Position',[120 120 W H],'Resize','on');
g = uigridlayout(fig,[1 2]); g.ColumnWidth = {340,'1x'}; g.RowHeight = {'1x'};

left = uipanel(g,'Title','Settings');
sp = uigridlayout(left,[12 2]);
sp.ColumnWidth = {180,'1x'}; sp.RowHeight = repmat({26},1,12);
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

r=r+1; lbl(sp,r,'Min lifetime for stats');
ud.efLife = uieditfield(sp,'numeric','Value',gf(fp,'MinLifetimeForStats',3),'Limits',[1 Inf], ...
    'Tooltip','Ignore filopodia shorter than this many frames');
ud.efLife.Layout.Row=r; ud.efLife.Layout.Column=2;

r=r+1; lbl(sp,r,'Export CSV');
ud.cbCSV = uicheckbox(sp,'Text','','Value',logical(gf(fp,'ExportCSV',true)));
ud.cbCSV.Layout.Row=r; ud.cbCSV.Layout.Column=2;

r=r+1; lbl(sp,r,'Make figures');
ud.cbFig = uicheckbox(sp,'Text','','Value',logical(gf(fp,'MakeFigures',true)));
ud.cbFig.Layout.Row=r; ud.cbFig.Layout.Column=2;

r=r+1; lbl(sp,r,'Preview metric');
ud.ddMetric = uidropdown(sp,'Items',{'length','lifetime','tip force','tip talin'},'Value','length');
ud.ddMetric.Layout.Row=r; ud.ddMetric.Layout.Column=2;

r=r+1;
ud.btnPreview = uibutton(sp,'Text','Update preview','ButtonPushedFcn',@(~,~)doPreview(fig));
ud.btnPreview.Layout.Row=r; ud.btnPreview.Layout.Column=[1 2];

r=r+1;
ud.lblSummary = uilabel(sp,'Text','','WordWrap','on','VerticalAlignment','top');
ud.lblSummary.Layout.Row=[r r+2]; ud.lblSummary.Layout.Column=[1 2];

r=r+3;
ud.btnDone = uibutton(sp,'Text','Done (save to this movie)','ButtonPushedFcn',@(~,~)doApply(fig));
ud.btnDone.Layout.Row=r; ud.btnDone.Layout.Column=1;
ud.btnCancel = uibutton(sp,'Text','Cancel','ButtonPushedFcn',@(~,~)delete(fig));
ud.btnCancel.Layout.Row=r; ud.btnCancel.Layout.Column=2;

right = uipanel(g,'Title','Stats preview (saved result)');
rg = uigridlayout(right,[1 1]); rg.Padding=[6 6 6 6];
ud.ax = uiaxes(rg);
title(ud.ax,'Click "Update preview" (needs a saved P6 run)');

fig.UserData = ud;
ud.ddChannel.ValueChangedFcn = @(s,~) doPreview(fig);
ud.ddMetric.ValueChangedFcn = @(s,~) doPreview(fig);
fig.UserData = ud;

doPreview(fig);
if nargout>0, varargout{1}=fig; end
end

% ===================================================================
function doPreview(fig)
ud = fig.UserData;
iChan = ud.ddChannel.Value;
[st, ok] = loadStats(ud.MD, iChan);
cla(ud.ax,'reset');
if ~ok || isempty(st)
    title(ud.ax,'No saved P6 result yet: set options, Done, then run P6.');
    ud.lblSummary.Text = ''; return;
end

switch ud.ddMetric.Value
    case 'length',    v = field2(st.pooled,'Lmean_nm');   nm='length (nm)';
    case 'lifetime',  v = field2(st.pooled,'lifetime_s'); nm='lifetime (s)';
    case 'tip force', v = field2(st.pooled,'tipForceMean');nm='tip traction (Pa)';
    case 'tip talin', v = field2(st.pooled,'tipTalinMean');nm='tip talin (a.u.)';
end
v = v(isfinite(v));
if isempty(v), title(ud.ax,'No data for this metric'); return; end
histogram(ud.ax, v, max(10,round(sqrt(numel(v)))));
xlabel(ud.ax,nm); ylabel(ud.ax,'# filopodia'); grid(ud.ax,'on');
title(ud.ax, sprintf('%s   (n=%d filopodia)', nm, numel(v)));

ud.lblSummary.Text = sprintf(['Filopodia: %d\nMedian length: %.0f nm\n' ...
    'Median lifetime: %.1f s\nMedian tip force: %.0f Pa'], ...
    st.nFilopodia, nanmedian(field2(st.pooled,'Lmean_nm')), ...
    nanmedian(field2(st.pooled,'lifetime_s')), nanmedian(field2(st.pooled,'tipForceMean')));
end

% ===================================================================
function v = field2(s,n), if isfield(s,n), v=s.(n)(:); else, v=[]; end, end
function m = nanmedian(x), x=x(isfinite(x)); if isempty(x), m=NaN; else, m=median(x); end, end

% ===================================================================
function [st, ok] = loadStats(MD, iChan)
st=[]; ok=false;
try
    ip = MD.getProcessIndex('FilopodiaStatisticsProcess',1,0);
    if isempty(ip), return; end
    pr = MD.processes_{ip};
    f = pr.outFilePaths_{1,iChan};
    if isempty(f)||exist(f,'file')~=2, f = fullfile(pr.funParams_.OutputDirectory,'filoStats.mat'); end
    if exist(f,'file')~=2, return; end
    S = load(f,'stats'); st = S.stats; ok = true;
catch
end
end

% ===================================================================
function p = collectParams(ud)
p = ud.proc.funParams_;
p.ChannelIndex        = ud.ddChannel.Value;
p.MinLifetimeForStats = ud.efLife.Value;
p.ExportCSV           = ud.cbCSV.Value;
p.MakeFigures         = ud.cbFig.Value;
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
    fprintf(2,'[filopodiaStatisticsProcessGUI] processGUI_ApplyFcn failed: %s\n', ME.message);
end
if ~appliedOK
    try, h=guidata(ud.mainFig); packageGUI('refreshPackage_Callback', h.figure1, [], h);
    catch ME, fprintf(2,'[filopodiaStatisticsProcessGUI] refresh fallback failed: %s\n', ME.message); end
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
