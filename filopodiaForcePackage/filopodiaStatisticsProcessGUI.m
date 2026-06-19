function varargout = filopodiaStatisticsProcessGUI(varargin)
%FILPODIASTATISTICSPROCESSGUI  Settings GUI for Process 6 (Statistics).
if numel(varargin)<3, error('Call from packageGUI only: crtProcGUI(''mainFig'',fig,procID)'); end
mainFig = varargin{2};   % handles.figure1
procID  = varargin{3};   % process index
pd(1) = def('MinLifetimeForStats','Min lifetime for stats (frames)','edit','Ignore filopodia shorter than this when computing statistics');
pd(2) = def('ExportCSV','Export CSV','checkbox','Write per-filopodium stats to a .csv file');
pd(3) = def('MakeFigures','Make figures','checkbox','Automatically generate summary plots after running');
fig = filoGUI_helper('init', mainFig, procID, 'FilopodiaStatisticsProcess', pd);
if nargout>0, varargout{1}=fig; end
end
function d = def(name,label,type,tip), d=struct('name',name,'label',label,'type',type,'tooltip',tip); end
