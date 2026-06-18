function filopodiaStatisticsProcessGUI(varargin)
%FILPODIASTATISTICSPROCESSGUI  Settings GUI for Process 6 (Statistics).
ip=inputParser; ip.addOptional('mainFig',[],@ishandle);
ip.addOptional('procID',[],@isscalar);
ip.KeepUnmatched=true; ip.parse(varargin{:});
mainFig=ip.Results.mainFig; procID=ip.Results.procID;
if isempty(mainFig)||isempty(procID), error('Call from packageGUI only.'); end
pd(1) = def('MinLifetimeForStats','Min lifetime for stats (frames)','edit','Ignore filopodia shorter than this when computing statistics');
pd(2) = def('ExportCSV','Export CSV','checkbox','Write per-filopodium stats to a .csv file');
pd(3) = def('MakeFigures','Make figures','checkbox','Automatically generate summary plots after running');
filoGUI_helper('init', mainFig, procID, 'FilopodiaStatisticsProcess', pd);
end
function d = def(name,label,type,tip), d=struct('name',name,'label',label,'type',type,'tooltip',tip); end
