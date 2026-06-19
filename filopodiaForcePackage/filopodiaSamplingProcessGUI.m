function varargout = filopodiaSamplingProcessGUI(varargin)
%FILOPODIASAMPLINGSPROCESSGUI  Settings GUI for Process 5 (Force/Intensity Sampling).
ip=inputParser; ip.addOptional('mainFig',[],@ishandle);
ip.addOptional('procID',[],@isscalar);
ip.KeepUnmatched=true; ip.parse(varargin{:});
mainFig=ip.Results.mainFig; procID=ip.Results.procID;
if isempty(mainFig)||isempty(procID), error('Call from packageGUI only.'); end
pd(1) = def('ShaftSampleStep','Shaft sample step (px)','edit','Arc-length step along tip->base shaft for force/talin profiles');
pd(2) = def('SampleRadius','Sample radius (px)','edit','Local averaging radius for talin intensity');
pd(3) = def('ForceProcessName','TFM force process name','editstr','Class name of the TFMPackage process providing traction force');
fig = filoGUI_helper('init', mainFig, procID, 'FilopodiaSamplingProcess', pd);
if nargout>0, varargout{1}=fig; end
end
function d = def(name,label,type,tip), d=struct('name',name,'label',label,'type',type,'tooltip',tip); end
