function varargout = filopodiaSamplingProcessGUI(varargin)
%FILOPODIASAMPLINGSPROCESSGUI  Settings GUI for Process 5 (Force/Intensity Sampling).
if numel(varargin)<3, error('Call from packageGUI only: crtProcGUI(''mainFig'',fig,procID)'); end
mainFig = varargin{2};   % handles.figure1
procID  = varargin{3};   % process index
pd(1) = def('ShaftSampleStep','Shaft sample step (px)','edit','Arc-length step along tip->base shaft for force/talin profiles');
pd(2) = def('SampleRadius','Sample radius (px)','edit','Local averaging radius for talin intensity');
pd(3) = def('ForceProcessName','TFM force process name','editstr','Class name of the TFMPackage process providing traction force');
fig = filoGUI_helper('init', mainFig, procID, 'FilopodiaSamplingProcess', pd);
if nargout>0, varargout{1}=fig; end
end
function d = def(name,label,type,tip), d=struct('name',name,'label',label,'type',type,'tooltip',tip); end
