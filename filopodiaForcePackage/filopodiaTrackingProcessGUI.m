function varargout = filopodiaTrackingProcessGUI(varargin)
%FILPODIATRACKINGPROCESSGUI  Settings GUI for Process 3 (Tracking).
if numel(varargin)<3, error('Call from packageGUI only: crtProcGUI(''mainFig'',fig,procID)'); end
mainFig = varargin{2};   % handles.figure1
procID  = varargin{3};   % process index
pd(1) = def('MaxLinkDist','Max link distance (px/frame)','edit','Maximum distance to link detections between frames');
pd(2) = def('MaxGapFrames','Max gap frames','edit','Maximum number of frames to bridge over a gap');
pd(3) = def('MinTrackLength','Min track length (frames)','edit','Discard tracks shorter than this');
fig = filoGUI_helper('init', mainFig, procID, 'FilopodiaTrackingProcess', pd);
if nargout>0, varargout{1}=fig; end
end
function d = def(name,label,type,tip), d=struct('name',name,'label',label,'type',type,'tooltip',tip); end
