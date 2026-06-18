function filopodiaTrackingProcessGUI(varargin)
%FILPODIATRACKINGPROCESSGUI  Settings GUI for Process 3 (Tracking).
ip=inputParser; ip.addOptional('mainFig',[],@ishandle);
ip.addOptional('procID',[],@isscalar);
ip.KeepUnmatched=true; ip.parse(varargin{:});
mainFig=ip.Results.mainFig; procID=ip.Results.procID;
if isempty(mainFig)||isempty(procID), error('Call from packageGUI only.'); end
pd(1) = def('MaxLinkDist','Max link distance (px/frame)','edit','Maximum distance to link detections between frames');
pd(2) = def('MaxGapFrames','Max gap frames','edit','Maximum number of frames to bridge over a gap');
pd(3) = def('MinTrackLength','Min track length (frames)','edit','Discard tracks shorter than this');
filoGUI_helper('init', mainFig, procID, 'FilopodiaTrackingProcess', pd);
end
function d = def(name,label,type,tip), d=struct('name',name,'label',label,'type',type,'tooltip',tip); end
