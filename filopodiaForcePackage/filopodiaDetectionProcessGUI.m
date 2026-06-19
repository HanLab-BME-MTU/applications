function varargout = filopodiaDetectionProcessGUI(varargin)
%FILPODIADETECTIONPROCESSGUI  Settings GUI for Process 2 (Detection).
ip=inputParser; ip.addOptional('mainFig',[],@ishandle);
ip.addOptional('procID',[],@isscalar);
ip.KeepUnmatched=true; ip.parse(varargin{:});
mainFig=ip.Results.mainFig; procID=ip.Results.procID;
if isempty(mainFig)||isempty(procID), error('Call from packageGUI only.'); end
pd(1) = def('PSFsigma','PSF / tip sigma (px)','edit','Tip blob scale; use ~2.6 for talin-GFP tips (larger than PSF)');
pd(2) = def('Alpha','Detection alpha','edit','Significance threshold for pointSourceDetection; lower = more detections');
pd(3) = def('TipMaxDistFromBody','Max tip distance from body (px)','edit','Adhesions farther than this are discarded');
pd(4) = def('BaseInsideBand','Base inside-body band (px)','edit','Include adhesions this far inside the body edge');
pd(5) = def('UseRidgeTips','Use ridge-tip augmentation','checkbox','Add skeleton endpoints as tip candidates (usually OFF)');
fig = filoGUI_helper('init', mainFig, procID, 'FilopodiaDetectionProcess', pd);
if nargout>0, varargout{1}=fig; end
end
function d = def(name,label,type,tip), d=struct('name',name,'label',label,'type',type,'tooltip',tip); end
