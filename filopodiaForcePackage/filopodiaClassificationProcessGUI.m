function varargout = filopodiaClassificationProcessGUI(varargin)
%FILPODIACLASSIFICATIONPROCESSGUI  Settings GUI for Process 4 (Classification).
ip=inputParser; ip.addOptional('mainFig',[],@ishandle);
ip.addOptional('procID',[],@isscalar);
ip.KeepUnmatched=true; ip.parse(varargin{:});
mainFig=ip.Results.mainFig; procID=ip.Results.procID;
if isempty(mainFig)||isempty(procID), error('Call from packageGUI only.'); end
pd(1)  = def('MinTipLifetime','Min tip lifetime (frames)','edit','Well-tracked tip must persist at least this long');
pd(2)  = def('MinLinearFrac','Min trajectory linearity','edit','PCA principal-axis variance fraction; 0.85 = mostly straight');
pd(3)  = def('MinTipDist','Min tip dist from body (px)','edit','Tip must reach at least this far outside the body');
pd(4)  = def('ShaftBand','Shaft adhesion band (px)','edit','Adhesion within this perp. distance of tip->base line = shaft adhesion');
pd(5)  = def('BodyMaxAngle','Body max angle (deg)','edit','Shaft must be within this angle of the body-ward direction');
pd(6)  = def('SweepRange','Direction sweep range (deg)','edit','Candidate shaft directions swept +/- this from seed');
pd(7)  = def('WShaft','Shaft adhesion reward (WShaft)','edit','Score bonus per shaft adhesion along the ray (encourages deep base)');
pd(8)  = def('WLen','Length penalty (WLen)','edit','Score penalty proportional to shaft length (prefer radial)');
pd(9)  = def('WPrior','Neighbor prior weight (WPrior)','edit','Angular deviation from already-fixed neighbors; higher = smoother direction field');
pd(10) = def('WOverlap','Shaft overlap penalty (WOverlap)','edit','Penalty if shaft crosses an already-fixed shaft');
pd(11) = def('WBaseSep','Base separation penalty (WBaseSep)','edit','Penalty if base is near an already-fixed base');
pd(12) = def('MinBaseSep','Min base separation (px)','edit','Bases closer than this get the WBaseSep penalty');
pd(13) = def('MinReachFrac','Min body-reach fraction','edit','Accept tip track only if it reaches the body in >= this fraction of frames');
fig = filoGUI_helper('init', mainFig, procID, 'FilopodiaClassificationProcess', pd);
if nargout>0, varargout{1}=fig; end
end
function d = def(name,label,type,tip), d=struct('name',name,'label',label,'type',type,'tooltip',tip); end
