function varargout = filopodiaSegmentationProcessGUI(varargin)
%FILOPODIASEGMENTATIONPROCESSGUI  Settings GUI for Process 1 (Segmentation).
if numel(varargin)<3, error('Call from packageGUI only: crtProcGUI(''mainFig'',fig,procID)'); end
mainFig = varargin{2};   % handles.figure1
procID  = varargin{3};   % process index
pd(1)  = def('GaussianBlurSigma','Body blur sigma (px)','edit','Gaussian blur applied before body thresholding; larger = smoother edge');
pd(2)  = def('BodyThreshold','Body threshold method','editstr','''otsu'' | ''rosin'' | numeric value');
pd(3)  = def('BodyOpenRadius','Body open radius (px)','edit','Morphological opening to remove filopodia roots from body mask');
pd(4)  = def('BodyClosingRadius','Body close radius (px)','edit','Morphological closing to smooth body edge');
pd(5)  = def('SteerableOrder','Steerable filter order','edit','Even integer (2 or 4); 4 gives sharper ridge response');
pd(6)  = def('SigmaArray','Steerable sigma(s) (px)','edit','Scale(s) for multi-scale steerable filter; e.g. [1 2]');
pd(7)  = def('HysteresisHigh','Ridge high threshold','edit','Upper threshold for shaft mask hysteresis; empty = auto');
pd(8)  = def('HysteresisLow','Ridge low threshold','edit','Lower threshold for shaft mask hysteresis; empty = auto');
fig = filoGUI_helper('init', mainFig, procID, 'FilopodiaSegmentationProcess', pd);
if nargout>0, varargout{1}=fig; end
end
function d = def(name,label,type,tip), d=struct('name',name,'label',label,'type',type,'tooltip',tip); end
