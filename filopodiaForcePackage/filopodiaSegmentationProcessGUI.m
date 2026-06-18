function filopodiaSegmentationProcessGUI(varargin)
%FILOPODIASEGMENTATIONPROCESSGUI  Settings GUI for Process 1 (Segmentation).
ip=inputParser; ip.addOptional('mainFig',[],@ishandle);
ip.addOptional('procID',[],@isscalar);
ip.KeepUnmatched=true; ip.parse(varargin{:});
mainFig=ip.Results.mainFig; procID=ip.Results.procID;
if isempty(mainFig)||isempty(procID), error('Call from packageGUI only.'); end
pd(1)  = def('GaussianBlurSigma','Body blur sigma (px)','edit','Gaussian blur applied before body thresholding; larger = smoother edge');
pd(2)  = def('BodyThreshold','Body threshold method','editstr','''otsu'' | ''rosin'' | numeric value');
pd(3)  = def('BodyOpenRadius','Body open radius (px)','edit','Morphological opening to remove filopodia roots from body mask');
pd(4)  = def('BodyClosingRadius','Body close radius (px)','edit','Morphological closing to smooth body edge');
pd(5)  = def('SteerableOrder','Steerable filter order','edit','Even integer (2 or 4); 4 gives sharper ridge response');
pd(6)  = def('SigmaArray','Steerable sigma(s) (px)','edit','Scale(s) for multi-scale steerable filter; e.g. [1 2]');
pd(7)  = def('HysteresisHigh','Ridge high threshold','edit','Upper threshold for shaft mask hysteresis; empty = auto');
pd(8)  = def('HysteresisLow','Ridge low threshold','edit','Lower threshold for shaft mask hysteresis; empty = auto');
filoGUI_helper('init', mainFig, procID, 'FilopodiaSegmentationProcess', pd);
end
function d = def(name,label,type,tip), d=struct('name',name,'label',label,'type',type,'tooltip',tip); end
