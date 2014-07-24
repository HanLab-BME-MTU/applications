function features = analyzeSTORMmovie(stack,settings,varargin)
%ANALYZESTORMMOVIE analyzes an artifical STORM movie from generateSTORMmovie 
%
%   required input arguments:
%          stack -> stack of images
%       settings -> struct containing settings of generateSTORMmovie
%
%   optional input arguments:
%         alpha -> alpha value for maxima detection, default 0.05
%         doMMF -> perform mixture model fit, default false
%       display -> show results, default false
%
%   output:
%       features -> feature information as output from pointSourceDetection
%
%   US 2012/11/15
%

%warning off;

ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=false;

ip.addRequired('stack',@isnumeric);
ip.addRequired('settings',@isstruct);

ip.addOptional('display',false,@islogical);
ip.addOptional('alpha',0.05,@isscalar);
ip.addOptional('doMMF',false,@islogical);

ip.parse(stack,settings,varargin{:});

stack=ip.Results.stack;
settings=ip.Results.settings;
alpha=ip.Results.alpha;
display=ip.Results.display;
doMMF=ip.Results.doMMF;


% extract variables from settings, not all are necessary
%lambda=settings.lambda;
%kon=settings.kon;
%koff=settings.koff;
%kb=settings.kb;
imSize=settings.imSize;
pxSize=settings.pxSize;
%nEmi=settings.nEmi;
pos=settings.pos;
%tAct=settings.tAct;
%tExp=settings.tExp;
nCycles=settings.nCycles;
nRep=settings.nRep;
%ngamma=settings.ngamma;
%gain=settings.gain;
%bg=settings.bg;
%kth=settings.kth;
sigmaPSF=settings.sigmaPSF;
%aveAmp=settings.aveAmp;

nFrames=min(size(stack,3),nCycles*nRep);
features=cell(nFrames,1);

for k=1:nFrames
    
    img=double(stack(:,:,k));
    sigma=sigmaPSF/pxSize;
    features{k}=pointSourceDetection(img,sigma,'alpha',alpha,'FitMixtures',doMMF);
    
    progressText(k/nFrames,'All work and no play makes Jack a dull boy');
    
end

fprintf(1,'\n');

% display results and compare with true values
if display
    
    % p: estimated localizations
    p=extractF(features);
    p=p(:,1:2);
    % pos: true localizations
    plot(pos(:,1),pos(:,2),'ro',p(:,1),p(:,2),'kx');
    set(gca,'DataAspectRatio',[1 1 1]);
    set(gca,'XLim',[0 imSize],'YLim',[0 imSize]);
end

end

