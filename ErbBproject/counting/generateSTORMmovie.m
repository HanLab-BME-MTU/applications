function [stack,emState,settings]=generateSTORMmovie(lambda,kon,koff,kb,varargin)
%GENERATESTORMMOVIE creates an artificial STORM experiment
%   
%   required input arguments:
%       lambda -> emission wavelength in nanometer
%          kon -> rate of activation by UV light [1/sec]
%         koff -> rate of deactivation of emitter [1/sec]
%           kb -> rate of photobleaching [1/sec]
%
%   optional input arguments:
%        imSize -> image dimension in X and Y in pixels, default 256
%        pxSize -> pixel size in nanometer, default 64
%           pos -> position of emitters, can be subpixel, default random
%          nEmi -> number of emitters, default 250
%          tAct -> length of UV activation pulse, default 0.01 sec
%          tExp -> exposure time, default 0.1 sec
%       nCycles -> number of activation cycles, default 250
%          nRep -> number of frames in activation cycle, default 10
%        ngamma -> average number of photons per activation, default 6000
%          gain -> gain of the camera, default 4
%            bg -> background noise, default 100
%
%   output:
%         stack -> imSize-by-imSize-by-nCycles*nRep image array
%       emState -> matrix monitoring emitter state over time
%      settings -> paramter settings for later analysis
%
%   US 2012/11/14
%

ip=inputParser;

ip.CaseSensitive=false;
ip.StructExpand=true;

% set required arguments
ip.addRequired('lambda',@isscalar);
ip.addRequired('kon',@isscalar);
ip.addRequired('koff',@isscalar);
ip.addRequired('kb',@isscalar);

% define optional arguments
ip.addOptional('imSize',256,@isscalar);
ip.addOptional('pxSize',64.0,@isscalar);
ip.addOptional('pos',[],@isnumeric);
ip.addOptional('nEmi',0,@isscalar);
ip.addOptional('tAct',0.01,@isscalar);
ip.addOptional('tExp',0.1,@isscalar);
ip.addOptional('nCycles',250,@isscalar);
ip.addOptional('nRep',10,@isscalar);
ip.addOptional('ngamma',6000,@isscalar);
ip.addOptional('gain',4,@isscalar);
ip.addOptional('bg',100,@isscalar);

% parse arguments
ip.parse(lambda,kon,koff,kb,varargin{:});

% required arguments
lambda=ip.Results.lambda;
kon=ip.Results.kon;
koff=ip.Results.koff;
kb=ip.Results.kb;

% optional arguments
imSize=ip.Results.imSize;
pxSize=ip.Results.pxSize;
nEmi=ip.Results.nEmi;
pos=ip.Results.pos;
tAct=ip.Results.tAct;
tExp=ip.Results.tExp;
nCycles=ip.Results.nCycles;
nRep=ip.Results.nRep;
ngamma=ip.Results.ngamma;
gain=ip.Results.gain;
bg=ip.Results.bg;

if kon*tAct >= 1.0 || koff*tExp >= 1.0 || tExp*kb >= 1.0
    emsg='probability from rate constants greater than 1';
    error(emsg);
end

% sigma of PSF: sigma=FWHM*2*sqrt[2*ln(2)]
sigmaPSF=(lambda/2)/2.35482;

% patch size of Gaussian spot
sigmaPX=sigmaPSF/pxSize;
w=4*floor(sigmaPX);
wg=-w:w;
[xg,yg]=meshgrid(wg,wg);

% how many emitters??
if nEmi == 0
    if isempty(pos)
        % if nothing is specified, one emitter per square micron
        % remote from the border of the image
        nEmi=floor((imSize-2*(w+1))^2*(pxSize/1000)^2);
        pos=rand(nEmi,2)*(imSize-2*(w+1))+w+1;
    else
        nEmi=numel(pos(:,1));
    end
else
    if isempty(pos)
        pos=rand(nEmi,2)*(imSize-2*(w+1))+w+1;
    else
        if nEmi ~= numel(pos(:,1))
            nEmi=min(nEmi,numel(pos(:,1)));
            pos=pos(1:nEmi,:);
        end
    end
end

% thermal activation rate is much smaller than the UV induced rate
kth=kon/1000;

% 3D array for movie output
stack=zeros(imSize,imSize,nCycles*nRep);

% vector for monitoring state of emitters (OFF=0, ON=1, BLEACHED=2)
emState=zeros(nEmi,1);

% average amplitude from ngamma photons
aveAmp=gain*ngamma/(2*pi*(sigmaPSF/pxSize)^2);

% save settings as return variable
settings=ip.Results;
settings.pos=pos;
settings.nEmi=nEmi;
settings.kth=kth;
settings.sigmaPSF=sigmaPSF;
settings.aveAmp=aveAmp;

% randomize initial state of emitters (this can of course be altered)
for ne=1:nEmi
    eta=rand();
    if eta < 0.01
        emState(ne)=1;
    end
end

% matrix monitoring state of emitters over time, will be returned
emStateM=NaN(nEmi,nCycles*nRep);

% begin imaging
nFrame=1;
for nc=1:nCycles
    
    % activation pulse
    for ne=1:nEmi
        switch( emState(ne) )
            % emitter is in OFF state
            case 0
                eta=rand();
                if eta < kon*tAct
                    emState(ne)=1;
                end
        end
    end
    
    % step through imaging cycle
    for nr=1:nRep
        
        emStateM(:,nFrame)=emState;
        
        frame=zeros(imSize,imSize);
        
%         idp= emState == 1;
%         px=pos(idp,1);
%         py=pos(idp,2);
%         
%         np=sum(idp);
%         amp=aveAmp+sqrt(aveAmp)*randn(np,1);
%         
%         frame=simGaussianSpots(imSize,imSize,sigmaPX,'x',px,'y',py,'A',amp,'Border','truncated');
        for ne=1:nEmi
            switch( emState(ne) )
                % emitter is in ON state
                case 1
                    % create Gaussian spot at emitter position
                    xi=round(pos(ne,1));
                    yi=round(pos(ne,2));
                    x=pos(ne,1)-xi;
                    y=pos(ne,2)-yi;
                    % lower/upper bounds in x and y
                    lbx=max(xi-w,1);
                    ubx=min(imSize,xi+w);
                    lby=max(yi-w,1);
                    uby=min(imSize,yi+w);
                    
                    % xa=max(1,xi-w):min(imSize,xi+w);
                    % ya=max(1,yi-w):min(imSize,yi+w);
                    
                    % xaa=xa-xa(1)+1;
                    % yaa=ya-ya(1)+1;
                    
                    wx=(lbx:ubx)-xi;
                    wy=(lby:uby)-yi;
                    
                    [xg,yg]=meshgrid(wx,wy);
                    
                    amp=aveAmp+sqrt(aveAmp)*randn();
                    g = amp*exp(-((xg-x).^2+(yg-y).^2) / (2*sigmaPX^2));
                    
                    xa = lbx:ubx;
                    ya = lby:uby;
                    frame(ya,xa) = frame(ya,xa) + g;
                    
                    % does emitter turn off or bleach?
                    eta=rand();
                    if eta < tExp*(koff+kb)
                        zeta=rand();
                        if zeta < koff/(koff+kb)
                            emState(ne)=0;
                        else
                            emState(ne)=2;
                        end
                    end
                % thermal activation from OFF state???
                case 0
                    eta=rand();
                    if eta < tExp*kth
                        emState(ne)=1;
                    end
            end
        end
        % add background and noise
        frame=frame+bg;
        frame=poissrnd(frame);
        stack(:,:,nFrame)=frame;
        nFrame=nFrame+1;
    end        
end

emState=emStateM;
stack=uint16(stack);

end

