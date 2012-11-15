function stack=generateSTORMmovie(pxSize,lambda,kon,koff,kb,varargin)
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
%      nPhotons -> number of emitted photons per blink, default 1000
%
%   output:
%       stack -> imSize-by-imSize-by-nCycles*nRep array
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
ip.addOptional('nPhotons',1000,@isscalar);

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
nPhotons=ip.Results.nPhotons;

if kon*tAct >= 1.0 || koff*tExp >= 1.0 || tExp*kb >= 1.0
    emsg='probability from rate constants greater than 1';
    error(emsg);
end

% how many emitters??
if nEmi == 0
    if isempty(pos)
        nEmi=floor(imSize^2*(pxSize/1000)^2);
        pos=rand(nEmi,2)*imSize;
    else
        nEmi=numel(pos(:,1));
    end
else
    if isempty(pos)
        pos=rand(nEmi,2)*imSize;
    end
end

nEmi=numel(pos(:,1));

% thermal activation rate is much smaller than the UV induced rate
kth=kon/1000;

% 3D array for movie output
stack=zeros(imSize,imSize,nCycles*nRep);

% vector for monitoring state of emitters (ON=1, OFF=0, BLEACHED=2)
emState=zeros(nEmi,1);

% average spread of localizations around true mean
sigmaLoc=(lambda/2)/sqrt(nPhotons)

% spread of point spread function: sigma=FWHM*2*sqrt[2*ln(2)]
sigmaPSF=(672/2)/2.35482

% avergae amplitude
aveAmp=nPhotons/sqrt(2*pi*(sigmaLoc/pxSize)^2);

% randomize initial state of emitters (this can of course be altered)
for ne=1:nEmi
    xi=rand();
    if xi < 0.01
        emState(ne)=1;
    end
end

% begin imaging
nFrame=1;
for nc=1:nCycles
    
    % activation pulse
    for ne=1:nEmi
        switch( emState(ne) )
            % emitter is in OFF state
            case 0
                xi=rand();
                if xi < kon*tAct
                    emState(ne)=1;
                end
        end
    end
    
    % step through imaging cycle
    for nr=1:nRep
        for ne=1:nEmi
            switch( emState(ne) )
                % emitter is in ON state
                case 1
                    % create Gaussian spot at emitter position
                    x=pos(ne,1);
                    y=pos(ne,2);
                    amp=aveAmp+sqrt(aveAmp)*randn();
                    signal=simGaussianSpot(
                    % does emitter turn off or bleach?
                    xi=rand()
                    if xi < tExp*(koff+kb)
                        zeta=rand();
                        if zeta < koff/(koff+kb)
                            emState(ne)=0;
                        else
                            emState(ne)=2;
                        end
                    end
                % thermal activation from OFF state???
                case 0
                    xi=rand();
                    if xi < tExp*kth
                        emState(ne)=1;
                    end
            end
        end
    end        
end

end

