function ps=psf(NA,M,la,zd,a,b,c);
% generate aberrated PSF object
% SYNOPSIS : psf=abpsf(NA,M,la,a,b,c)
%
% INPUT:   NA = objective numerical aperture
%          M = objective magnification
%          la = fluorescent wavelenght
%          zd = distance obj-detector
% 
%
% OUTPUT: Psf object 
%
% constructor:
%                  abpsf(tImg,sImg)
%
% public methods:
%                 psdata=makePSF(ps,dim)
%                 val=psfValue(x,y,z)
%                  
% private methods:
%                   psf = integ(psf)

% c: 13/09/00 dT


%init psf data
ps.NA=NA;
ps.M=M;
ps.zd=zd;
ps.wvl=la;
ps.parm1=a;
ps.parm2=b;
ps.parm3=c;
ps.norm=0;

% default values


ps = class(ps,'psf');
ps.norm=abs(quad('integ',0,1,[],[],ps,0,0,0))^2;