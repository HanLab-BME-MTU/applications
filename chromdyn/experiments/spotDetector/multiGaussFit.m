function [gaussFit, gaussgrad] = multiGaussFit(fSze,parms,dataProperties)
%MULTIGAUSSFIT try to fit two gausssian distributions to image
%
% SYNOPSIS [fit,parms,parms] = multiGaussFit(fSze,parms)
%
% INPUT fSze   : size
%       parms  : 
% OUTPUT parms     :  center coordinates 
%        gaussgrad :  goodness of fit

% c: 29/05/01	dT

%CONST DEFINITIONS
FT_SIGMA=dataProperties.FT_SIGMA;


background=parms(end);

%only cengter of gaussian
GAUSS_PARMS=4; 
%other parms
GLOB_PARMS=1;

nGs=(length(parms)-GLOB_PARMS)/GAUSS_PARMS;

if nGs ~= round(nGs)
    error('wrong number of parameters');
end;

gs_parms=[];

for i=1:nGs
    st=(i-1)*GAUSS_PARMS+1;
    center=parms(st:st+2);
    weight=parms(st+3);
    gs_parms=[gs_parms center FT_SIGMA weight];
end;

[gaussFit gaussgrad]=multiGauss(fSze,gs_parms);

%snip gradient
tmpgrad=zeros(prod(fSze),nGs*GAUSS_PARMS);

for i=1:nGs
    % *7 because there are 7 parms in multiGauss
    st=(i-1)*7+1;
    tmpgrad(:,(i-1)*GAUSS_PARMS+1:i*GAUSS_PARMS)=[gaussgrad(:,st:st+2) gaussgrad(:,st+6)];
end;

% gradient of bg
sca=ones(length(gaussFit(:)),1);
gaussgrad=[tmpgrad  sca];

% add background
gaussFit=gaussFit+background;
