function [gaussFit, gaussgrad] = multiGauss(fSze,parms)
%MULTIGAUSS create superpos of gausssian distributions
%
% SYNOPSIS [gaussFit, gaussgrad] = multiGauss(fSze,parms)
%
% INPUT fSze   : size
%       parms  : 
% OUTPUT gaussFit  :  image of gaussians
%        gaussgrad :  image gradient

% c: 29/05/01	dT

GAUSS_PARMS=7;

nGs=length(parms)/GAUSS_PARMS;

if nGs ~= round(nGs)
    error('wrong number of parameters');
end;


hSze=floor(fSze/2);
gaussFit=zeros(fSze);
gaussgrad=zeros(prod(fSze),nGs*GAUSS_PARMS);

for i=1:nGs
    st=(i-1)*GAUSS_PARMS+1;
    center = parms(st:st+2);
    sigma = parms(st+3:st+5);
    
    weight = parms(st+6);
    
    gm = GaussMask3D(sigma,fSze,center,0);
    
    %place the gaussians
    gaussFit=gaussFit+weight*gm;
    
    if nargout > 1

        % compute parm gradient
        x1=([-floor(fSze(1)/2):floor(fSze(1)/2)]-center(1))./sigma(1);
        y1=([-floor(fSze(2)/2):floor(fSze(2)/2)]-center(2))./sigma(2);
        z1=([-floor(fSze(3)/2):floor(fSze(3)/2)]-center(3))./sigma(3);

        sca=ones(length(gaussFit(:)),1);
        tg=[weight*gm(:)*ones(1,7)];

        [yc, xc,zc]=meshgrid(y1,x1,z1);

        gaussgrad(:,(i-1)*GAUSS_PARMS+1:i*GAUSS_PARMS)=...
            [xc(:)/sigma(1), yc(:)/sigma(2), zc(:)/sigma(3),...
            xc(:).^2/sigma(1) +yc(:).^2/sigma(2), zc(:).^2/sigma(3), sca/weight ].*tg;

    end
end

