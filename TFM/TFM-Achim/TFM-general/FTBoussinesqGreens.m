function [G]=FTBoussinesqGreens(i,j,kx,ky,E)
v=0.5;

k=sqrt(kx.^2+ky.^2);
preFactor=2*(1+v)./(E.*k.^3);

if i==1 && j==1
    G=preFactor.*((1-v).*k.^2+v.*ky.^2);
elseif (i==1 && j==2) || (i==2 && j==1) 
    G=-v*preFactor.*kx.*ky;
elseif i==2 && j==2 
    G=preFactor.*((1-v).*k.^2+v.*kx.^2);
else
    'something went wrong'
end

%remove the NaN if the Greensfunction has been evaluated at zero.
nanMat=isnan(G);
if sum(nanMat(:))>0
    %display('The Boussinesq Greensfunction has been evaluated at zero')
    %display('To resolve this, the maximal value has been set')
    %G(nanMat)=max(max(G));
    
    %display('To resolve this, this value has been set to zero')
    G(nanMat)=0;
end