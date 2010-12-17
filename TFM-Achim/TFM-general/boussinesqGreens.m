%This Greensfunction is only valid for v=0.5;
function [G]=boussinesqGreens(i,j,x,y,E)
v=0.5;

r=sqrt(x.^2+y.^2);
preFactor=(1+v)./(pi.*E.*r.^3);

if i==1 && j==1 
    G=preFactor.*((1-v).*r.^2+v.*x.^2);
elseif (i==1 && j==2) || (i==2 && j==1) 
    G=v*preFactor.*x.*y;
elseif i==2 && j==2 
    G=preFactor.*((1-v).*r.^2+v.*y.^2);
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