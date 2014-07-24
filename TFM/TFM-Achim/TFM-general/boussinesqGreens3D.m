%This Greensfunction is only valid for v=0.5;
function [G]=boussinesqGreens3D(i,j,x,y,z,E,v)
% This function returns the i,j component of the Greensfunction for the
% most general 3D problem of a force acting on the top surface. z gives the
% distance from this surface.

if nargin < 7 || isempty(v)
    v=0.5;
end

if isempty(z)
    z=0;
end

r=sqrt(x.^2+y.^2+z.^2);
preFactor=(1+v)./(2*pi.*E.*r.^3);
a = (2*(1-v)*r+z)./(r+z);
b = (2*r.*(v*r+z)+z.^2)./(r+z).^2;
c = (1-2.*v)./(r+z);
if i==1 && j==1 
    G=preFactor.*(a.*r.^2 + b.*x.^2);
elseif (i==1 && j==2) || (i==2 && j==1) 
    G=preFactor.*(b.*x.*y);
elseif i==1 && j==3 
    G=preFactor.*(x.*z-c*r.^2.*x);
elseif i==2 && j==2
    G=preFactor.*(a.*r.^2 + b.*y.^2);
elseif i==2 && j==3
    G=preFactor.*(y.*z-c*r.^2.*y);
elseif i==3 && j==1
    G=preFactor.*(x.*z+c*r.^2.*x);
elseif i==3 && j==2
    G=preFactor.*(y.*z+c*r.^2.*y);
elseif i==3 && j==3
    G=preFactor.*(z.^2+2*(1-v)*r.^2);
else
    'something went wrong'
end

%remove the NaN if the Greensfunction has been evaluated at zero.
nanMat=isnan(G);
if sum(sum(nanMat))>0
    %display('The Boussinesq Greensfunction has been evaluated at zero')
    %display('To resolve this, the maximal value has been set')
    %G(nanMat)=max(max(G));
    
    %display('To resolve this, this value has been set to zero')
    G(nanMat)=0;
end