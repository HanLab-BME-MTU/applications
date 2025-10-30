%This Greensfunction is only valid for v=0.5 if not specified;
function [G]=boussinesqGreens3D(i,j,x,y,z)
E=8000;
v=0.5;

r=sqrt(x.^2+y.^2+z.^2);
preFactor=2*(1+v)./(pi.*E);

if i==1 && j==1 
    G=preFactor.*((2.*(1-v)*(r-z)./r.*(r-z))+(2.*r.*(v.*r-z)+z.^2).*x.^2./r.^3.*(r-z).^2);
elseif (i==1 && j==2) || (i==2 && j==1)
    G=preFactor.*((2.*r.*(v.*r-z)+z.^2).*x.*y./((r.^3).*(r-z).^2));
elseif i==1 && j==3
    G=preFactor.*(x.*z./r.^3-((1-2.*v).*x)./r.*(r-z));
elseif i==2 && j==2
    G=preFactor.*((2.*(1-v)*(r-z)./r.*(r-z))+(2.*r.*(v.*r-z)+z.^2).*y.^2./r.^3.*(r-z).^2);
elseif i==2 && j==3
    G=preFactor.*(y.*z./r.^3-((1-2.*v).*y)./r.*(r-z));
elseif i==3 && j==1
    G=preFactor.*(x.*z./r.^3+((1-2.*v).*x)./r.*(r-z));
elseif i==3 && j==2
    G=preFactor.*(y.*z./r.^3-((1-2.*v).*y)./r.*(r-z));
elseif i==3 && j==3
    G=preFactor.*((2.*(1-v)./r)+(z.^2./r.^3));
else
    dis('something went wrong')
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

    
