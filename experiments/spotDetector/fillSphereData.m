function m=fillSphereData(radius,rto);
%FILLSPHEREDATA is a helper function to create a mask with ones inside an ellipsoid and zeros everywhere else.

if nargin==1
    rto=1;
end;
radius=ceil(radius);
l=2*radius-1;
zradius=round(rto*radius);
lz=2*zradius-1;
steps=1/(1*(lz-1));
m=zeros(l,l,lz);
for v=0:steps:1
    for w=0:steps:1
        x=(radius-1)*v;
        y=(radius-1)*w;
        z=round(rto*sqrt((radius-1/rto)^2-x^2-y^2));
        x=round(x);
        y=round(y);
        if(isreal(z))        
            f=z;
            m(radius-x:radius+x,radius-y:radius+y,zradius-f:zradius+f)=1;
        end;
    end;
end;
