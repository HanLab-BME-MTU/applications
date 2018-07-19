function val=psfValue(ps,x,y,z)
%integrand for psf
x=x*ps.M;
y=y*ps.M;
%z=z*ps.M;
val=quad('integ',0,1,[],[],ps,x,y,z);
val=abs(val)^2/ps.norm;
