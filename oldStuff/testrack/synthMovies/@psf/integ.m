function intP=integ(rho,ps,x,y,z)
%integrand for psf
arho=0;
%  for i=1:length(ps.parm1)
% arho=arho+ps.parm1(i)*rho.^(i-1);
% end;

%----- temporary 14/6/2002:
arho=rho.^0;
%--------------

a=ps.zd*ps.NA/sqrt(ps.M^2-ps.NA^2);
intP=besselj(0,2*pi/ps.wvl*a*rho*sqrt(x^2+y^2)/ps.zd).*arho.*exp(j*abPhase(ps,rho,z)).*rho;
