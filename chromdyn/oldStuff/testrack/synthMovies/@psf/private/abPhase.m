function abPh=abPhase(ps,rho,zv)
%phase aberration
% brho=0;
% crho=0;
% for i=1:length(ps.parm2)
%    brho=brho+ps.parm2(i)*rho.^(2*(i-1));
% end;
% for i=1:length(ps.parm3)
%    crho=crho+ps.parm3(i)*rho.^(2*(i-1));
% end;
% 
% abPh=zv*crho+brho;
no=sqrt((1-(ps.NA/ps.parm3(1)*rho).^2));
abPh=2*pi/ps.wvl *(zv*ps.parm3(1)*no + ps.parm2(1)*100*(sqrt((1-(ps.NA/ps.parm2(1)*rho).^2))-ps.parm3(1)^2/ps.parm2(1)^2*no ));