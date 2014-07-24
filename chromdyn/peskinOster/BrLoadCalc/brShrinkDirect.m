function vel=brShrinkDirect(x,omega,delta);
%BRSHRINKDIRECT%direct vel load relationship shrinking MT according to
%peskin & oster
%input:
%         x      : paramtes of the polymerization (beta gamma p)
%         omega  : load 
%         delta  : monomer size divided for 2 for actin and 13 for MT
%OUTPUT
%         vel    : velocity for the load and the param.


beta=x(1);
gama=x(2);
p=x(3);

velUp=delta*gama*exp(-omega/2).*(p*gama*(exp(omega/2)-exp(-omega/2))+beta);
velDown=gama*(exp(omega/2)-p*exp(-omega/2))+beta;

vel=velUp./velDown;