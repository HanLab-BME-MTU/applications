function vel=brShrinkDirect(x,omega,delta);
% direct vel load relationship shrinking MT according to peskin & oster


beta=x(1);
gama=x(2);
p=x(3);

velUp=delta*gama*exp(-omega/2).*(p*gama*(exp(omega/2)-exp(-omega/2))+beta);
velDown=gama*(exp(omega/2)-p*exp(-omega/2))+beta;

vel=velUp./velDown;