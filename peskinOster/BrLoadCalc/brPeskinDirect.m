function vel=brPeskinDirect(x,omega,delta)
%BRPESKINDIRECT load-velocity relation from peskin for polymerization
%input:
%         x      : paramtes of the polymerization (alpha beta)
%         omega  : load 
%         delta  : monomer size divided for 2 for actin and 13 for MT
%OUTPUT
%         vel    : velocity for the load and the param.

alpha = x(1);
beta  = x(2);

vel=delta*(alpha*exp(-omega)-beta);