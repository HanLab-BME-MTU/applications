function vel=brPeskinDirect(x,omega,delta)


alpha = x(1);
beta  = x(2);

vel=delta*(alpha*exp(-omega)-beta);