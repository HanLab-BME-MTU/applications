function val=acceptMatch(s0,sImg,ndeg)

NOISELEV_PROB=[0.05 0.95];
t=(s0/sImg)^2;
pr=fcdf(t,ndeg,ndeg);

val=pr>NOISELEV_PROB(1) & pr<NOISELEV_PROB(2);
%val=1