function [c,ceq]=threshold_app2dcon(thr,N,r,alpha_desired);

t0=thr(1);
y=thr(2);

[lambda,v]=fminsearch(@threshold_prob2d,1,[],N,r,t0,y);

c=abs(v-alpha_desired)/alpha_desired;
c=c*1e1;
ceq=[];