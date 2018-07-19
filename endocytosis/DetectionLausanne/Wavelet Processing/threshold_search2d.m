%
% Automatically look for optimum 2D
%

function [TW,TS]=threshold_search2d(N,r,alphaB)

alpha_desired=alphaB;

t0=real(sqrt(-my_LambertW(-1,-alpha_desired.^2*2*pi)));
y=1./t0;

[thr v]=fmincon(@threshold_app2d,[t0 y],[],[],[],[],[t0/8 0.01],[4*t0 1.0],@threshold_app2dcon,optimset('Display','off','TolX',1e-6,'TolFun',1e-4,'MaxFunEvals',1e6,'LargeScale','off','DiffMaxChange',1e-3),N,r,alpha_desired);

TW=thr(1);
TS=thr(2);

