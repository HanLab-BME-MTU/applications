function v=threshold_app2d(thr,N,r,alpha_desired);

t0=thr(1);
y=thr(2);

if t0<0,
    v=Inf;
elseif y>t0,
    v=t0+y;
else
    [lambda,v]=fminsearch(@threshold_prob2d,10,optimset('MaxFunEvals',500),N,r,t0,y);
    v=t0+y;
end;
