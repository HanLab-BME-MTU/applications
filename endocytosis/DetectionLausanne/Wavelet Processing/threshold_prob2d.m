function v=threshold_prob2d(lambda,N,r,t0,y);

if lambda<0,
  lambda=0.1;
end;

if isnan(t0),
  alpha_t0=1;
else
  alpha_t0=1-tcdf(t0,r);
end;

part1=(mygammainc(r./(2*(lambda*y).^2),r/2)-sqrt(2/r)*lambda*y.*mygammainc(r./(2*(lambda*y).^2),(r+1)/2))/gamma(r/2);

part2=...
  alpha_t0+...
  lambda/(sqrt(2*pi)*(1+t0^2/r)^(r/2))-...
  lambda*y/2*sqrt(r/(2*pi))*mybetainc(1/(1+t0^2/r),(r+1)/2,1/2);

part3=alpha_t0;

v=part1+part2+part3;
