function f = spBndTF(s,sp)

if isnumeric(sp)
   f = sp*ones(size(s));
   return;
end

sz = size(s);
s  = reshape(s,1,prod(sz));

smin = min(sp.knots);
smax = max(sp.knots);

f   = zeros(size(s));
ind = find(s>smin & s<smax);

if ~isempty(ind)
   f(ind) = fnval(sp,s(ind)).';
end

f = reshape(f,sz);
