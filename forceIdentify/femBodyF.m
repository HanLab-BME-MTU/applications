
function f = femBodyF(x,y,fs,u)

if isempty(fs)
   f = u*ones(size(x));
   return;
end

name = fs.dim;
sz = size(x);
x  = reshape(x,1,prod(sz));
y  = reshape(y,1,prod(sz));

[f,pe] = postinterp(fs,name,[x;y],'u',u);
f(pe) = 0;

f = reshape(f,sz);
