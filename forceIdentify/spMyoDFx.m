function fx = spMyoDFx(x,y,fs,u)

if isempty(fs)
   fx = u*ones(size(x));
   return;
end

name = fs.dim;
sz = size(x);
x  = reshape(x,1,prod(sz));
y  = reshape(y,1,prod(sz));

[fx,pe] = postinterp(fs,name{1},[x;y],'u',[u;zeros(size(u))]);
fx(pe) = 0;

fx = reshape(fx,sz);
