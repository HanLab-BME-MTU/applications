function fy = spMyoDFy(x,y,fs,u)

if isempty(fs)
   fy = u*ones(size(x));
   return;
end

name = fs.dim;
sz = size(x);
x  = reshape(x,1,prod(sz));
y  = reshape(y,1,prod(sz));

[fy,pe] = postinterp(fs,name{2},[x;y],'u',[zeros(size(u));u]);
fy(pe) = 0;

fy = reshape(fy,sz);
