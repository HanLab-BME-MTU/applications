
function f = femBodyF(x,y,fs,u,varargin)

roiXi = [];
roiYi = [];

if nargin > 4
   roiXi = varargin{1};
   roiYi = varargin{2};
end

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

if ~isempty(roiXi)
   in = inpolygon(x,y,roiXi,roiYi);
   outID = find(in == 0);
   f(outID) = 0;
end

f = reshape(f,sz);
