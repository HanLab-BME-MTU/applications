
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
   if ~iscell(roiXi)
      in = inpolygon(x,y,roiXi,roiYi);
      iOut = find(in ~= 1);
   else
      iOut = 1:length(x);
      for kk = 1:length(roiXi)
         in  = inpolygon(x,y,roiXi{kk},roiYi{kk});
         iIn = find(in==1);
         iOut(iIn) = 0;
      end
      iOut(find(iOut==0)) = [];
   end
   f(iOut) = 0;
end

f = reshape(f,sz);
