function msh = rectmesh(curvL,curvT,curvR,curvB,hIn,vIn)
%RECTMESH: This function provides anisothtropic mesh for deformed rectangle
%          whose 4 sides can be any curves.
%
% The 4 side curves are specified by 'curvL' (Left), 'curvB' (Bottom), 'curvT'
% (Top) and 'curvR' (Right).
%
% SYNOPSIS :
%    msh = rectMesh(curvL,curvB,curvT,curvR,hIn,vIn);
%
% INPUT :
%    curvL : 2-by-m numerical matrix where m>=2 is the number of data points
%       on the Left-side curve. 'curvL(1,:)' are the x-coordinates and
%       'curvL(2,:)' are the y-coordinates. 
%       'curvB', 'curvT' and 'curvR' are defined in the same way. The number
%       of data points on the Left and Right sides (or the Bottom and Top
%       sides) have to be the same. But between the Left-Right and the
%       Top-Bottom pairs, they can be different.
%       The data points on the 4 side curves have to be given in the 
%       clockwise order.
%    vIn : It is either the number of meshing points in the vertical direction 
%      (on the Left and Right sides) or a vector of parameter values (between 
%      0 and 1) whose length is > 2. The parameter values specify how the left
%      and right sides are divided unevenly. It can also be an m-by-2 matrix 
%      of parameter values (between 0 and 1) whose first column is for
%      dividing the left sides and whose right column is for dividing the
%      right sides.
%    hIn : The same as 'vIn' but for the horizontal direction (on the
%       Bottom and Top sides). If it is m-by-2 matrix, the first column is for
%       Bottom and the second column is for Top.
%
% OUTPUT :
%    msh : The FEMLAB mesh structure. It can be used to define the geometry
%       using the FEMLAB concept, mesh as geometry.

%Label the 4 side curves according to the edge labels of the femlab 
% rectangular mesh.
curve{1} = curvL;
curve{2} = curvT;
curve{3} = curvR;
curve{4} = curvB;

%Compute the curve parameters for the data points on each curve using the
% chord length between each data points. The curve parameter runs from 0 to 1.
% The 4 curves will be described by spline functions given in piecewise
% polynomial form: 'ppx' and 'ppy'

ss = cell(4,1); %Curve parameters.
cL = ones(4,1); %Curve length.
for j = 1:4
   ss{j} = [0 sqrt(diff(curve{j}(1,:)).^2+diff(curve{j}(2,:)).^2)];
   for k = 3:length(ss{j})
      ss{j}(k) = ss{j}(k-1)+ss{j}(k);
   end
   cL(j) = ss{j}(end);

   %Normalize the curve parameter so that it runs from 0 to 1.
   ss{j} = ss{j}/cL(j);
end

%Our task is to establish a one to one map between a square of side length 
% 'hScale-by-vScale' to the deformed rectangle bounded by the 4 input curves. 
% We shall think of the square
% as a piece of elastic membrane that can be stretched so that the 4 sides of
% the square match the 4 boundary curves of the deformed rectangle. This
% physical process can be numerically simulated by solving the elastic
% equation for the membrane.

%Set the side length of the rectangle.
hScale = max(cL(2),cL(4));
vScale = max(cL(1),cL(3));

%The four sides of the rectangle.
rect{1} = [zeros(1,size(curvL,2));ss{1}*vScale];
rect{2} = [ss{2}*hScale;vScale*ones(1,size(curvT,2));];
rect{3} = [hScale*ones(1,size(curvR,2));(1-ss{3})*vScale];
rect{4} = [(1-ss{4})*hScale;zeros(1,size(curvB,2));];

%The displacements of the four sides stretched to match the four curves.
dispx = cell(4,1);
dispy = cell(4,1);
for j = 1:4
   dispx{j} = curve{j}(1,:)-rect{j}(1,:);
   dispy{j} = curve{j}(2,:)-rect{j}(2,:);

   %Interpolating the displacements.
   ppx{j} = spline(ss{j},dispx{j});
   ppy{j} = spline(ss{j},dispy{j});
end

ind = 1;
bndInd = [1 2 3 4];

options = elOptionsSet('EPType','YModulPRatio','BCType', ...
   {'Dirichlet','Dirichlet','Dirichlet','Dirichlet'});

fn.YModul = 1;
fn.PRatio = 0.3;

%fn.VDragCoef = 0;
%fn.TimeStep  = 1;

msg = ['Something is wrong with the last two arguments, ' ...
   '''hIn'' and ''vIn''.'];
[m,n] = size(hIn);
if m == 1 & n == 1
   dH = hScale/(hIn-1); 
   smplH1 = linspace(0,hScale,hIn);
   smplH2 = linspace(dH/2,hScale-dH/2,hIn-1);
elseif min(m,n) == 1
   if max(m,n) == 2
      error(msg);
   end

   smplH1 = hIn*hScale;
   smplH2 = (hIn(1:end-1)+hIn(2:end))/2*hScale;
else
   error(msg);
end

[m,n] = size(vIn);
if m == 1 & n == 1
   dV = vScale/(vIn-1);
   smplV1 = linspace(0,vScale,vIn);
   smplV2 = linspace(dV/2,vScale-dV/2,vIn-1);
elseif min(m,n) == 1
   if max(m,n) == 2
      error(msg);
   end

   smplV1 = vIn*vScale;
   smplV2 = (vIn(1:end-1)+vIn(2:end))/2*vScale;
else
   error(msg);
end

[gridH1,gridV1] = meshgrid(smplH1,smplV1);
[gridH2,gridV2] = meshgrid(smplH2,smplV2);

msh.p = [reshape(gridH1,1,prod(size(gridH1))) ...
   reshape(gridH2,1,prod(size(gridH2))); ...
   reshape(gridV1,1,prod(size(gridV1))) ...
   reshape(gridV2,1,prod(size(gridV2)))];

msh = meshenrich(msh);

%arcLen = zeros(4,1);
%We need to map each edge element on the four sides of the rectangular mesh to 
% a straight line on the four curves (now represented by 'ppx' and 'ppy').
% Doing it this way will keep the interior nodes of the rectangular mesh,
% 'msh' interior after it is mapped to the curved rectangular.
for k = 1:4
   bndEI = find(msh.e(5,:)==k);
   edgeS1 = [0 msh.e(4,bndEI)];
   [edgeS1,sortI] = sort(edgeS1);
   arcLen = edgeS1(end);

   edgeUx1 = ppval(ppx{k},edgeS1/arcLen);
   edgeUy1 = ppval(ppy{k},edgeS1/arcLen);

   edgeS2  = (edgeS1(1:end-1)+edgeS1(2:end))/2;
   edgeUx2 = (edgeUx1(1:end-1)+edgeUx1(2:end))/2;
   edgeUy2 = (edgeUy1(1:end-1)+edgeUy1(2:end))/2;

   edgeS  = zeros(1,length(edgeS1)+length(edgeS2));
   edgeUx = zeros(1,length(edgeUx1)+length(edgeUx2));
   edgeUy = zeros(1,length(edgeUy1)+length(edgeUy2));

   edgeS(1:2:end)    = edgeS1;
   edgeS(2:2:end-1)  = edgeS2;
   edgeUx(1:2:end)   = edgeUx1;
   edgeUx(2:2:end-1) = edgeUx2;
   edgeUy(1:2:end)   = edgeUy1;
   edgeUy(2:2:end-1) = edgeUy2;

   ppx{k} = spline(edgeS/arcLen,edgeUx);
   ppy{k} = spline(edgeS/arcLen,edgeUy);
end

[gridH3,gridV3] = meshgrid(smplH2,smplV1);
[gridH4,gridV4] = meshgrid(smplH1,smplV2);

msh2.p = [msh.p [reshape(gridH3,1,prod(size(gridH3))) ...
   reshape(gridH4,1,prod(size(gridH4))); ...
   reshape(gridV3,1,prod(size(gridV3))) ...
   reshape(gridV4,1,prod(size(gridV4)))]];

msh2 = meshenrich(msh2);

for j = 4:-1:1
   k = 5-j;
   bndEI = find(msh2.e(5,:)==k);
   edgeS1 = [0 msh2.e(4,bndEI)];
   [edgeS1,sortI] = sort(edgeS1);
   arcLen = edgeS1(end);

   edgeUx1 = ppval(ppx{j},1-edgeS1/arcLen);
   edgeUy1 = ppval(ppy{j},1-edgeS1/arcLen);

   edgeS2  = (edgeS1(1:end-1)+edgeS1(2:end))/2;
   edgeUx2 = (edgeUx1(1:end-1)+edgeUx1(2:end))/2;
   edgeUy2 = (edgeUy1(1:end-1)+edgeUy1(2:end))/2;

   edgeS  = zeros(1,length(edgeS1)+length(edgeS2));
   edgeUx = zeros(1,length(edgeUx1)+length(edgeUx2));
   edgeUy = zeros(1,length(edgeUy1)+length(edgeUy2));

   edgeS(1:2:end)    = edgeS1;
   edgeS(2:2:end-1)  = edgeS2;
   edgeUx(1:2:end)   = edgeUx1;
   edgeUx(2:2:end-1) = edgeUx2;
   edgeUy(1:2:end)   = edgeUy1;
   edgeUy(2:2:end-1) = edgeUy2;

   pp2x{k} = spline(edgeS,edgeUx);
   pp2y{k} = spline(edgeS,edgeUy);

   fn.BndDispx{k} = 'bndDisp';
   fn.BndDispy{k} = 'bndDisp';
   fp.BndDispx{k} = {{'s'} {pp2x{k}}};
   fp.BndDispy{k} = {{'s'} {pp2y{k}}};
end

fem = elModelAssemble([],msh2,options,fn,fp,ind,bndInd);
fem = elasticSolve(fem,[]);
[dispU1,dispU2] = postinterp(fem,'u1','u2',msh.p);

msh.p = msh.p+[dispU1;dispU2];
%msh.e = mesh.e;
%msh.t = mesh.t;
%msh.v = mesh.v;

%msh   = meshenrich(msh);


