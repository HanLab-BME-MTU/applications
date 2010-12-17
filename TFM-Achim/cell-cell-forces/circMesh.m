function [p,e,t]=circMesh(curve,numRef,hmaxVal)
% curve Nx2 Matrix: first column has to be x-coordinate, second column has
%                   to be the y-coordinate.

% Run along the whole curve, and determine the curve length up to each
% point and calculate the arc length for each point, normalized to 2*pi:

numPts=length(curve);
cumLength=zeros(numPts,1);
for i=1:length(curve)
    if i>1
        cumLength(i)=cumLength(i-1)+sqrt((curve(i,1)-curve(i-1,1))^2+(curve(i,2)-curve(i-1,2))^2);
    end
end
% The arc length of each point is then given by:
% b=2*pi*cumLength/cumLength(end);

gd=zeros(2*numPts+2,1);

gd(1)=2; % That means the geometry is a polygon;
gd(2)=numPts; % That is the number of segments along the boundary.

gd(3:numPts+2)=curve(:,1); % These are the x-coordinates
gd(numPts+3:2*numPts+2)=curve(:,2); % These are the y-coordinates

dl=decsg(gd);

%figure(1)
%pdegplot(dl), axis equal

if nargin<3 || isempty(hmaxVal)
    [p,e,t]=initmesh(dl,'Hgrad',1.01);%,1.0001);
    display('Hgrad might be set to lower values ~1.0001?!');
else
    [p,e,t]=initmesh(dl,'Hmax',hmaxVal);
end

if nargin<3
    for i=1:numRef
        [p,e,t]=refinemesh(dl,p,e,t);
    end
end

% All edges defined by initmesh get there on boundary ID number which is
% given in the 5th row of the edge array 'e'. Here we set this to '1' by
% hand! This means that the boundary conditions are the same for all edges!
% (Otherwise boundary conditions had to be given for each edge in e):
e(5,:)=1;

figure(2)
pdemesh(p,e,t), axis equal


%pdegplot(defineGeom), axis equal 
%[p,e,t]=initmesh(defineGeom); 
%pdemesh(p,e,t), axis equal

end



% function [x,y]=defineGeom(bs,s)
%     nbs=4; 
% 
%     if nargin==0  
%       x=nbs;   
%       return; 
%     end
% 
%     % The following should be changed for a rectangular mesh. Then the intervall
%     % will depend on the length of each edge!
% 
%     dl=[  0      pi/2   pi       3*pi/2
%           pi/2   pi     3*pi/2   2*pi
%           1      1      1        1
%           0      0      0        0];
% 
%     if nargin==1   
%       x=dl(:,bs);   
%       return;
%     end 
% 
%     x=zeros(size(s)); 
%     y=zeros(size(s)); 
% 
%     s_new=pdearcl(b,curve,s,0,2*pi);
% 
%     x(:) = interp1(b,curve(:,1),s_new,'linear');
%     y(:) = interp1(b,curve(:,2),s_new,'linear');
% 
% end

%create the following data as an example:

% numVert=100;
% alpha=linspace(0,2*pi,numVert);
% R=1+0*rand(size(alpha));
% circ_curve_x=R.*cos(alpha);
% circ_curve_y=R.*sin(alpha);
% curve=[circ_curve_x' circ_curve_y'];
% curve=curve(1:end-1,:);
% 
% circMesh(curve,0)

