function H=arrow3x(p1,p2,spCol,w,h,ip)
%ARROW3 cretes a 3D arrow
%   ARROW3(P1,P2) will draw vector lines (2D/3D) from P1 to P2 with
%   arrowheads, where P1 and P2 can be either nx2 matrices (for 2D),
%   or nx3 matrices (for 3D).
%
%   ARROW3(P1,P2,S,W,H,IP) can be used to specify properties of the
%   line and arrowhead; W=width of arrowhead, H=height of arrowhead.
%   If a value is provided for IP, an initial point marker will be
%   plotted with width IP; use IP = 0 for default width.  
%
%   Copyright(c)2002, Version 3.00
%       Jeff Chang <pchang@mtu.edu>
%       Tom Davis  <tdavis@eng.usf.edu>

%===============================================================
% Error Checking
if ishold, restore=1; else restore=0; end
if nargin<2
    error([upper(mfilename),' requires at least two input arguments'])
end
[r1,c1]=size(p1); [r2,c2]=size(p2);
if r1~=r2, error('P1 and P2 must have same number of rows'), end
if c1~=c2, error('P1 and P2 must have same number of columns'), end
if c1==1 | (nargin>5 & length(ip)>1)
    error(['Invalid input, type HELP ',upper(mfilename),' for usage examples'])
end
if c1==2, p1=[p1,zeros(r1,1)]; p2=[p2,zeros(r1,1)]; end
if ~restore, clf, hold on, xys=0; end, view(c1), F=gca;
xys=0;
lw=.5;
n=r1;

if p1==p2 %no distance between spots => no arrow
    H=[];
    return
end
%===============================================================
% Line
if lw>0
        P=zeros(3*n,3); i=1:n;
        P(3*i-2,:)=p1(i,:); P(3*i-1,:)=p2(i,:); P(3*i,:)=NaN;
        H1=plot3(P(:,1),P(:,2),P(:,3),'-','Color',spCol,'LineWidth',lw);
end

%===============================================================
% Scale
ar=get(F,'DataAspectRatio'); ar=sqrt(3)*ar/norm(ar);
%set(F,'DataAspectRatioMode','manual')
if nargin<4 | isempty(w)              % width
    if xys, w=1;
    else
        xr=get(F,'xlim'); yr=get(F,'ylim'); zr=get(F,'zlim');
        w=norm([xr(2)-xr(1),yr(2)-yr(1),zr(2)-zr(1)])/70;
    end
end
if xys, w=w/70; end
if nargin<5 | isempty(h), h=3*w; end  % height

%===============================================================
% Arrowhead
if xys, p1=q1; p2=q2; end
W=(p1-p2)./repmat(ar,n,1);            % new z direction
W=W./repmat(sqrt(sum(W.*W,2)),1,3);   % unit vector
U=[-W(:,2),W(:,1),zeros(n,1)];        % new x direction
N=sqrt(sum(U.*U,2));                  % norm(U)
i=find(N<eps); j=length(i);
U(i,:)=repmat([1,0,0],j,1); N(i)=ones(j,1);
U=U./repmat(N,1,3);                   % unit vector
V=cross(W,U,2);                       % new y direction

m1=30; m2=10; num=200;                % max, min grid spacing, and threshold limits
if n<num, m=round((m2-m1)*n/num+m1);  % adjust grid spacing automatically
else m=m2; end                        % to speed up when matrix size>num

[x,y,z]=cylinder([0,w/2],m);
G=surf(x,y,h*z);
X=get(G,'XData'); Y=get(G,'YData'); Z=get(G,'ZData');
[j,k]=size(X);
c=repmat(spCol,n,1); 
H2=zeros(n,1);
for i=1:n                             % translate, rotate, and scale
    H2(i)=copyobj(G,F);
    newxyz=[X(:),Y(:),Z(:)]*[U(i,:);V(i,:);W(i,:)]*diag(ar)+...
        repmat(p2(i,:),j*k,1);
    newx=reshape(newxyz(:,1),j,k);
    newy=reshape(newxyz(:,2),j,k);
    newz=reshape(newxyz(:,3),j,k);
    if xys
        newx=newx*dx+xr(1); newy=newy*dy+yr(1);
        if xs, newx=10.^newx; end
        if ys, newy=10.^newy; end
    end
    set(H2(i),'XData',newx,'YData',newy,'ZData',newz,...
        'FaceColor',c(i,:),'EdgeColor','None')
end
delete(G)

%===============================================================
% Finish
if ~restore, hold off, end
%if c1==3                              % good for 3D view
 %   set(gcf,'Renderer','OpenGL')
 %end

if nargout, H=[H1(:);H2(:)]; end




