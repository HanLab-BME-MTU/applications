function [u1,u2,exx,eyy,exy,sxx,syy,sxy,u1x,u2x,u1y,u2y]=postProcSol(u,p,t)
% To calculate the stress we will neeed the Poisson ratio and the LOCAL
% Young's modulus at the node points p:
E=gfun_youngs_mod(p(1,:)',p(2,:)');
v=gfun_poisson_ratio;


np=length(p);
nt=length(t);

% determine the x- and y-component of the solution:
u1=u(1:np);
u2=u(np+1:end);

% calculate the gradient. The gardient is returnd at the center of each
% triangle:
[u1x,u1y]=pdegrad(p,t,u1);
[u2x,u2y]=pdegrad(p,t,u2);

% first interpolate to nodes:
u1x=pdeprtni(p,t,u1x);
u1y=pdeprtni(p,t,u1y);
u2x=pdeprtni(p,t,u2x);
u2y=pdeprtni(p,t,u2y);

% calculate the strain:
exx=u1x;
eyy=u2y;
exy=1/2*(u1y+u2x);

% calculate the stress:
sxx=E./(1-v.^2).*(  exx+v.*eyy          );
syy=E./(1-v.^2).*(v*exx+   eyy          );
sxy=E./(1-v.^2).*(            +(1-v)*exy);


% Alternatively we could interpolate u to the triangle midpoints and return
% all functions evaluated at triangular midpoints. The midpoints could be
% obtained by:
% xpos_n=p(1,:);
% ypos_n=p(2,:);
% xpos_t=pdeintrp(p,t,xpos_n');
% ypos_t=pdeintrp(p,t,ypos_n');
% plot(xpos_t,ypos_t,'.r')
% u1_t=pdeintrp(p,t,u1);
% u2_t=pdeintrp(p,t,u2);


% interpolate between functions defined at triangle nodes and functions
% defined at triangle midpoints: 
%pdeintrp 
%pdeprtni  

% interpolates a functions from a triangular mesh to a rectangular grid:
%tri2grid 

% compute gradients of the solution:
%pdegrad 
%pdecgrad 

% large number of options for plotting the solution:
%pdeplot 

% convenient shorthands for pdeplot:
%pdecont 
%pdesurf 