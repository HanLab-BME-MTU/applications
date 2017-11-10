function [fx fy x_out y_out sol_coef]=calcSolFromFwdMap(M,u,forceMesh,L,x_out,y_out)
% Normalization of basis function will be important when taking the norm!!!
% This has not been considered yet! 
% If M has size mxn then u is a column vector of length m. u is the 
% displacement data on the mesh, that has been used to determine the
% forward map M.

% X = A\B is the solution to the equation AX = B
sol_coef=(L*eye(2*forceMesh.numNodes)+M'*M)\(M'*u);

%if no points of interests are specified, i.e x_out, y_out, then forces are 
%calculated on the nodes of the force mesh:                                                      
if nargin<6
    x_out=forceMesh.p(:,1);
    y_out=forceMesh.p(:,2);
end

fx=0*x_out;
fy=0*y_out;

%reconstructed forces are obtained by multiplying the ceof with the
%resepctive basis.
%This can be done much quicker, since most of the functions evaluate to zero 
%e.g. fx(x0,y0) is non-zero only for 1:forceMesh.numNodes and base functions
%base(j).f_intp_x which have support that comprises (x0,y0)
for j=1:2*forceMesh.numNodes
        fx = fx+sol_coef(j)*forceMesh.base(j).f_intp_x(x_out,y_out);
        fy = fy+sol_coef(j)*forceMesh.base(j).f_intp_y(x_out,y_out);
end