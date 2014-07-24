function [fx fy x_out y_out M pos_u u] = BEM_force_reconstruction3D(planes,forceMesh,E,v,L,x_out,y_out,method)
%this should be generalized to:
%function [fx fy] = BEM_force_reconstruction(x_in,y_in,ux_in,uy_in,forceMesh,E,L,x_out,y_out)


% This has to be done for every z-layer. First check how many there are:
M=[];
u=[];
pos_u=[];
numLayers=length(planes);
for layer=1:numLayers
    x=planes(layer).pos(:,1);
    y=planes(layer).pos(:,2);
    z=planes(layer).z;
    ux=planes(layer).vec(:,1);
    uy=planes(layer).vec(:,2);

    [~, cols]=size(x);

    if cols>1
        x_vec=reshape(x,[],1);
        y_vec=reshape(y,[],1);
        ux_vec=reshape(ux,[],1);
        uy_vec=reshape(uy,[],1);
        u=vertcat(u,ux_vec,uy_vec);
    else
        x_vec=x;
        y_vec=y;
        u=vertcat(u,ux,uy);
    end
    pos_u=vertcat(pos_u,horzcat(x_vec,y_vec));

    %construction of forward map, this takes a long time!

    M_out=calcFwdMapFastBEM3D(x_vec, y_vec, z, forceMesh, E, v);
    M=vertcat(M,M_out);
end


% normalization of basis function will be important when taking the norm!!!
% This has not been considered yet! 
% X = A\B is the solution to the equation AX = B
if nargin == 8 && strcmp(method,'fast')
    sol_coef=(L*eye(2*forceMesh.numBasis)+M'*M)\(M'*u);
else
    sol_coef=(L*eye(2*forceMesh.numNodes)+M'*M)\(M'*u);
end

%if no points of interests are specified, i.e x_out, y_out, then forces are 
%calculated on the nodes of the force mesh:                                                      
if isempty(x_out) || isempty(y_out)
    x_out=forceMesh.p(:,1);
    y_out=forceMesh.p(:,2);
end


%Evaluation of the solution:
if nargin == 8 && strcmp(method,'fast')
    [fx fy x_out y_out]=calcForcesFromCoef(forceMesh,sol_coef,x_out,y_out,'new');
else
    for j=1:2*forceMesh.numNodes
        fx = fx+sol_coef(j)*forceMesh.base(j).f_intp_x(x_out,y_out);
        fy = fy+sol_coef(j)*forceMesh.base(j).f_intp_y(x_out,y_out);
    end
end

