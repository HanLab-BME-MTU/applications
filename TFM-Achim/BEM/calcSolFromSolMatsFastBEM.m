function [pos_f,force,sol_coef]=calcSolFromSolMatsFastBEM(M,sol_mats,u,forceMesh,L,x_out,y_out)

% If M has size mxn then u is a column vector of length m. u is the 
% displacement data on the mesh, that has been used to determine the
% forward map M.

%if no points of interests are specified, i.e x_out, y_out, then forces are 
%calculated on the nodes of the force mesh:                                                      
if nargin<9 || isempty(x_out)
    x_out=forceMesh.p(:,1);
    y_out=forceMesh.p(:,2);
end

% See BEM_force_reconstruction for a nice explanation of the next
% steps:
if strcmp(sol_mats.tool,'svd')
    U=sol_mats.U;
    s=sol_mats.s;
    V=sol_mats.V;
    [sol_coef,~,~] = tikhonov(U,s,V,u,sqrt(L));
elseif strcmp(sol_mats.tool,'gsvd')
    % gSVD takes about twice as long as SVD
    U =sol_mats.U;
    sm=sol_mats.sm;
    X =sol_mats.X;
    [sol_coef,~,~] = tikhonov(U,sm,X,u,sqrt(L));
elseif strcmp(sol_mats.tool,'backslash')
    % This matrix multiplication takes most of the time. Therefore we
    % store it for later use:
    [normWeights]=getNormWeights(forceMesh);
    eyeWeights =diag(normWeights);
    MpM=sol_mats.MpM;
    sol_coef=(L*eyeWeights+ MpM)\(M'*u);
else
    error('No solution matrices have been found')
end

%Evaluation of the solution:
[fx fy x_out y_out]=calcForcesFromCoef(forceMesh,sol_coef,x_out,y_out,'new');
pos_f=horzcat(x_out,y_out);
force=horzcat(   fx,   fy);