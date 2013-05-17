function [pos_f,force,sol_coef]=calcSolFromSolMatsFastBEM(M,sol_mats,u,forceMesh,L,x_out,y_out)

% If M has size mxn then u is a column vector of length m. u is the 
% displacement data on the mesh, that has been used to determine the
% forward map M.

%if no points of interests are specified, i.e x_out, y_out, then forces are 
%calculated on the nodes of the force mesh:                                                      
if nargin<6 || isempty(x_out)
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
elseif strcmp(sol_mats.tool,'QR')
%     [normWeights]=getNormWeights(forceMesh);
%     eyeWeights = sol_mats.eyeWeights;
    sol_nW=sol_mats.nW;
    sol_L =sol_mats.L;
    % check that regularization parameter and weights have not changed
    % (since Q,R have been calculated for a certain set of reg. par. and
    % weights!). But this should always be the case:
    if sol_L==L %&& sum(sol_nW~=normWeights)==0
        Q=sol_mats.Q;
        R=sol_mats.R;
        sol_coef=R\(Q'*(M'*u));
    else
        error('Weights or regularization parameter have been changed. QR cannot be reused!')
    end
elseif strcmp(sol_mats.tool,'1NormReg')
    eyeWeights =sol_mats.eyeWeights;
    MpM=sol_mats.MpM;
    M=sol_mats.M;
    L=sol_mats.L;
    maxIter = sol_mats.maxIter;
    tolx = sol_mats.tolx;
    tolr = sol_mats.tolr;
    sol_coef = iterativeL1Regularization(M,MpM,u,eyeWeights,L,maxIter,tolx,tolr); 
elseif strcmpi(sol_mats.tool,'1NormRegLaplacian')
    % Now, perform the sparse deconvolution.
    Lap = sol_mats.Lap;
    % plot the solution for the corner
    MpM=sol_mats.MpM;
    maxIter = sol_mats.maxIter ;
    tolx = sol_mats.tolx;
    tolr = sol_mats.tolr;
    L=sol_mats.L;
%         [sol_coef,L] = calculateLfromLcurveSparse(M,MpM,u,Lap,maxIter,tolx,tolr,solMethodBEM);
    sol_coef = iterativeL1Regularization(M,MpM,u,L,-Lap,maxIter,tolx,tolr); %400=maximum iteration number
elseif strcmpi(sol_mats.tool,'LaplacianReg')
    % second order tikhonov regularization (including diagonal)
    % make Lap matrix
%         nBeads = round(size(M,1)/2);
    Lap = sol_mats.Lap;
    MpM=sol_mats.MpM;
    M=sol_mats.M;
    L=sol_mats.L;
    % For L-curve
%         [sol_coef,L] = calculateLfromLcurve(M,MpM,u,Lap,solMethodBEM);
    sol_coef=(-L*Lap+ MpM)\(M'*u);
elseif strcmp(sol_mats.tool,'backslash')
    % This matrix multiplication takes most of the time. Therefore we
    % store it for later use:
%     [normWeights]=getNormWeights(forceMesh);
%     eyeWeights =diag(normWeights);
    eyeWeights = sol_mats.eyeWeights;
    MpM=sol_mats.MpM;
    sol_coef=(L*eyeWeights+ MpM)\(M'*u);
else
    error('No solution matrices have been found')
end

%Evaluation of the solution:
[fx fy x_out y_out]=calcForcesFromCoef(forceMesh,sol_coef,x_out,y_out,'new');
pos_f=horzcat(x_out,y_out);
force=horzcat(   fx,   fy);