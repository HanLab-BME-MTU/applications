function [fx fy x_out y_out M pos_u u sol_coef] = refine_BEM_force_reconstruction(x,y,ux,uy,M_old,refinedForceMesh,E,L,meshPtsFwdSol,x_out,y_out)

[~, cols]=size(x);

if cols>1
    x_vec=reshape(x,[],1);
    y_vec=reshape(y,[],1);
    ux_vec=reshape(ux,[],1);
    uy_vec=reshape(uy,[],1);
    u=vertcat(ux_vec,uy_vec);
else
    x_vec=x;
    y_vec=y;
    u=vertcat(ux,uy);
end
pos_u=horzcat(x_vec,y_vec);

% refine the forward map, this might take a long time!
% here, I only support the fast BEM:
M=refineFwdMapFastBEM(x_vec, y_vec, M_old, refinedForceMesh, E, meshPtsFwdSol);

% In case of L>0 and if the basis function have differently sized support
% then, the basis functions have the be appropriately weighted in the
% regularization term according to their unit volume. E.g. denser sampling
% in region of constant stress must not change the result. Thus L*|F|^2
% should be the same!
% Get the unit volume for each basisfunction:
weights  =vertcat(refinedForceMesh.basis(:).unitVolume);
% We have the same weights for the x and y components of the basis functions:
repWeights=repmat(weights(:),2,1);
% Normalize the weights with the largest base-function:
normWeights=repWeights/max(repWeights);
% Create the diagonal matrix with the weights, if the support is the same
% for all basis functions, then this will end up to be the identity matrix
% as has been used before [eye(2*refinedForceMesh.numBasis)]:
eyeWeights=diag(normWeights);

% X = A\B is the solution to the equation AX = B
tic
display('Solve for the coefficients, this requires a lot of memory... ')
sol_coef=(L*eyeWeights+M'*M)\(M'*u);
toc

% Evaluation of the solution:
tic
display('Evaluation of the solution... ')
[fx fy x_out y_out]=calcForcesFromCoef(refinedForceMesh,sol_coef,x_out,y_out,'new');
toc