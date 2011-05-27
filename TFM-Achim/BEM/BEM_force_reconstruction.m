function [fx fy x_out y_out M pos_u u sol_coef] = BEM_force_reconstruction(x,y,ux,uy,forceMesh,E,L,x_out,y_out,method,meshPtsFwdSol)
% Input :  x,y,ux,uy have to be in the same units, namely
%         pixels. 
% Output: The output fx,fy is actually a surface stress with the same units
%         as the input E! In particular, the unit of the output force is
%         independent of the units of the input x,y,ux,uy.
%         The reason for this is essentially that the elastic stress is
%         only dependent on the non-dimensional strain which is given by
%         spatial derivatives of the displacements, that is du/dx. If u and
%         dx (essentially cluster_size) are in the same units, then the
%         resulting force has the same dimension as the input E.

%         u: is the measured displacement! (not the model u!)     

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

%construction of forward map, this takes a long time!

display('2.) Building up forward map:...');
tic;
if nargin >= 10 && strcmp(method,'fast')
    M=calcFwdMapFastBEM(x_vec, y_vec, forceMesh, E, meshPtsFwdSol);    
else
    M=calcFwdMap(x_vec, y_vec, forceMesh, E, meshPtsFwdSol);
end
toc;
display('Done: forward map!');



% X = A\B is the solution to the equation AX = B
tic;
display('3.) Solve for coefficients, this is memory intensive [~5min]:... ')
if nargin >= 10 && strcmp(method,'fast')
    % If there are more than one basis function class, then correct
    % weighting of basis function according to their volume is of course
    % important when taking the norm!!! If there is only one basis function
    % class, then the weights are all one (e.g. square lattice where
    % boundary nodes are skipped) 
    % See refine_BEM_force_reconstruction for a nice explanation of the next
    % steps:
    weights    =vertcat(forceMesh.basis(:).unitVolume); % volume of the basis function
    repWeights =repmat(weights(:),2,1); % the basis function for x/y comp. have the same weight
    normWeights=repWeights/max(repWeights); % normalize it with max value.
    eyeWeights =diag(normWeights);    
    weightList=unique(normWeights);
    
    
    if length(weightList)==1
        doSVD =1;
        doGSVD=0;
    else
        doSVD =0;
        doGSVD=1;
    end
    % Checked (g)SVD against Matlab inversion! Found no advantage of
    % (g)SVD over Matlab inversion but (g)SVD is much slower. Differences
    % seem to arise only for irregular meshes or at mesh boundaries.
    % For size(M'*M)=7688*7688 I find:
    % M'*M\ =  19sec
    % csvd  = 340sec
    % cgsvd = 561sec
    % Therefore we force the Matlab back slash operator:
    forceBackSlash=1;
    
    if doSVD && ~forceBackSlash        
        tic;
        [U,s,V] = csvd(M);
        [sol_coef,rho,eta] = tikhonov(U,s,V,u,sqrt(L));
        toc;
    elseif doGSVD && ~forceBackSlash
        % gSVD takes about twice as long as SVD
        tic;
        [U,sm,X,V] = cgsvd(M,eyeWeights);
        [sol_coef,rho,eta] = tikhonov(U,sm,X,u,sqrt(L));
        toc;
    else
        sol_coef=(L*eyeWeights+M'*M)\(M'*u);
    end
    % Here we use the identity matrix (all basis classes have equal weight):
    % sol_coef=(L*eye(2*forceMesh.numBasis)+M'*M)\(M'*u);
else
    % normalization of basis function will be important when taking the norm!!!
    % This has not been considered yet! 
    sol_coef=(L*eye(2*forceMesh.numNodes)+M'*M)\(M'*u);
end
toc;
display('Done: coefficients!');


%if no points of interests are specified, i.e x_out, y_out, then forces are 
%calculated on the nodes of the force mesh:                                                      
if nargin<9 || isempty(x_out)
    x_out=forceMesh.p(:,1);
    y_out=forceMesh.p(:,2);
end

%Evaluation of the solution:
display('4.) Evaluate solution:... ')
tic;
if nargin >= 10 && strcmp(method,'fast')
    [fx fy x_out y_out]=calcForcesFromCoef(forceMesh,sol_coef,x_out,y_out,'new');
else
    for j=1:2*forceMesh.numNodes
        fx = fx+sol_coef(j)*forceMesh.base(j).f_intp_x(x_out,y_out);
        fy = fy+sol_coef(j)*forceMesh.base(j).f_intp_y(x_out,y_out);
    end
end
toc;
display('Done: solution!')