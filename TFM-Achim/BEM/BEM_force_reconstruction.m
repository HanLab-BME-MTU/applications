function [fx, fy, x_out, y_out, M, pos_u, u, sol_coef, sol_mats] = BEM_force_reconstruction(x,y,ux,uy,forceMesh,E,L,x_out,y_out,method,varargin)
% Synopsis  [fx fy x_out y_out M pos_u u sol_coef sol_mats] = BEM_force_reconstruction(x,y,ux,uy,forceMesh,E,L,x_out,y_out,method,meshPtsFwdSol,solMethodBEM)
%
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

% Input check
ip =inputParser;
ip.addRequired('x',@isnumeric);
ip.addRequired('y',@isnumeric);
ip.addRequired('ux',@isnumeric);
ip.addRequired('uy',@isnumeric);
ip.addRequired('forceMesh',@isstruct);
ip.addRequired('E',@isscalar);
ip.addRequired('L',@isscalar);
ip.addRequired('x_out',@isnumeric);%@(x)isscalar(x)||isempty(x));
ip.addRequired('y_out',@isnumeric);%@(x)isscalar(x)||isempty(x));
ip.addRequired('method',@(x)ischar(x)||isempty(x)); % updated in case for BEM
ip.addOptional('meshPtsFwdSol',[],@(x)isscalar(x) ||isempty(x));
ip.addOptional('solMethodBEM','QR',@ischar);
ip.addParamValue('basisClassTblPath','',@ischar);
ip.addParamValue('wtBar',-1,@isscalar);
ip.addParamValue('imgRows',[],@isscalar);
ip.addParamValue('imgCols',[],@isscalar);
ip.parse(x,y,ux,uy,forceMesh,E,L,x_out,y_out,method,varargin{:});
meshPtsFwdSol=ip.Results.meshPtsFwdSol;
solMethodBEM=ip.Results.solMethodBEM;
basisClassTblPath=ip.Results.basisClassTblPath;
wtBar=ip.Results.wtBar;
imgRows = ip.Results.imgRows;
imgCols = ip.Results.imgCols;

if nargin < 12 || isempty(solMethodBEM)
    solMethodBEM='QR';
end

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
    if nargin >= 15
        M=calcFwdMapFastBEM(x_vec, y_vec, forceMesh, E, meshPtsFwdSol,...
            'basisClassTblPath',basisClassTblPath,'wtBar',wtBar,'imgRows',imgRows,'imgCols',imgCols);    
    else
        M=calcFwdMapFastBEM(x_vec, y_vec, forceMesh, E, meshPtsFwdSol,...
            'basisClassTblPath',basisClassTblPath,'wtBar',wtBar);
    end
        
else
    span = 1:length(forceMesh.bounds);
    M=calcFwdMap(x_vec, y_vec, forceMesh, E, span, meshPtsFwdSol);
end
toc;
display('Done: forward map!');



% x = A\B is the solution to the equation Ax = B. The equation we have to
% solve is:
% (L*eyeWeights+ MpM)*sol_coef=(M'*u);

% For L-curve
% [~,n] = size(M);
% 
% B0 = M'*M;
% b  = M'*u;
% 
% bfSigmaRange2 = L*[1e-5 1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4 1e5];
% bfSigmaRange1 = bfSigmaRange2/5;
% bfSigmaRange3 = bfSigmaRange2*5;
% bfSigmaCandidate = [bfSigmaRange1 bfSigmaRange2 bfSigmaRange3];
% bfSigmaCandidate = sort(bfSigmaCandidate);
% coefNorm    = zeros(1,length(bfSigmaCandidate));
% residueNorm = zeros(1,length(bfSigmaCandidate));
% for kk = 1:length(bfSigmaCandidate)
%  B_i = B0 + bfSigmaCandidate(kk)*eye(n);
%  coef_i = B_i\b;
%  coefNorm(kk)    = sqrt(sum(coef_i.^2));
%  residueNorm(kk) = sqrt(sum((M*coef_i-u).^2));
% end
% 
% %Plot the L-curve.
% figure; hold on;
% plot(coefNorm,residueNorm,'.');
% for kk = 1:length(bfSigmaCandidate)
%  text(coefNorm(kk),residueNorm(kk),num2str(bfSigmaCandidate(kk)));
% end
      
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
    
    [normWeights,listNormWeights]=getNormWeights(forceMesh);
    eyeWeights =diag(normWeights);
    
    if length(listNormWeights)==1
        needGSVD =0;
    else
        needGSVD =1;
    end
    % Checked (g)SVD against Matlab inversion! Found no advantage of
    % (g)SVD over Matlab inversion but (g)SVD is much slower. Differences
    % seem to arise only for irregular meshes or at mesh boundaries.
    % For size(M'*M)=7688*7688 I find:
    % M'*M\ =  19sec
    % csvd  = 340sec
    % cgsvd = 561sec
    % Therefore we force the Matlab back slash operator:
%     forceQR=1;
%     forceBackSlash=0;
    
    if ~needGSVD && strcmpi(solMethodBEM,'svd')
        [U,s,V] = csvd(M);
        [sol_coef,~,~] = tikhonov(U,s,V,u,sqrt(L));
        % store these matrices for next frames:
        sol_mats.U=U;
        sol_mats.s=s;
        sol_mats.V=V;
        sol_mats.tool='svd';
    elseif strcmpi(solMethodBEM,'svd') || strcmpi(solMethodBEM,'gsvd')
        % gSVD takes about twice as long as SVD
        [U,sm,X,~] = cgsvd(M,eyeWeights);
        [sol_coef,~,~] = tikhonov(U,sm,X,u,sqrt(L));
        % store these matrices for next frames:
        sol_mats.U =U;
        sol_mats.sm=sm;
        sol_mats.X =X;
        sol_mats.tool='gsvd';
    elseif strcmpi(solMethodBEM,'QR')
        % for a force field with 2*6400 basis function, the residual
        % between the QR-sol and the sol obtained from the backslash
        % operator was: 2.0057e-06 for a mean force magnitude of
        % 85.7. Thus they seem to be numerical identical!
        % accounting for badly scaled linear system - Sangyoon 02/20/13
        % sol_coef = (M'*M+L*D^2)\(M'*u_sol); where D = scaling matrix
        % reference: p11 in Neumaier, Solving ill-conditioned and singular
        % linear systems: a tutorial on regularization
        MpM=M'*M;
        D = diag(sqrt(diag(MpM))); % scaling diagonal matrix
        [Q,R] = qr((MpM+L*D^2));
        u_sol = u;
        sol_coef=R\(Q'*(M'*u_sol));
        sol_mats.Q=Q;
        sol_mats.R=R;
        sol_mats.L=L;
        sol_mats.nW=normWeights;        
        sol_mats.tool='QR';
    elseif strcmpi(solMethodBEM,'backslash') || strcmpi(solMethodBEM,'\')
        % This matrix multiplication takes most of the time. Therefore we
        % store it for later use:
        MpM=M'*M;
        sol_coef=(L*eyeWeights+ MpM)\(M'*u);
        % sol_coef=(L*eyeWeights+M'*M)\(M'*u);
        % store these matrices for next frames:
        sol_mats.MpM=MpM;
        sol_mats.tool='backslash';
    else
        error(['I don''t understand the input for the solution method: ',solMethodBEM])
    end
    % Here we use the identity matrix (all basis classes have equal weight):
    % sol_coef=(L*eye(2*forceMesh.numBasis)+M'*M)\(M'*u);
else
    % normalization of basis function will be important when taking the norm!!!
    % This has not been considered yet!
    MpM=M'*M;
    sol_coef=(L*eye(2*forceMesh.numNodes)+MpM)\(M'*u);
    % store these matrices for next frames:
    sol_mats.MpM=MpM;
    sol_mats.tool='backslash';
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
    fx = zeros(size(x_out));
    fy = zeros(size(y_out));
    for j=1:2*forceMesh.numNodes
        fx = fx+sol_coef(j)*forceMesh.base(j).f_intp_x(x_out,y_out);
        fy = fy+sol_coef(j)*forceMesh.base(j).f_intp_y(x_out,y_out);
    end
end
toc;
display('Done: solution!')