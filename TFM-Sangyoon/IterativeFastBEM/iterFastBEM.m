function [fx, fy, x_out, y_out, M, pos_u, u, sol_coef, sol_mats] = iterFastBEM(x,y,ux,uy,forceMesh,E,L,x_out,y_out,method,varargin)
% Synopsis  [fx fy x_out y_out M pos_u u sol_coef sol_mats] = iterFastBEM(x,y,ux,uy,forceMesh,E,L,x_out,y_out,method,meshPtsFwdSol,solMethodBEM)
% This function uses sparse matrix for forward map to be able to take care
% of dense grid.
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
% Sangyoon Han (March 2013)

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
ip.addOptional('fwdMap',[],@isnumeric);
ip.parse(x,y,ux,uy,forceMesh,E,L,x_out,y_out,method,varargin{:});
meshPtsFwdSol=ip.Results.meshPtsFwdSol;
solMethodBEM=ip.Results.solMethodBEM;
basisClassTblPath=ip.Results.basisClassTblPath;
wtBar=ip.Results.wtBar;
imgRows = ip.Results.imgRows;
imgCols = ip.Results.imgCols;
M = ip.Results.fwdMap;

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

display('2.) Building up sparse forward map:...');
tic;
if nargin >= 10 && strcmp(method,'fast') && isempty(M)
    if nargin >= 15
        M=calcFwdMapFastBEMSparse(x_vec, y_vec, forceMesh, E, meshPtsFwdSol,...
            'basisClassTblPath',basisClassTblPath,'wtBar',wtBar,'imgRows',imgRows,'imgCols',imgCols);    
    else
        M=calcFwdMapFastBEMSparse(x_vec, y_vec, forceMesh, E, meshPtsFwdSol,...
            'basisClassTblPath',basisClassTblPath,'wtBar',wtBar);
    end
elseif isempty(M)
    span = 1:length(forceMesh.bounds);
    M=calcFwdMap(x_vec, y_vec, forceMesh, E, span, meshPtsFwdSol);
else
    display('Using input Forward Map');
end
toc;
display('Done: forward map!');



% x = A\B is the solution to the equation Ax = B. The equation we have to
% solve is:
% (L*eyeWeights+ MpM)*sol_coef=(M'*u);

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
%         D = diag(sqrt(diag(MpM))); % scaling diagonal matrix
%         [Q,R] = qr((MpM+L*D^2));
%         L = chooseRegParamFromLCurve(M,MpM,u,L,eyeWeights);
        [Q,R] = qr((MpM+L*eyeWeights));
        sol_coef=R\(Q'*(M'*u));
        sol_mats.Q=Q;
        sol_mats.R=R;
        sol_mats.L=L;
        sol_mats.nW=normWeights;        
        sol_mats.tool='QR';
    elseif strcmpi(solMethodBEM,'QRscaled')
        % accounting for badly scaled linear system - Sangyoon 02/20/13
        % sol_coef = (M'*M+L*D^2)\(M'*u_sol); where D = scaling matrix
        % reference: p11 in Neumaier, Solving ill-conditioned and singular
        % linear systems: a tutorial on regularization
        MpM=M'*M;
        D = diag(sqrt(diag(MpM))); % scaling diagonal matrix
%         D = D./normest(D);
        % For L-curve
%         L = chooseRegParamFromLCurve(M,MpM,u,L,D);

        [Q,R] = qr((MpM+L*D^2));
        sol_coef=R\(Q'*(M'*u));
        sol_mats.Q=Q;
        sol_mats.R=R;
        sol_mats.L=L;
        sol_mats.nW=normWeights;        
        sol_mats.tool='QRscaled';
    elseif strcmpi(solMethodBEM,'LaplacianRegularization')
        % second order tikhonov regularization (including diagonal)
        % make Lap matrix
%         nBeads = round(size(M,1)/2);
        nBasisx = size(unique(forceMesh.p(:,1)),1);
        nBasisy = size(unique(forceMesh.p(:,2)),1);
        nBasis = nBasisx*nBasisy;
        k=1;
        % this is for assuring the boundary to be considered for being
        % penalized for regularization
        display('Building Laplacian Map...')
        tic;
        Lap = sparse(2*nBasis,2*nBasis);
        for ii=1:nBasisx
            tempLapx = zeros(nBasisy,nBasis); %for parfor constraints
            tempLapy = zeros(nBasisy,nBasis); %for parfor constraints
            parfor jj=1:nBasisy
                Lap2D = zeros(nBasisy,nBasisx);
                if ii==1 || jj==1 || ii==nBasisx || jj==nBasisy
%                     Lap2D(ii,jj)=1; % 0th order regularization. Potentially 
%                     % this can be improved with basis function convolution
                    if ii==1 
                        if jj==1
                            Lap2D(ii,jj)=-2;
                            Lap2D(ii+1,jj)=1;
                            Lap2D(ii,jj+1)=1;
                        elseif jj==nBasisy
                            Lap2D(ii,jj)=-2;
                            Lap2D(ii+1,jj)=1;
                            Lap2D(ii,jj-1)=1;
                        else
                            Lap2D(ii,jj)=-2;
                            Lap2D(ii+1,jj)=1;
                            Lap2D(ii,jj+1)=1;
                        end
                    elseif ii==nBasisx
                        if jj==1
                            Lap2D(ii,jj)=-2;
                            Lap2D(ii-1,jj)=1;
                            Lap2D(ii,jj+1)=1;
                        elseif jj==nBasisy
                            Lap2D(ii,jj)=-2;
                            Lap2D(ii-1,jj)=1;
                            Lap2D(ii,jj-1)=1;
                        else
                            Lap2D(ii,jj)=-2;
                            Lap2D(ii-1,jj)=1;
                            Lap2D(ii,jj-1)=1;
                        end
                    elseif jj==1
                        if ii~=nBasisx && ii~=1
                            Lap2D(ii,jj)=-2;
                            Lap2D(ii+1,jj)=1;
                            Lap2D(ii,jj+1)=1;
                        end
                    elseif jj==nBasisy                        
                        if ii~=nBasisx && ii~=1
                            Lap2D(ii,jj)=-2;
                            Lap2D(ii+1,jj)=1;
                            Lap2D(ii,jj-1)=1;
                        end
                    end                        
                else
                    % diagonal laplacian
                    Lap2D(ii,jj) = -6;
                    Lap2D(ii,jj+1) = 1;
                    Lap2D(ii,jj-1) = 1;
                    Lap2D(ii+1,jj) = 1;
                    Lap2D(ii-1,jj) = 1;
                    Lap2D(ii+1,jj+1) = 0.5;
                    Lap2D(ii-1,jj+1) = 0.5;
                    Lap2D(ii+1,jj-1) = 0.5;
                    Lap2D(ii-1,jj-1) = 0.5;
    %                 % orthogonal laplacian
    %                 Lap2D(ii,jj) = -4;
    %                 Lap2D(ii,jj+1) = 1;
    %                 Lap2D(ii,jj-1) = 1;
    %                 Lap2D(ii+1,jj) = 1;
    %                 Lap2D(ii-1,jj) = 1;
                end
                tempLapx(jj,:) = reshape(Lap2D,nBasis,1)';
%                 Lap(k,nBasis+1:2*nBasis) = reshape(Lap2D,nBasis,1)';
%                 Lap(k+(nBasisx-2)*(nBasisy-2),1:nBasis) = reshape(Lap2D,nBasis,1)';
                tempLapy(jj,:) = reshape(Lap2D,nBasis,1)';
            end
            Lap((ii-1)*nBasisy+1:(ii-1)*nBasisy+nBasisx,1:nBasis) = tempLapx;
            Lap((ii-1)*nBasisy+nBasis+1:(ii-1)*nBasisy+nBasis+nBasisx,nBasis+1:2*nBasis) = ...
                tempLapy;
        end
        toc
%         Lap = zeros((nBasisx-2)*(nBasisy-2)*2,2*nBeads);
%         for ii=2:nBasisx-1
%             for jj=2:nBasisy-1
%                 Lap2D = zeros(nBasisy,nBasisx);
%                 % diagonal laplacian
%                 Lap2D(ii,jj) = -6;
%                 Lap2D(ii,jj+1) = 1;
%                 Lap2D(ii,jj-1) = 1;
%                 Lap2D(ii+1,jj) = 1;
%                 Lap2D(ii-1,jj) = 1;
%                 Lap2D(ii+1,jj+1) = 0.5;
%                 Lap2D(ii-1,jj+1) = 0.5;
%                 Lap2D(ii+1,jj-1) = 0.5;
%                 Lap2D(ii-1,jj-1) = 0.5;
% %                 % orthogonal laplacian
% %                 Lap2D(ii,jj) = -4;
% %                 Lap2D(ii,jj+1) = 1;
% %                 Lap2D(ii,jj-1) = 1;
% %                 Lap2D(ii+1,jj) = 1;
% %                 Lap2D(ii-1,jj) = 1;
%                 Lap(k,1:nBasis) = reshape(Lap2D,nBasis,1)';
% %                 Lap(k,nBasis+1:2*nBasis) = reshape(Lap2D,nBasis,1)';
% %                 Lap(k+(nBasisx-2)*(nBasisy-2),1:nBasis) = reshape(Lap2D,nBasis,1)';
%                 Lap(k+(nBasisx-2)*(nBasisy-2),nBasis+1:2*nBasis) = reshape(Lap2D,nBasis,1)';
%                 k=k+1;
%             end
%         end
        MpM=M'*M;

        % For L-curve
%         L = chooseRegParamFromLCurve(M,MpM,u,L,Lap);
        
        [Q,R] = qr((MpM+L*(Lap'*Lap)));
        sol_coef=R\(Q'*(M'*u));
        sol_mats.Q=Q;
        sol_mats.R=R;
        sol_mats.L=L;
        sol_mats.nW=normWeights;
        sol_mats.tool='LaplacianRegularization';
    elseif strcmpi(solMethodBEM,'1NormReg')
        % Now, perform the sparse deconvolution.
        disp('Performing sparse deconvolution; adoped from Aster et. al.')

        % plot the solution for the corner
        sol_coef = iterativeL1Regularization(M,u,eyeWeights,L,400); %400=maximum iteration number
        sol_mats.nW=normWeights;
        sol_mats.tool='1Norm-0thReg';
    elseif strcmpi(solMethodBEM,'CG')
        % Conjugate Gradient
        ith = 30;
        [X,rho,eta]=cgls(M,u,100); 
        sol_coef=X(:,ith); % 
        display(['residual norm is ' num2str(rho(ith)) ', and coeff norm is ' num2str(eta(ith)) ' for ' num2str(ith) 'th iterations.'])
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
%     [normWeights,~]=getNormWeights(forceMesh);
%     eyeWeights =diag(normWeights);

    MpM=M'*M;
    sol_coef=(L*eye(2*forceMesh.numNodes)+MpM)\(M'*u);
%     sol_coef=(L*eyeWeights+MpM)\(M'*u);
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
    [fx,fy,x_out,y_out]=calcForcesFromCoef(forceMesh,sol_coef,x_out,y_out,'new');
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

function Lout = chooseRegParamFromLCurve(M,MpM,u,L,Lap)
display('Calculating L-curve ...')
b  = M'*u;
bfSigmaRange2 = L*[1e-5 1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4 1e5];
bfSigmaRange1 = bfSigmaRange2/5;
bfSigmaRange3 = bfSigmaRange2*5;
bfSigmaCandidate = [bfSigmaRange1 bfSigmaRange2 bfSigmaRange3];
bfSigmaCandidate = sort(bfSigmaCandidate);
coefNorm    = zeros(1,length(bfSigmaCandidate));
residueNorm = zeros(1,length(bfSigmaCandidate));
Lap2 = Lap'*Lap;
for kk = 1:length(bfSigmaCandidate)
 B_i = MpM + bfSigmaCandidate(kk)*Lap2;
 coef_i = B_i\b;
 coefNorm(kk)    = sqrt(sum((Lap*coef_i).^2)); % Solution norm ||Lm||
 residueNorm(kk) = sqrt(sum((M*coef_i-u).^2)); % Residual norm ||Gm-d||
end

%Plot the L-curve.
figure; hold on;
plot(residueNorm,coefNorm,'.'); xlabel('Residual norm ||M*coef-u||'); 
ylabel('Solution seminorm ||L*coef||');
for kk = 1:length(bfSigmaCandidate)
 text(residueNorm(kk),coefNorm(kk),num2str(bfSigmaCandidate(kk)));
end

options.WindowStyle='normal';
answer = inputdlg('Please identify the corner:','Input for corner',1,{num2str(L)},options);
Lout = str2double(answer{1});

function mreg=iterativeL1Regularization(G,d,L,alpha,maxiter,tolx,tolr)
% adopted from Aster et al.
% Default for tolr=1.0e-6
if (nargin < 7)
  tolr=1.0e-6;
end

% Default for tolx=1.0e-4;
if (nargin < 6)
  tolx=1.0e-4;
end

% Default for maxiter=100
if (nargin < 5)
  maxiter=100;
end

% unchanging constants in the system that is solved repeatedly
GTG=G'*G;
GTd=G'*d;
% Start with an initial unweighted solution.
iter=1;
display(['L1 Norm regularization (iteration:' num2str(iter) ')'])
tic
m=(2*GTG+alpha*(L'*L))\(2*G'*d);
toc
% iterate until maxiter or we converge and return
while (iter < maxiter)
  iter=iter+1;

  % get get the magnitude of Lm, but don't let any element be less than tolr
  absLm=abs(L*m);
  absLm(absLm<tolr)=tolr;

  % build the diagonal weighting matrix for this iteration
  R=diag(1./absLm);

  mold=m;
  
%   % use LSQR to get m that minimizes argmin||[            M              ]f - [d]||2
%   %                                           sqrt(alpha/2)*sqrt(absLm)*L      0
%   display('LSQR...')
%   tic
%   A_bottom = sqrt(alpha/2)*sqrt(R)*L;
%   A = [G;A_bottom];
%   b = [d;zeros(size(A_bottom,1),1)];
%   m = lsqr(A,b,tolx,maxiter);
%   toc
%   % get the new iterate and check for convergance
%   display('QR...')
%   tic
%   [Q,Rr] = qr(2*GTG+alpha*L'*R*L);
%   m=Rr\(Q'*(2*GTd));
%   toc
%   display(norm(m))
  display(['backslash... (iteration:' num2str(iter) ')'])
  tic
  m=(2*GTG+alpha*L'*R*L)\(2*GTd);
  toc
%   display(norm(m))

  if (norm(m-mold)/(1+norm(mold)) < tolx)
    mreg=m;
    return
  end
end

% Give a warning, if desired, but return best solution.
%warning('irlslreg1 maximum iterations exceeded.');
mreg=m;

function [X,rho,eta]=cgls(G,d,niter)
% Parameter Estimation and Inverse Problems, 2nd edition, 2011 
% by R. Aster, B. Borchers, C. Thurber
%
% [X,rho,eta]=cgls(G,d,niter)
%
% Performs niter iterations of the CGLS algorithm on the least
% squares problem
% 
%   min norm(G*m-d)
%
% The iterates 1, 2, ..., niter are returned in the columns of the
% matrix X.  For each iterate we also compute rho(i)=norm(G*m-d)
% and eta(i)=norm(m).

% Figure out problem size.
[nrows,ncols]=size(G);
if (length(d) ~= nrows)
  error('G and d do not match in size.');
end

% Setup space for the results.
X=zeros(ncols,niter);
rho=zeros(niter,1);
eta=zeros(niter,1);

% Setup for the first iteration.
m=zeros(ncols,1);
p=zeros(ncols,1);
beta=0;
s=d;
r=G'*s;

% Main loop- perform CGLS iterations.
for k=0:niter-1
  % We'll precompute r'*r since it's used in several places.
  rtr=r'*r;

  %  Update beta.
  if (k>0)
    beta=rtr/(prevr'*prevr);
  end

  %  Update p
  p=r+beta*p;

  % Compute the new alpha.  To avoid doing the matrix vector
  % multiplication repeatedly, we store G*p in Gp
  Gp=G*p;
  alpha=rtr/(Gp'*Gp);

  % Update m.
  m=m+alpha*p;

  % Update s.
  s=s-alpha*Gp;

  % Save r for the next iteration, and then update it.
  prevr=r;
  r=G'*s;

  % Store the new iterate.
  X(:,k+1)=m;
  rho(k+1)=norm(s);
  eta(k+1)=norm(m);
end