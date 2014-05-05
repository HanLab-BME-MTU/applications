function [fx, fy, x_out, y_out, M, pos_u, u, sol_coef, sol_mats] = ...
    BEM_force_reconstruction(x,y,ux,uy,forceMesh,E,L,x_out,y_out,method,varargin)
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
% Achim Besser 2011
% Sangyoon Han 2013

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
ip.addParamValue('LcurveDataPath','',@ischar);
ip.addParamValue('LcurveFigPath','',@ischar);
ip.addParamValue('LcurveFactor','',@isscalar);
ip.addParamValue('wtBar',-1,@isscalar);
ip.addParamValue('imgRows',[],@isscalar);
ip.addParamValue('imgCols',[],@isscalar);
ip.addOptional('fwdMap',[],@isnumeric);
ip.addParamValue('thickness',472,@isscalar); % default assuming 34 um with 72 nm/pix resolution
ip.addParamValue('useLcurve',false,@islogical); % default assuming 34 um with 72 nm/pix resolution
ip.addParamValue('paxImg',[],@ismatrix);
ip.parse(x,y,ux,uy,forceMesh,E,L,x_out,y_out,method,varargin{:});
meshPtsFwdSol=ip.Results.meshPtsFwdSol;
solMethodBEM=ip.Results.solMethodBEM;
basisClassTblPath=ip.Results.basisClassTblPath;
LcurveDataPath=ip.Results.LcurveDataPath;
LcurveFigPath=ip.Results.LcurveFigPath;
LcurveFactor=ip.Results.LcurveFactor;
wtBar=ip.Results.wtBar;
imgRows = ip.Results.imgRows;
imgCols = ip.Results.imgCols;
M = ip.Results.fwdMap;
thickness = ip.Results.thickness;    
paxImage = ip.Results.paxImg;
useLcurve = ip.Results.useLcurve;    

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
    ux_vec=ux;
    uy_vec=uy;
    u=vertcat(ux,uy);
end
pos_u=horzcat(x_vec,y_vec);

%construction of forward map, this takes a long time!

display('2.) Building up forward map:...');
tic;
if nargin >= 10 && strcmp(method,'fast') && isempty(M)
    if nargin >= 15
        M=calcFwdMapFastBEM(x_vec, y_vec, forceMesh, E, meshPtsFwdSol,...
            'basisClassTblPath',basisClassTblPath,'wtBar',wtBar,'imgRows',imgRows,'imgCols',imgCols,'thickness',thickness);    
    else
        M=calcFwdMapFastBEM(x_vec, y_vec, forceMesh, E, meshPtsFwdSol,...
            'basisClassTblPath',basisClassTblPath,'wtBar',wtBar);
    end
elseif isempty(M)
    span = 1:length(forceMesh.bounds);
    M=calcFwdMap(x_vec, y_vec, forceMesh, E, span, meshPtsFwdSol,'conv_free');
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
    % steps: <-Achim's comment
    
    % It is more correct to consider the weight to be Cholesky factor of /int hi(x) hj(x) dx
    
%     [normWeights,listNormWeights]=getNormWeights(forceMesh);
%     eyeWeights =diag(normWeights);
    
%     if length(listNormWeights)==1
        needGSVD =0;
%     else
%         needGSVD =1;
%     end
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
        [eyeWeights,~] =getGramMatrix(forceMesh);
        % make matrix with paxImage at the basis nodes
        MpM=M'*M;
        if ~isempty(paxImage)
            paxWeights = getPaxWeights(forceMesh,paxImage,x_vec,y_vec,ux_vec,uy_vec);
            [Q,R] = qr((MpM+L*eyeWeights.*paxWeights));
            sol_coef=R\(Q'*(M'*u));
            sol_mats.eyeWeights=eyeWeights;
            sol_mats.L=L;
        else
            if useLcurve
                [~,L] = calculateLfromLcurve(M,MpM,u,eyeWeights,LcurveDataPath,LcurveFigPath);
            end
            [Q,R] = qr((MpM+L*eyeWeights));
            sol_coef=R\(Q'*(M'*u));
            sol_mats.Q=Q;
            sol_mats.R=R;
            sol_mats.L=L;
        end            
        sol_mats.tool='QR';
    elseif strcmpi(solMethodBEM,'LaplacianReg')
        % second order tikhonov regularization (including diagonal)
        % make Lap matrix
        Lap = buildLaplacian(forceMesh);
        MpM=M'*M;
        % For L-curve
%         [sol_coef,L] = calculateLfromLcurve(M,MpM,u,Lap,solMethodBEM);
        sol_coef=(-L*Lap+ MpM)\(M'*u);
        % store these matrices for next frames:
        sol_mats.M=M;
        sol_mats.MpM=MpM;
        sol_mats.Lap=Lap;
        sol_mats.L=L;
        sol_mats.tool='LaplacianReg';
    elseif strcmpi(solMethodBEM,'1NormReg')
        % Now, perform the sparse deconvolution.
        disp('Performing sparse deconvolution; adoped from Aster et. al.')

        [eyeWeights,~] =getGramMatrix(forceMesh);
        % plot the solution for the corner
        tolx =  sqrt(forceMesh.numBasis)*5e-3; % This will make tolx sensitive to overall number of nodes. (rationale: the more nodes are, 
        % the larger tolerance should be, because misfit norm can be larger out of more nodes).
        disp(['tolerance value: ' num2str(tolx)])
        MpM=M'*M;
        maxIter = 10;
        tolr = 1e-7;
        if useLcurve
            disp('L-curve ...')
            [sol_coef,L] = calculateLfromLcurveSparse(L,M,MpM,u,eyeWeights,maxIter,tolx,tolr,solMethodBEM,LcurveDataPath,LcurveFigPath,LcurveFactor);
        else
            sol_coef = iterativeL1Regularization(M,MpM,u,eyeWeights,L,maxIter,tolx,tolr); 
        end
%         sol_mats.nW=normWeights;
        sol_mats.eyeWeights=eyeWeights;
        sol_mats.L=L;
        sol_mats.M = M;
        sol_mats.MpM = MpM;
        sol_mats.maxIter = maxIter;
        sol_mats.tolx = tolx;
        sol_mats.tolr = tolr;
        sol_mats.tool='1NormReg';
    elseif strcmpi(solMethodBEM,'1NormRegLaplacian')
        % Now, perform the sparse deconvolution.
        disp('Performing sparse deconvolution with laplacian')

        disp('Building laplacian operator')
        Lap = buildLaplacian(forceMesh);
        % plot the solution for the corner
        MpM=M'*M;

        maxIter = 10;
        tolx = 2e-2;
        tolr = 1e-7;
%         [sol_coef,L] = calculateLfromLcurveSparse(M,MpM,u,Lap,maxIter,tolx,tolr,solMethodBEM);
        sol_coef = iterativeL1Regularization(M,MpM,u,L,-Lap,maxIter,tolx,tolr); %400=maximum iteration number
        sol_mats.L=L;
        sol_mats.Lap = Lap;
        sol_mats.M = M;
        sol_mats.MpM = MpM;
        sol_mats.maxIter = maxIter;
        sol_mats.tolx = tolx;
        sol_mats.tolr = tolr;
        sol_mats.tool='1NormRegLaplacian';
    elseif strcmpi(solMethodBEM,'backslash') || strcmpi(solMethodBEM,'\')
        % This matrix multiplication takes most of the time. Therefore we
        % store it for later use:
        [eyeWeights,~] =getGramMatrix(forceMesh);
        MpM=M'*M;
        [sol_coef,L] = calculateLfromLcurve(M,MpM,u,eyeWeights,LcurveDataPath,LcurveFigPath);
%         sol_coef=(L*eyeWeights+ MpM)\(M'*u);
        % store these matrices for next frames:
        sol_mats.eyeWeights=eyeWeights;
        sol_mats.MpM=MpM;
        sol_mats.L=L;
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

function Lap = buildLaplacian(forceMesh)
nBasis = forceMesh.numBasis;
basisx=zeros(nBasis,1);
basisy=zeros(nBasis,1);
for k=1:nBasis
    basisx(k) = forceMesh.basis(k).node(1);
    basisy(k) = forceMesh.basis(k).node(2);
end
nBasisx = size(unique(basisx),1);
nBasisy = size(unique(basisy),1);

% this is for assuring the boundary to be considered for being
% penalized for regularization
display('Building Laplacian Map...')
tic;
Lap = zeros(2*nBasis,2*nBasis);
k=1;
for ii=1:nBasisx
%     tempLapx = zeros(nBasisy,nBasis); %for parfor constraints
%     tempLapy = zeros(nBasisy,nBasis); %for parfor constraints
    for jj=1:nBasisy
        Lap2D = zeros(nBasisy,nBasisx);
        if ii==1 || jj==1 || ii==nBasisx || jj==nBasisy
            % this can be improved with basis function convolution
            if ii==1 
                if jj==1
                    Lap2D(jj,ii)=-2;
                    Lap2D(jj+1,ii)=1;
                    Lap2D(jj,ii+1)=1;
                elseif jj==nBasisy
                    Lap2D(jj,ii)=-2;
                    Lap2D(jj-1,ii)=1;
                    Lap2D(jj,ii+1)=1;
                else
                    Lap2D(jj,ii)=-2;
                    Lap2D(jj,ii+1)=1;
                    Lap2D(jj+1,ii)=1;
                end
            elseif ii==nBasisx
                if jj==1
                    Lap2D(jj,ii)=-2;
                    Lap2D(jj+1,ii)=1;
                    Lap2D(jj,ii-1)=1;
                elseif jj==nBasisy
                    Lap2D(jj,ii)=-2;
                    Lap2D(jj,ii-1)=1;
                    Lap2D(jj-1,ii)=1;
                else
                    Lap2D(jj,ii)=-2;
                    Lap2D(jj,ii-1)=1;
                    Lap2D(jj-1,ii)=1;
                end
            elseif jj==1
                if ii~=nBasisx && ii~=1
                    Lap2D(jj,ii)=-2;
                    Lap2D(jj+1,ii)=1;
                    Lap2D(jj,ii+1)=1;
                end
            elseif jj==nBasisy                        
                if ii~=nBasisx && ii~=1
                    Lap2D(jj,ii)=-2;
                    Lap2D(jj-1,ii)=1;
                    Lap2D(jj,ii+1)=1;
                end
            end                        
        else
            % diagonal laplacian
            Lap2D(jj,ii) = -6;
            Lap2D(jj,ii+1) = 1;
            Lap2D(jj,ii-1) = 1;
            Lap2D(jj+1,ii) = 1;
            Lap2D(jj-1,ii) = 1;
            Lap2D(jj+1,ii+1) = 0.5;
            Lap2D(jj-1,ii+1) = 0.5;
            Lap2D(jj+1,ii-1) = 0.5;
            Lap2D(jj-1,ii-1) = 0.5;
        end
        Lap(k,1:nBasis) = reshape(Lap2D,nBasis,1)';
        Lap(k+nBasis,nBasis+1:2*nBasis) = reshape(Lap2D,nBasis,1)';
%         Lap(k+nBasis,1:nBasis) = reshape(Lap2D,nBasis,1)';
%         Lap(k,nBasis+1:2*nBasis) = reshape(Lap2D,nBasis,1)';
        k=k+1;
%         tempLapx(jj,:) = reshape(Lap2D,nBasis,1)';
%         tempLapy(jj,:) = reshape(Lap2D,nBasis,1)';
    end
%     Lap((ii-1)*nBasisy+1:(ii-1)*nBasisy+nBasisy,1:nBasis) = tempLapx;
%     Lap((ii-1)*nBasisy+nBasis+1:(ii-1)*nBasisy+nBasis+nBasisy,nBasis+1:2*nBasis) = tempLapy;
end
toc

function [sol_coef,reg_corner] = calculateLfromLcurve(M,MpM,u,eyeWeights,LcurveDataPath,LcurveFigPath)
%examine a logarithmically spaced range of regularization parameters
alphas=10.^(-9:.125:-4.5);
rho=zeros(length(alphas),1);
eta=zeros(length(alphas),1);
mtik=zeros(size(M,2),length(alphas));
for i=1:length(alphas);
  mtik(:,i)=(MpM+alphas(i)*eyeWeights)\(M'*u);
  rho(i)=norm(M*mtik(:,i)-u);
  eta(i)=norm(mtik(:,i),1);
end

% Find the corner of the Tikhonov L-curve
% [reg_corner,ireg_corner,~]=l_curve_corner(rho,eta,alphas);
[reg_corner,ireg_corner,~]=regParamSelecetionLcurve(rho,eta,alphas);

% Plot the sparse deconvolution L-curve.
hLcurve = figure;
set(hLcurve, 'Position', [50 300 200 200])

loglog(rho,eta,'k-');
xlabel('Residual Norm ||Gm-d||_{2}');
ylabel('Solution Norm ||Lm||_{2}');
hold on
% mark and label the corner
if mod(ireg_corner,1)>0 % if ireg_corner is interpolated
    rho_corner = rho(floor(ireg_corner))+mod(ireg_corner,1)*(rho(floor(ireg_corner)+1)-rho(floor(ireg_corner)));
    eta_corner = eta(floor(ireg_corner))+mod(ireg_corner,1)*(eta(floor(ireg_corner)+1)-eta(floor(ireg_corner)));
else
    rho_corner = rho(ireg_corner);
    eta_corner = eta(ireg_corner);
end    
H=loglog(rho_corner,eta_corner,'ro');
set(H,'markersize',6)
H=text(rho_corner,1.1*eta_corner,...
    ['    ',num2str(reg_corner,'%5.1e')]);
set(H,'Fontsize',7);
% axis([1e-2 100 0.001 1e8])
disp('Printing L-curve...')
% print -deps2 nameSave
print(hLcurve,'Lcurve.eps','-depsc')
saveas(hLcurve,LcurveFigPath);

save(LcurveDataPath,'rho','eta','reg_corner','ireg_corner','alphas','mtik','-v7.3');

if mod(ireg_corner,1)>0 % if ireg_corner is interpolated
    disp(['L-corner regularization parmater L = ' num2str(reg_corner) '... final solution calculation ...'])
    sol_coef=(MpM+reg_corner*eyeWeights)\(M'*u);
else
    sol_coef = mtik(:,ireg_corner);
end

%old code
% display('Calculating L-curve ...')
% b  = M'*u;
% bfSigmaRange2 = L*[1e-5 1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4 1e5];
% bfSigmaRange1 = bfSigmaRange2/5;
% bfSigmaRange3 = bfSigmaRange2*5;
% bfSigmaCandidate = [bfSigmaRange1 bfSigmaRange2 bfSigmaRange3];
% bfSigmaCandidate = sort(bfSigmaCandidate);
% coefNorm    = zeros(1,length(bfSigmaCandidate));
% residueNorm = zeros(1,length(bfSigmaCandidate));
% Lap2 = Lap'*Lap;
% for kk = 1:length(bfSigmaCandidate)
%  B_i = MpM + bfSigmaCandidate(kk)*Lap2;
%  coef_i = B_i\b;
%  coefNorm(kk)    = sqrt(sum((Lap*coef_i).^2)); % Solution norm ||Lm||
%  residueNorm(kk) = sqrt(sum((M*coef_i-u).^2)); % Residual norm ||Gm-d||
% end
% 
% %Plot the L-curve.
% figure; hold on;
% plot(residueNorm,coefNorm,'.'); xlabel('Residual norm ||M*coef-u||'); 
% ylabel('Solution seminorm ||L*coef||');
% for kk = 1:length(bfSigmaCandidate)
%  text(residueNorm(kk),coefNorm(kk),num2str(bfSigmaCandidate(kk)));
% end
% 
% options.WindowStyle='normal';
% answer = inputdlg('Please identify the corner:','Input for corner',1,{num2str(L)},options);
% Lout = str2double(answer{1});

function [sol_coef,reg_corner] = calculateLfromLcurveSparse(L,M,MpM,u,eyeWeights,maxIter,tolx,tolr,nameSave,LcurveDataPath,LcurveFigPath,LcurveFactor)
%examine a logarithmically spaced range of regularization parameters
alphas=10.^(log10(L)-2.5:1.25/LcurveFactor:log10(L)+2);
rho=zeros(length(alphas),1);
eta=zeros(length(alphas),1);
msparse=zeros(size(M,2),length(alphas));
for i=1:length(alphas);
    disp(['testing L = ' num2str(alphas(i)) '... '])
    msparse(:,i)=iterativeL1Regularization(M,MpM,u,eyeWeights,alphas(i),maxIter,tolx,tolr);
    rho(i)=norm(M*msparse(:,i)-u);
    eta(i)=norm(msparse(:,i),1);
end

% Find the L-corner
% [reg_corner,ireg_corner,~]=l_curve_corner(rho,eta,alphas);
[reg_corner,ireg_corner,~]=regParamSelecetionLcurve(rho,eta,alphas);

% Plot the sparse deconvolution L-curve.
hLcurve = figure;
set(hLcurve, 'Position', [100 100 500 500])

loglog(rho,eta,'k-');
xlabel('Residual Norm ||Gm-d||_{2}');
ylabel('Solution Norm ||m||_{1}');
hold on
% mark and label the corner
if mod(ireg_corner,1)>0 % if ireg_corner is interpolated
    rho_corner = rho(floor(ireg_corner))+mod(ireg_corner,1)*(rho(floor(ireg_corner)+1)-rho(floor(ireg_corner)));
    eta_corner = eta(floor(ireg_corner))+mod(ireg_corner,1)*(eta(floor(ireg_corner)+1)-eta(floor(ireg_corner)));
else
    rho_corner = rho(ireg_corner);
    eta_corner = eta(ireg_corner);
end    
H=loglog(rho_corner,eta_corner,'ro');
set(H,'markersize',6)
H=text(rho_corner,1.1*eta_corner,...
    ['    ',num2str(reg_corner,'%5.1e')]);
set(H,'Fontsize',7);
% axis([1e-2 100 0.001 1e8])
disp('Displaying the 1-norm L-curve')
% print -deps2 nameSave
print(hLcurve,strcat(nameSave,'.eps'),'-depsc')
saveas(hLcurve,LcurveFigPath);
save(LcurveDataPath,'rho','eta','reg_corner','ireg_corner','alphas','rho_corner','eta_corner','msparse','-v7.3');

if mod(ireg_corner,1)>0 % if ireg_corner is interpolated
    disp(['L-corner regularization parmater L = ' num2str(reg_corner) '... final solution calculation ...'])
    sol_coef=iterativeL1Regularization(M,MpM,u,eyeWeights,reg_corner,maxIter,tolx,tolr);
else
    sol_coef = msparse(:,ireg_corner);
end


