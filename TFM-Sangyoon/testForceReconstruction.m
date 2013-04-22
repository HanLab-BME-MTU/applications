%% calculateForceFineMeshImgScaled_fn
function [] = testForceReconstruction(forceType,percentNoise,savePath,forceMeshFastBEMfine,M_FastBEMfine,basisClassTablePathFine)
% calculateForceFineMeshImgScaled_fn(forceType,percentNoise,savePath,forceMesh,M,M_FastBEM)
% This function compares the accuracy and performance of each force
% reconstruction technique using artificially given displacment with given
% noise level.
% example: testForceReconstruction('groupForce',1,
% '/home/sh268/orchestra/home/Documents/TFM-simulation/n1s',forceMesh,forceMeshFastBEM,M,M_FastBEM,
% '/home/sh268/orchestra/home/Documents/TFM-simulation/basisFunction5050.mat')
% % n1s means noise= 1 percent and s = smooth

% input:
%       imScale:        image scale from 50x50 pixel (10 means to increase
%                       image size to 500x500 pixels)
%       forceType:      path to the folder
%       percentNoise:   sec/frame
%       savePath:       image resolution, nm/pixel
%       forceMesh:      forceMesh
%       forceMeshFastBEM
%       M:              path and name of myosin image 
%       M_FastBEM:      total frame number
%       basisClassTablePath:     path to basisClassTable
% output:
%       images of forces (.fig and .tif):  stored in './img/'
%                   force maps, displacement field, x-sectional profile
%       data:       regularization parameter for each force method
%                   mean square deviation, peak force underestimation
%                   M and M_FastBEM
%                   stored in './data/'
% Sangyoon Han March 2013

%% parameter setup
E=8000;  %Young's modulus, unit: Pa
addNoise=1;
% percentNoise=1/100;
% savePath = 'xy0150fastBEM_10pNoise.mat';

s=0.5;  %Poisson's ratio, only needed for FTTC
imScale=1;
meshPtsFwdSol=2^10;
% L=0;
numPoints_u=50*imScale;       %(must be even number)
% numPoints_out=50;

xmin =1;
xmax =50*imScale;
ymin =1;
ymax =50*imScale;

imgPath = [savePath '/img/'];
dataPath = [savePath '/data/'];
if (~exist(imgPath,'dir') || ~exist(dataPath,'dir'))
    mkdir(imgPath);
    mkdir(dataPath);
end
cd(imgPath);
%% Mesh generation and artificial force generation
[x_mat_u, y_mat_u]=meshgrid(linspace(xmin,xmax,numPoints_u) , linspace(ymin,ymax,numPoints_u));
x_vec_u=reshape(x_mat_u,[],1);
y_vec_u=reshape(y_mat_u,[],1); % fine mesh for original force

% % Hexagonal force mesh
% displField.pos = [x_vec_u y_vec_u];
% [xvec,yvec]=createHexGridFromDisplField(displField,1);

% assumes um/pixel = 0.1, i.e. 10 pixel = 1 um.
% I need force distribution with multiple sources, unit: Pa -Sangyoon 013113
% I need to create 4 types of forces (large FA (20 pix), large force (2kPa) - group1
%                                     large FA, small force (500-1000 Pa) - group2
%                                     small FA (5 pix), large force - group3
%                                     small FA, small force - group4)
display('constructing original force field...')
tic;
fpos = [10*imScale,10*imScale;
        12*imScale,11*imScale;
        15*imScale,10*imScale;
        11*imScale,14*imScale; % 1st cluster
        
        15*imScale,26*imScale;
        17*imScale,31*imScale;
        28*imScale,27*imScale; % large adhesion
        21*imScale,20*imScale; % 2nd cluster
        
        27*imScale,39*imScale; 
        31*imScale,42*imScale; 
        37*imScale,43*imScale;
        31*imScale,36*imScale; %
        38*imScale,40*imScale;
        40*imScale,33*imScale; % 3rd cluster,  large adhesion

        33*imScale,9*imScale; 
        39*imScale,18*imScale; 
        43*imScale,24*imScale;
        38*imScale,12*imScale;
        42*imScale,17*imScale]; % 4th cluster
fvec = [1200,-100;
        1800,-160;
        1500,-130;
        1300,-110;
        
        1100,-380;
        1300,-480;
        850,-280; % large adhesion, small force
        400,-10;

        500,-600;
        450,-700;
        250,-900;
        800,-1200; 
        280,-800;
        700,-2000; % large adhesion, large force

        -450,500;
        -400,550;
        -350,400;
        -350,500;
        -150,600];

force_x1 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(1,1),fpos(1,2),fvec(1,1),fvec(1,2),forceType,'smallFA');
force_y1 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(1,1),fpos(1,2),fvec(1,1),fvec(1,2),forceType,'smallFA'); %group1
force_x2 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(2,1),fpos(2,2),fvec(2,1),fvec(2,2),forceType,'smallFA');
force_y2 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(2,1),fpos(2,2),fvec(2,1),fvec(2,2),forceType,'smallFA'); %group1
force_x3 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(3,1),fpos(3,2),fvec(3,1),fvec(3,2),forceType,'smallFA');
force_y3 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(3,1),fpos(3,2),fvec(3,1),fvec(3,2),forceType,'smallFA'); %group1
force_x4 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(4,1),fpos(4,2),fvec(4,1),fvec(4,2),forceType,'smallFA');
force_y4 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(4,1),fpos(4,2),fvec(4,1),fvec(4,2),forceType,'smallFA'); %group2
force_x5 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(5,1),fpos(5,2),fvec(5,1),fvec(5,2),forceType,'smallFA');
force_y5 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(5,1),fpos(5,2),fvec(5,1),fvec(5,2),forceType,'smallFA'); %group2
force_x6 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(6,1),fpos(6,2),fvec(6,1),fvec(6,2),forceType,'smallFA');
force_y6 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(6,1),fpos(6,2),fvec(6,1),fvec(6,2),forceType,'smallFA'); %group2
force_x7 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(7,1),fpos(7,2),fvec(7,1),fvec(7,2),forceType,'smallFA');
force_y7 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(7,1),fpos(7,2),fvec(7,1),fvec(7,2),forceType,'smallFA'); %group2, large one
force_x8 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(8,1),fpos(8,2),fvec(8,1),fvec(8,2),forceType,'smallFA');
force_y8 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(8,1),fpos(8,2),fvec(8,1),fvec(8,2),forceType,'smallFA'); %group2
force_x9 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(9,1),fpos(9,2),fvec(9,1),fvec(9,2),forceType,'smallFA');
force_y9 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(9,1),fpos(9,2),fvec(9,1),fvec(9,2),forceType,'smallFA'); %group2
force_x10 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(10,1),fpos(10,2),fvec(10,1),fvec(10,2),forceType,'smallFA');
force_y10 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(10,1),fpos(10,2),fvec(10,1),fvec(10,2),forceType,'smallFA'); %group1


force_x11 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(11,1),fpos(11,2),fvec(11,1),fvec(11,2),forceType,'smallFA');
force_y11 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(11,1),fpos(11,2),fvec(11,1),fvec(11,2),forceType,'smallFA'); %group1
force_x12 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(12,1),fpos(12,2),fvec(12,1),fvec(12,2),forceType,'smallFA');
force_y12 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(12,1),fpos(12,2),fvec(12,1),fvec(12,2),forceType,'smallFA'); %group1
force_x13 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(13,1),fpos(13,2),fvec(13,1),fvec(13,2),forceType,'smallFA');
force_y13 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(13,1),fpos(13,2),fvec(13,1),fvec(13,2),forceType,'smallFA'); %group1
force_x14 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(14,1),fpos(14,2),fvec(14,1),fvec(14,2),forceType,'smallFA');
force_y14 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(14,1),fpos(14,2),fvec(14,1),fvec(14,2),forceType,'smallFA'); %group2, large one
force_x15 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(15,1),fpos(15,2),fvec(15,1),fvec(15,2),forceType,'smallFA');
force_y15 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(15,1),fpos(15,2),fvec(15,1),fvec(15,2),forceType,'smallFA'); %group2
force_x16 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(16,1),fpos(16,2),fvec(16,1),fvec(16,2),forceType,'smallFA');
force_y16 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(16,1),fpos(16,2),fvec(16,1),fvec(16,2),forceType,'smallFA'); %group2
force_x17 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(17,1),fpos(17,2),fvec(17,1),fvec(17,2),forceType,'smallFA');
force_y17 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(17,1),fpos(17,2),fvec(17,1),fvec(17,2),forceType,'smallFA'); %group2
force_x18 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(18,1),fpos(18,2),fvec(18,1),fvec(18,2),forceType,'smallFA');
force_y18 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(18,1),fpos(18,2),fvec(18,1),fvec(18,2),forceType,'smallFA'); %group2
force_x19 = assumedForceShifted(1,x_mat_u,y_mat_u,fpos(19,1),fpos(19,2),fvec(19,1),fvec(19,2),forceType,'smallFA');
force_y19 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(19,1),fpos(19,2),fvec(19,1),fvec(19,2),forceType,'smallFA'); %group2

force_x = force_x1+force_x2+force_x3+force_x4+force_x5+force_x6+...
        force_x7+force_x8+force_x9+force_x10+force_x11+force_x12+...
        force_x13+force_x14+force_x15+force_x16+force_x17+force_x18+force_x19;
force_y = force_y1+force_y2+force_y3+force_y4+force_y5+force_y6+...
        force_y7+force_y8+force_y9+force_y10+force_y11+force_y12+...
        force_y13+force_y14+force_y15+force_y16+force_y17+force_y18+force_y19;
toc
%% Forward solution
display('forward solution for fine solution...')
tic
[ux, uy]=fwdSolution(x_mat_u,y_mat_u,E,[],[],[],[],...
    @(x,y) assumedForceShifted(1,x,y,10*imScale,10*imScale,1200,-100,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,12*imScale,11*imScale,1800,-160,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,15*imScale,10*imScale,1500,-130,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,11*imScale,14*imScale,1300,-110,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,15*imScale,26*imScale,1100,-380,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,17*imScale,31*imScale,1300,-480,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,28*imScale,27*imScale,850,-280,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,21*imScale,20*imScale,400,-10,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,27*imScale,39*imScale,500,-600,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,31*imScale,42*imScale,450,-700,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,37*imScale,43*imScale,250,-900,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,31*imScale,36*imScale,800,-1200,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,38*imScale,40*imScale,280,-800,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,40*imScale,33*imScale,700,-1000,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,33*imScale,9*imScale,-450,-800,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,39*imScale,18*imScale,-400,550,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,43*imScale,24*imScale,-350,400,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,38*imScale,12*imScale,-350,500,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,42*imScale,17*imScale,-150,600,forceType,'smallFA'),...
    @(x,y) assumedForceShifted(2,x,y,10*imScale,10*imScale,1200,-100,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,12*imScale,11*imScale,1800,-160,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,15*imScale,10*imScale,1500,-130,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,11*imScale,14*imScale,1300,-110,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,15*imScale,26*imScale,1100,-380,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,17*imScale,31*imScale,1300,-480,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,28*imScale,27*imScale,850,-280,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,21*imScale,20*imScale,400,-10,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,27*imScale,39*imScale,500,-600,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,31*imScale,42*imScale,450,-700,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,37*imScale,43*imScale,250,-900,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,31*imScale,36*imScale,800,-1200,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,38*imScale,40*imScale,280,-800,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,40*imScale,33*imScale,700,-1000,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,33*imScale,9*imScale,-450,-800,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,39*imScale,18*imScale,-400,550,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,43*imScale,24*imScale,-350,400,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,38*imScale,12*imScale,-350,500,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,42*imScale,17*imScale,-150,600,forceType,'smallFA'),'fft',[],meshPtsFwdSol);
toc

display(['adding noise (' num2str(percentNoise*100) '%) ...'])
if addNoise==1
    max_u=max([max(max(ux)) max(max(uy))]); 
    noise_xu=normrnd(0,percentNoise*max_u,numPoints_u,numPoints_u);
    noise_yu=normrnd(0,percentNoise*max_u,numPoints_u,numPoints_u);
    ux=ux+noise_xu;
    uy=uy+noise_yu;
end

u_mat(:,:,1)=ux;
u_mat(:,:,2)=uy;
ux_vec=reshape(ux,[],1);
uy_vec=reshape(uy,[],1);
% u=vertcat(ux_vec,uy_vec);

%% grid setting
pix_durch_my=1; %seem to be important only for the energy
% grid_mat_u(:,:,1)=x_mat_u;
% grid_mat_u(:,:,2)=y_mat_u;
% i_max = size(grid_mat_u,1);
% j_max = size(grid_mat_u,2);
% cluster_size = grid_mat_u(1,1,1) - grid_mat_u(2,2,1);

% % [pos,vec,force_FTTC, fnorm_FTTC] = reg_fourier_TFM_used_till_2010_07_16(grid_mat_u,u_mat,E,s, pix_durch_my, cluster_size, i_max, j_max, 0);  %0.000020
% [~,~,force_FTTC, ~] = reg_fourier_TFM(grid_mat_u,u_mat,E,s, pix_durch_my, cluster_size, i_max, j_max, 0);  %0.000020
% 
% rec_force_FTTC(:,:,1)=reshape(force_FTTC(:,1),i_max,j_max);
% rec_force_FTTC(:,:,2)=reshape(force_FTTC(:,2),i_max,j_max);
%% FTTC-reconstruction with regularization with fine mesh
L= 8.8e-8;
grid_mat_u(:,:,1)=x_mat_u; %dense
grid_mat_u(:,:,2)=y_mat_u; %dense
i_max = size(grid_mat_u,1);
j_max = size(grid_mat_u,2);
cluster_size = grid_mat_u(1,1,1) - grid_mat_u(2,2,1);

[~,~,force_FTTCfine, ~] = reg_fourier_TFM(grid_mat_u,u_mat,E,s, pix_durch_my, cluster_size, i_max, j_max, L);  
fttcL = L;
rec_force_FTTCregfine(:,:,1)=reshape(force_FTTCfine(:,1),i_max,j_max);
rec_force_FTTCregfine(:,:,2)=reshape(force_FTTCfine(:,2),i_max,j_max);
%% FastBEM with fine mesh (For figure 3) 
% Mesh creation
if isempty(forceMeshFastBEMfine)
    display('creating fine mesh and basis function (1x1)')
    tic
    forceMeshFastBEMfine=createMeshAndBasisFastBEM(x_vec_u,y_vec_u,true,[],0);
    toc
end
L = 8.8e-8;
display(['FastBEM with fine mesh (' num2str(L) ')...'])
[fx_FastBEMregfine,fy_FastBEMregfine,~,~,M_FastBEMfine,~,~,~,sol_mats_regfine]...
    = BEM_force_reconstruction(x_mat_u,y_mat_u,ux,uy,forceMeshFastBEMfine,...
    E,L,[],[],'fast',meshPtsFwdSol,...
    'backslash','basisClassTblPath',basisClassTablePathFine,'imgRows',ymax,'imgCols',xmax,'fwdMap',M_FastBEMfine);
rec_force_FastBEMregfine(:,:,1)=fx_FastBEMregfine;
rec_force_FastBEMregfine(:,:,2)=fy_FastBEMregfine;
%% FastBEM-reconstruction with laplacian regularization with fine mesh
L = 7.1e-9; %sol_mats_reg.L;
display(['FastBEM with 2nd order regularization (' num2str(L) ')...'])
tic
[fx_FastBEMreg2fine,fy_FastBEMreg2fine,~,~,~,~,~,~,sol_mats_reg2fine]...
    = BEM_force_reconstruction(x_mat_u,y_mat_u,ux,uy,forceMeshFastBEMfine,...
    E,L,[],[],'fast',meshPtsFwdSol,...
    'LaplacianReg','basisClassTblPath',basisClassTablePathFine,...
    'imgRows',ymax,'imgCols',xmax,'fwdMap',M_FastBEMfine);
toc
rec_force_FastBEMreg2fine(:,:,1)=fx_FastBEMreg2fine;
rec_force_FastBEMreg2fine(:,:,2)=fy_FastBEMreg2fine;

%% FastBEM-reconstruction with 1Norm-0th order regularization with fine mesh
L = 4e-5; %sol_mats_reg.L;
display(['FastBEM with 1Norm-0th order regularization with fine mesh (' num2str(L) ')...'])
tic
[fx_FastBEM1nfine,fy_FastBEM1nfine,~,~,~,~,~,~,sol_mats_1nfine]...
    = BEM_force_reconstruction(x_mat_u,y_mat_u,ux,uy,forceMeshFastBEMfine,...
    E,L,[],[],'fast',meshPtsFwdSol,...
    '1NormReg','basisClassTblPath',basisClassTablePathFine,...
    'imgRows',ymax,'imgCols',xmax,'fwdMap',M_FastBEMfine);
toc
rec_force_FastBEM1nfine(:,:,1)=fx_FastBEM1nfine;
rec_force_FastBEM1nfine(:,:,2)=fy_FastBEM1nfine;

%% FastBEM-reconstruction with 1Norm-2th order regularization with fine mesh
L = 1e-4; %sol_mats_reg.L;
display(['FastBEM with 1Norm-2nd order regularization with fine mesh (' num2str(L) ')...'])
tic
[fx_FastBEM1n2fine,fy_FastBEM1n2fine,~,~,~,~,~,~,sol_mats_1n2fine]...
    = BEM_force_reconstruction(x_mat_u,y_mat_u,ux,uy,forceMeshFastBEMfine,...
    E,L,[],[],'fast',meshPtsFwdSol,...
    '1NormRegLaplacian','basisClassTblPath',basisClassTablePathFine,...
    'imgRows',ymax,'imgCols',xmax,'fwdMap',M_FastBEMfine);
toc
rec_force_FastBEM1n2fine(:,:,1)=fx_FastBEM1n2fine;
rec_force_FastBEM1n2fine(:,:,2)=fy_FastBEM1n2fine;

%% quick saving of M and forcemesh
save([dataPath '/MandforceMesh.mat'],'M_FastBEMfine','forceMeshFastBEMfine',...
    'rec_force_FTTCregfine','rec_force_FastBEMregfine','rec_force_FastBEMreg2fine','fy_FastBEMreg2fine','rec_force_FastBEM1nfine',...
    'rec_force_FastBEM1n2fine','-v7.3');
%% Heat map presentation for each force magnitude
pos = [x_vec_u y_vec_u]; %dense
disp_vec = [ux_vec uy_vec]; %dense

rec_force_org_vec = [reshape(force_x,[],1) reshape(force_y,[],1)];
[grid_matfine,tmat_org, ~, ~] = interp_vec2grid(pos+disp_vec, rec_force_org_vec,[],grid_mat_u); %1:cluster size
tnorm_org = (tmat_org(:,:,1).^2 + tmat_org(:,:,2).^2).^0.5; %this should be fine mesh
fmat_org(:,:,1) = force_x;
fmat_org(:,:,2) = force_y;
fnorm_org = (force_x.^2 + force_y.^2).^0.5; %this should be fine mesh

rec_force_FTTCregfine_vec = [reshape(rec_force_FTTCregfine(:,:,1),[],1) reshape(rec_force_FTTCregfine(:,:,2),[],1)];
[~,tmat_fttcregfine, ~, ~] = interp_vec2grid(pos+disp_vec, rec_force_FTTCregfine_vec,[],grid_mat_u); %1:cluster size
tnorm_fttcregfine = (tmat_fttcregfine(:,:,1).^2 + tmat_fttcregfine(:,:,2).^2).^0.5;
fnorm_fttcregfine = (rec_force_FTTCregfine(:,:,1).^2 + rec_force_FTTCregfine(:,:,2).^2).^0.5;

rec_force_FastBEMregfine_vec = [reshape(rec_force_FastBEMregfine(:,:,1),[],1) reshape(rec_force_FastBEMregfine(:,:,2),[],1)];
[~,tmat_FastBEMregfine, ~, ~] = interp_vec2grid(pos+disp_vec, rec_force_FastBEMregfine_vec,[],grid_mat_u); %1:cluster size
tnorm_FastBEMregfine = (tmat_FastBEMregfine(:,:,1).^2 + tmat_FastBEMregfine(:,:,2).^2).^0.5;
fmat_FastBEMregfine(:,:,1) = reshape(rec_force_FastBEMregfine(:,:,1),size(x_mat_u));
fmat_FastBEMregfine(:,:,2) = reshape(rec_force_FastBEMregfine(:,:,2),size(y_mat_u));
fnorm_FastBEMregfine = (fmat_FastBEMregfine(:,:,1).^2 + fmat_FastBEMregfine(:,:,2).^2).^0.5;

rec_force_FastBEMreg2fine_vec = [reshape(rec_force_FastBEMreg2fine(:,:,1),[],1) reshape(rec_force_FastBEMreg2fine(:,:,2),[],1)];
[~,tmat_FastBEMreg2fine, ~, ~] = interp_vec2grid(pos+disp_vec, rec_force_FastBEMreg2fine_vec,[],grid_mat_u); %1:cluster size
tnorm_FastBEMreg2fine = (tmat_FastBEMreg2fine(:,:,1).^2 + tmat_FastBEMreg2fine(:,:,2).^2).^0.5;
fmat_FastBEMreg2fine(:,:,1) = reshape(rec_force_FastBEMreg2fine(:,:,1),size(x_mat_u));
fmat_FastBEMreg2fine(:,:,2) = reshape(rec_force_FastBEMreg2fine(:,:,2),size(y_mat_u));
fnorm_FastBEMreg2fine = (fmat_FastBEMreg2fine(:,:,1).^2 + fmat_FastBEMreg2fine(:,:,2).^2).^0.5;

rec_force_FastBEM1nfine_vec = [reshape(rec_force_FastBEM1nfine(:,:,1),[],1) reshape(rec_force_FastBEM1nfine(:,:,2),[],1)];
[~,tmat_FastBEM1nfine, ~, ~] = interp_vec2grid(pos+disp_vec, rec_force_FastBEM1nfine_vec,[],grid_mat_u); %1:cluster size
tnorm_FastBEM1nfine = (tmat_FastBEM1nfine(:,:,1).^2 + tmat_FastBEM1nfine(:,:,2).^2).^0.5;
fmat_FastBEM1nfine(:,:,1) = reshape(rec_force_FastBEM1nfine(:,:,1),size(x_mat_u));
fmat_FastBEM1nfine(:,:,2) = reshape(rec_force_FastBEM1nfine(:,:,2),size(y_mat_u));
fnorm_FastBEM1nfine = (fmat_FastBEM1nfine(:,:,1).^2 + fmat_FastBEM1nfine(:,:,2).^2).^0.5;

rec_force_FastBEM1n2fine_vec = [reshape(rec_force_FastBEM1n2fine(:,:,1),[],1) reshape(rec_force_FastBEM1n2fine(:,:,2),[],1)];
[~,tmat_FastBEM1n2fine, ~, ~] = interp_vec2grid(pos+disp_vec, rec_force_FastBEM1n2fine_vec,[],grid_mat_u); %1:cluster size
tnorm_FastBEM1n2fine = (tmat_FastBEM1n2fine(:,:,1).^2 + tmat_FastBEM1n2fine(:,:,2).^2).^0.5;
fmat_FastBEM1n2fine(:,:,1) = reshape(rec_force_FastBEM1n2fine(:,:,1),size(x_mat_u));
fmat_FastBEM1n2fine(:,:,2) = reshape(rec_force_FastBEM1n2fine(:,:,2),size(y_mat_u));
fnorm_FastBEM1n2fine = (fmat_FastBEM1n2fine(:,:,1).^2 + fmat_FastBEM1n2fine(:,:,2).^2).^0.5;

%% Force at cross-sectional line for large FA, high peak
% For force 1
tnorm1D_org_g1 = fnorm_org(fpos(14,2),fpos(14,1)-10*imScale:fpos(14,1)+10*imScale);
tnorm1D_FTTCregfine_g1 = fnorm_fttcregfine(fpos(14,2),fpos(14,1)-10*imScale:fpos(14,1)+10*imScale);
tnorm1D_FastBEMregfine_g1 = fnorm_FastBEMregfine(fpos(14,2),fpos(14,1)-10*imScale:fpos(14,1)+10*imScale);
tnorm1D_FastBEMreg2fine_g1 = fnorm_FastBEMreg2fine(fpos(14,2),fpos(14,1)-10*imScale:fpos(14,1)+10*imScale);
tnorm1D_FastBEM1nfine_g1 = fnorm_FastBEM1nfine(fpos(14,2),fpos(14,1)-10*imScale:fpos(14,1)+10*imScale);
tnorm1D_FastBEM1n2fine_g1 = fnorm_FastBEM1n2fine(fpos(14,2),fpos(14,1)-10*imScale:fpos(14,1)+10*imScale);
x1Du_g1 = grid_mat_u(fpos(14,2),fpos(14,1)-10*imScale:fpos(14,1)+10*imScale,1);
% %% Force at cross-sectional line at force_4 for large FA, low peak
% % force_x4 = assumedForceShifted(1,x_mat_u,y_mat_u,37*10,45*10,300,-800,forceType,'largeFA');
% % force_y4 = assumedForceShifted(2,x_mat_u,y_mat_u,37*10,45*10,300,-800,forceType,'largeFA'); %group2
% tnorm1D_org_g2 = fnorm_org(fpos(7,2),fpos(7,1)-10*imScale:fpos(7,1)+10*imScale);
% tnorm1D_FTTCregfine_g2 = fnorm_fttcregfine(fpos(7,2),fpos(7,1)-10*imScale:fpos(7,1)+10*imScale);
% tnorm1D_FastBEMregfine_g2 = fnorm_FastBEMregfine(fpos(7,2),fpos(7,1)-10*imScale:fpos(7,1)+10*imScale);
% tnorm1D_FastBEMreg2fine_g2 = fnorm_FastBEMreg2fine(fpos(7,2),fpos(7,1)-10*imScale:fpos(7,1)+10*imScale);
% tnorm1D_FastBEMn1fine_g2 = fnorm_FastBEM1nfine(fpos(7,2),fpos(7,1)-10*imScale:fpos(7,1)+10*imScale);
% x1Du_g2 = grid_mat_u(fpos(7,2),fpos(7,1)-10*imScale:fpos(7,1)+10*imScale,1);
% %% Force at cross-sectional line at force_7 for small FA, high peak
% tnorm1D_org_g3 = fnorm_org(31*imScale,17*imScale-10*imScale:17*imScale+10*imScale);
% tnorm1D_FTTCregfine_g3 = fnorm_fttcregfine(31*imScale,17*imScale-10*imScale:17*imScale+10*imScale);
% tnorm1D_FastBEMregfine_g3 = fnorm_FastBEMregfine(31*imScale,17*imScale-10*imScale:17*imScale+10*imScale);
% tnorm1D_FastBEMreg2fine_g3 = fnorm_FastBEMreg2fine(31*imScale,17*imScale-10*imScale:17*imScale+10*imScale);
% tnorm1D_FastBEMn1fine_g3 = fnorm_FastBEM1nfine(31*imScale,17*imScale-10*imScale:17*imScale+10*imScale);
% x1Du_g3 = grid_mat_u(31*imScale,17*imScale-10*imScale:17*imScale+10*imScale,1);
% %% Force at cross-sectional line at force_10 for small FA, low peak
% % force_x10 = assumedForceShifted(1,x_mat_u,y_mat_u,12*10,20*10,600,-480,forceType,'smallFA');
% % force_y10 = assumedForceShifted(2,x_mat_u,y_mat_u,12*10,20*10,600,-480,forceType,'smallFA'); %group4
% tnorm1D_org_g4 = fnorm_org(20*imScale,21*imScale-10*imScale:21*imScale+10*imScale);
% tnorm1D_FTTCregfine_g4 = fnorm_fttcregfine(20*imScale,21*imScale-10*imScale:21*imScale+10*imScale);
% tnorm1D_FastBEMregfine_g4 = fnorm_FastBEMregfine(20*imScale,21*imScale-10*imScale:21*imScale+10*imScale);
% tnorm1D_FastBEMreg2fine_g4 = fnorm_FastBEMreg2fine(20*imScale,21*imScale-10*imScale:21*imScale+10*imScale);
% tnorm1D_FastBEMn1fine_g4 = fnorm_FastBEM1nfine(20*imScale,21*imScale-10*imScale:21*imScale+10*imScale);
% x1Du_g4 = grid_mat_u(20*imScale,21*imScale-10*imScale:21*imScale+10*imScale,1);
%% Boundary cutting for FastBEM
% tnorm_FastBEMregsc(1,:)=0;
% tnorm_FastBEMregsc(:,1)=0;
% tnorm_FastBEMregsc(end,:)=0;
% tnorm_FastBEMregsc(:,end)=0;
% tnorm_FastBEMreg2(1,:)=0;
% tnorm_FastBEMreg2(:,1)=0;
% tnorm_FastBEMreg2(end,:)=0;
% tnorm_FastBEMreg2(:,end)=0;

%% Precision in locating peak force
tnormd(:,:,1) = fnorm_fttcregfine;
tnormd(:,:,2) = fnorm_FastBEMregfine;
tnormd(:,:,3) = fnorm_FastBEMreg2fine;
tnormd(:,:,4) = fnorm_FastBEM1nfine;
tnormd(:,:,5) = fnorm_FastBEM1n2fine;
nMethods = 5;

% zeroI = [1 1]; %index of zero coordinates
% [~,locMaxI,~] = findMaxScoreI(tnorm_org,zeroI,1,0); %from dense data
locMaxI = fpos(:,2:-1:1);
nPeaks = size(locMaxI,1);
flocMax = zeros(nPeaks, nMethods);
flocMaxRatio = zeros(nPeaks, nMethods);
flocMaxAngle = zeros(nPeaks, nMethods);
flocMaxDTMS = zeros(nPeaks, nMethods);
pSR = zeros(nMethods,1); %peak stress ratio
pSRerr = zeros(nMethods,1);
DTA = zeros(nMethods,1); %deviation of traction angle
DTAerr = zeros(nMethods,1);
DTMS = zeros(nMethods,1); %deviation of traction magnitude in the surrounding
DTMSerr = zeros(nMethods,1);
flocMaxOrg = diag(fnorm_org(locMaxI(:,1),locMaxI(:,2)));
r = 1.1; %
backgroundMaskfine = sqrt((x_mat_u-(locMaxI(1,1))).^2+(y_mat_u-(locMaxI(1,2))).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(2,1))).^2+(y_mat_u-(locMaxI(2,2))).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(3,1))).^2+(y_mat_u-(locMaxI(3,2))).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(4,1))).^2+(y_mat_u-(locMaxI(4,2))).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(5,1))).^2+(y_mat_u-(locMaxI(5,2))).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(6,1)/5+1)).^2+(y_mat_u-(locMaxI(6,2)/5+1)).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(7,1)/5+1)).^2+(y_mat_u-(locMaxI(7,2)/5+1)).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(8,1)/5+1)).^2+(y_mat_u-(locMaxI(8,2)/5+1)).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(9,1)/5+1)).^2+(y_mat_u-(locMaxI(9,2)/5+1)).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(10,1)/5+1)).^2+(y_mat_u-(locMaxI(10,2)/5+1)).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(11,1)/5+1)).^2+(y_mat_u-(locMaxI(11,2)/5+1)).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(12,1)/5+1)).^2+(y_mat_u-(locMaxI(12,2)/5+1)).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(13,1)/5+1)).^2+(y_mat_u-(locMaxI(13,2)/5+1)).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(14,1)/5+1)).^2+(y_mat_u-(locMaxI(14,2)/5+1)).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(15,1)/5+1)).^2+(y_mat_u-(locMaxI(15,2)/5+1)).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(16,1)/5+1)).^2+(y_mat_u-(locMaxI(16,2)/5+1)).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(17,1)/5+1)).^2+(y_mat_u-(locMaxI(17,2)/5+1)).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(18,1)/5+1)).^2+(y_mat_u-(locMaxI(18,2)/5+1)).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(19,1)/5+1)).^2+(y_mat_u-(locMaxI(19,2)/5+1)).^2)>r;

for jj=1:nMethods 
    for k=1:nPeaks
        flocMax(k,jj) = tnormd(locMaxI(k,1),locMaxI(k,2),jj);
        flocMaxRatio(k,jj) = flocMax(k,jj)/fnorm_org(locMaxI(k,1),locMaxI(k,2)); % ratio based on original force norm (0~1)
        if jj==1
            recF_vecfine = [rec_force_FTTCregfine(locMaxI(k,1),locMaxI(k,2),1),...
                rec_force_FTTCregfine(locMaxI(k,1),locMaxI(k,2),2)];
        elseif jj==2
            recF_vecfine = [fmat_FastBEMregfine(locMaxI(k,1)+1,locMaxI(k,2),1),...
                fmat_FastBEMregfine(locMaxI(k,1),locMaxI(k,2),2)];
        elseif jj==3
            recF_vecfine = [fmat_FastBEMreg2fine(locMaxI(k,1),locMaxI(k,2),1),...
                fmat_FastBEMreg2fine(locMaxI(k,1),locMaxI(k,2),2)];
        elseif jj==4
            recF_vecfine = [fmat_FastBEM1nfine(locMaxI(k,1)+1,locMaxI(k,2),1),...
                fmat_FastBEM1nfine(locMaxI(k,1),locMaxI(k,2),2)];
        elseif jj==5
            recF_vecfine = [fmat_FastBEM1n2fine(locMaxI(k,1),locMaxI(k,2),1),...
                fmat_FastBEM1n2fine(locMaxI(k,1),locMaxI(k,2),2)];
        end
        orgF_vec = [fmat_org(locMaxI(k,1),locMaxI(k,2),1),...
                    fmat_org(locMaxI(k,1),locMaxI(k,2),2)];

        flocMaxAngle(k,jj) = acosd(recF_vecfine*orgF_vec'/...
            (norm(recF_vecfine)*norm(orgF_vec))); % in degree
        %indices for surroundings (2 um radius)
        if k==7 || k==14
            r = 1.1;%4.1;
        else
            r = 1.1;
        end
        ringMaskfine = sqrt((x_mat_u-locMaxI(k,1)).^2+(y_mat_u-locMaxI(k,2)).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(k,1))).^2+(y_mat_u-(locMaxI(k,2))).^2)<=r+1.1 ...
                    & backgroundMaskfine;
        neighborTnormfine = ringMaskfine.*tnormd(:,:,jj);
        meanNeiTnormfine = sum(sum(neighborTnormfine))/sum(sum(ringMaskfine));
        flocMaxDTMS(k,jj) = meanNeiTnormfine/flocMax(k,jj);
        % Peak force localization(position) match - skipped for now.
    end
    % statistic of flocMaxRatio
    pSR(jj) = mean(flocMaxRatio(:,jj));
    pSRerr(jj) = sqrt(std(flocMaxRatio(:,jj))/nPeaks);
    DTA(jj) = mean(flocMaxAngle(:,jj));
    DTAerr(jj) = sqrt(std(flocMaxAngle(:,jj))/nPeaks);
    DTMS(jj) = mean(flocMaxDTMS(:,jj));
    DTMSerr(jj) = sqrt(std(flocMaxDTMS(:,jj))/nPeaks);
end
%% Figures
%% Quiver plots of BEM results compared to original, FTTC, regularized FTTC 
h5 = figure;
forceScale=.5*sqrt(max(max(force_x1.^2+force_y1.^2)));
hold on
quiver(x_mat_u,y_mat_u,force_x/forceScale,force_y/forceScale,0,'k');
quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),rec_force_FTTCregfine(:,:,1)/forceScale,rec_force_FTTCregfine(:,:,2)/forceScale,0,'g');
quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),fmat_FastBEMregfine(:,:,1)/forceScale,fmat_FastBEMregfine(:,:,2)/forceScale,0,'b');
quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),fmat_FastBEMreg2fine(:,:,1)/forceScale,fmat_FastBEMreg2fine(:,:,2)/forceScale,0,'r');
quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),fmat_FastBEM1nfine(:,:,1)/forceScale,fmat_FastBEM1nfine(:,:,2)/forceScale,0,'c');
quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),fmat_FastBEM1n2fine(:,:,1)/forceScale,fmat_FastBEM1n2fine(:,:,2)/forceScale,0,'m');

hold off
xlim([xmin-1 xmax+1])
ylim([ymin-1 ymax+1])
% Make the text of the legend italic and color it brown
hleg = legend('Org Force','FTTC',...
    'L2Norm 0th','L2Norm 2nd',...
    'L1Norm 0th','L1Norm 2nd',...
              'Location','NorthEast');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1],'Fontsize',6)
title('reconstructed forces')
hgexport(h5,strcat(imgPath,'force vector field'),hgexport('factorystyle'),'Format','tiff')
hgsave(h5,strcat(imgPath,'force vector field'),'-v7.3')
% delete(h5)

%% Heatmaps
h1 = figure;
hold off
set(h1, 'Position', [100 100 1200 800])
%% Heatmaps
tmin = min(min(tnorm_org))-100;
tmax = max(max(tnorm_org))*1.1;
colormap('jet');
subplot(3,4,1), surf(grid_matfine(:,:,1), grid_matfine(:,:,2), tnorm_org,'FaceColor','interp',...
	'EdgeColor','none', 'FaceLighting','phong'), zlim([tmin-1 tmax+1]), ...
    view(0,90), title('Original force')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
% % hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
xlim([xmin xmax])
ylim([ymin ymax])

subplot(3,4,2), surf(grid_matfine(:,:,1), grid_matfine(:,:,2), tnorm_fttcregfine), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title(['FTTC regularization L=' num2str(fttcL)])
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
% hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
xlim([xmin xmax])
ylim([ymin ymax])

subplot(3,4,3), surf(grid_matfine(:,:,1), grid_matfine(:,:,2), tnorm_FastBEMregfine), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title(['FastBEM regularization L=' num2str(sol_mats_regfine.L)])
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
% hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
xlim([xmin xmax])
ylim([ymin ymax])

subplot(3,4,5), surf(grid_matfine(:,:,1), grid_matfine(:,:,2), tnorm_FastBEMreg2fine), zlim([tmin-1 tmax+1])
view(0,90), shading interp; title(['FastBEM 2nd order regularization L=' num2str(sol_mats_reg2fine.L)] )
set(gca, 'DataAspectRatio', [1,1,1]);
% hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
xlim([xmin xmax])
ylim([ymin ymax])
caxis([tmin-0.01 tmax+0.01]);

subplot(3,4,6), surf(grid_matfine(:,:,1), grid_matfine(:,:,2), tnorm_FastBEM1nfine), zlim([tmin-1 tmax+1])
view(0,90), shading interp; title(['FastBEM 1norm regularization L=' num2str(sol_mats_1nfine.L)])
set(gca, 'DataAspectRatio', [1,1,1]);
% hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
xlim([xmin xmax])
ylim([ymin ymax])
caxis([tmin-0.01 tmax+0.01]);

subplot(3,4,7), surf(grid_matfine(:,:,1), grid_matfine(:,:,2), tnorm_FastBEM1n2fine), zlim([tmin-1 tmax+1])
view(0,90), shading interp; title(['FastBEM 1norm Laplacian L=' num2str(sol_mats_1n2fine.L)])
set(gca, 'DataAspectRatio', [1,1,1]);
% hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
xlim([xmin xmax])
ylim([ymin ymax])
caxis([tmin-0.01 tmax+0.01]);

subplot(3,4,8), caxis([tmin-0.01 tmax+0.01]), axis off; colorbar('West')
%% profile plot
subplot(3,4,9), %for high peak
plot(x1Du_g1,tnorm1D_org_g1,'k','Linewidth',2)
hold on
plot(x1Du_g1,tnorm1D_FTTCregfine_g1,'c','Linewidth',2)
plot(x1Du_g1,tnorm1D_FastBEMregfine_g1,'g','Linewidth',2)
plot(x1Du_g1,tnorm1D_FastBEMreg2fine_g1,'r','Linewidth',2)
plot(x1Du_g1,tnorm1D_FastBEM1nfine_g1,'b','Linewidth',2)
plot(x1Du_g1,tnorm1D_FastBEM1n2fine_g1,'m','Linewidth',2)
xlabel('x')
ylabel('Traction Stress (Pa)')
ymax1 = max(tnorm1D_org_g1)*1.5;
ylim([0 ymax1]);
xlim([x1Du_g1(1) x1Du_g1(end)]);
title({'1D stress profile for','large FA, high peak'})

% subplot(3,4,12), %for high peak
% plot(x1Du_g2,tnorm1D_org_g2,'k','Linewidth',2)
% hold on
% plot(x1Du_g2,tnorm1D_FTTCregfine_g2,'g','Linewidth',2)
% plot(x1Du_g2,tnorm1D_FastBEMregfine_g2,'b','Linewidth',2)
% plot(x1Du_g2,tnorm1D_FastBEMreg2fine_g2,'r','Linewidth',2)
% plot(x1Du_g2,tnorm1D_FastBEMn1fine_g2,'m','Linewidth',2)
% xlabel('x')
% ylabel('Traction Stress (Pa)')
% ylim([0 ymax1]);
% title({'1D stress profile for','large FA, low peak'})
% 
% subplot(3,4,13), %for high peak
% plot(x1Du_g3,tnorm1D_org_g3,'k','Linewidth',2)
% hold on
% plot(x1Du_g3,tnorm1D_FTTCregfine_g3,'g','Linewidth',2)
% plot(x1Du_g3,tnorm1D_FastBEMregfine_g3,'b','Linewidth',2)
% plot(x1Du_g3,tnorm1D_FastBEMreg2fine_g3,'r','Linewidth',2)
% plot(x1Du_g3,tnorm1D_FastBEMn1fine_g3,'m','Linewidth',2)
% xlabel('x')
% ylabel('Traction Stress (Pa)')
% ylim([0 ymax1]);
% title({'1D stress profile for','small FA, high peak'})
% 
% subplot(3,4,14), %for low peak
% plot(x1Du_g4,tnorm1D_org_g4,'k','Linewidth',2)
% hold on
% plot(x1Du_g4,tnorm1D_FTTCregfine_g4,'g','Linewidth',2)
% plot(x1Du_g4,tnorm1D_FastBEMregfine_g4,'b','Linewidth',2)
% plot(x1Du_g4,tnorm1D_FastBEMreg2fine_g4,'r','Linewidth',2)
% plot(x1Du_g4,tnorm1D_FastBEMn1fine_g4,'m','Linewidth',2)
% xlabel('x')
% ylabel('Traction Stress (Pa)')
% ylim([0 ymax1]);
% title({'1D stress profile for','small FA, low peak'})
hleg = legend('Original','FTTC','L2 0th','L2 2nd','L1 0th','L1 2nd','Location','East');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1],'FontSize',6)

%% Peak Magnitude Ratio
subplot(3,4,10),
bar(pSR), hold on
errorbar(pSR,pSRerr)
title('Peak Stress ratio')
xlabel({'FTTC reg vs. ',' FastBEM reg vs. ',' FastBEM 2nd reg',' FastBEM 1norm',' 1norm 2nd'});
%% Deviation of traction magnitude surroundings (DTMS)
subplot(3,4,11),
bar(DTMS), hold on
errorbar(DTMS,DTMSerr)
title('Deviation of Traction Magnitude from Surroundings')
xlabel({'FTTC reg vs. ',' FastBEM reg vs. ',' FastBEM 2nd reg',' FastBEM 1norm',' 1norm 2nd'});

%% save
hgexport(h1,strcat(imgPath,'fig3 TFM accuracy with finer mesh'),hgexport('factorystyle'),'Format','tiff')
hgsave(h1,strcat(imgPath,'fig3 TFM accuracy with finer mesh'),'-v7.3')
% delete(h1)

%% Fig S1
h2 = figure;
%% Peak plot
subplot(1,2,1)
plot(flocMaxOrg/1000,flocMax(:,1)/1000,'g.'), hold on
plot(flocMaxOrg/1000,flocMax(:,2)/1000,'b.'), 
plot(flocMaxOrg/1000,flocMax(:,3)/1000,'r.'), 
plot(flocMaxOrg/1000,flocMax(:,4)/1000,'go'), 
plot(flocMaxOrg/1000,flocMax(:,5)/1000,'bo'), 
plot(0:.1:3.2,0:.1:3.2,'k--'),hold off
xlim([0 2.7]),ylim([0 3.3]);
xlabel('Original Peak Stress (kPa)')
ylabel('Reconstructed Peak Stress (kPa)')
title({'Reconstructed Peak Stress','w.r.t. Original Peak Stress'})
hleg2=legend('FTTCf','0thf','2ndf','1norm','1norm2nd','Location','NorthWest');
set(hleg2,'FontAngle','italic','TextColor',[.3,.2,.1],'FontSize',7)
%% Deviation of traction angle (DTA)
subplot(1,2,2)
bar(DTA), hold on
errorbar(DTA,DTAerr)
title('Deviation of Traction Angle')
xlabel({'FTTCf','0thf','2ndf','1norm','1norm2nd'});
%% save
hgexport(h2,strcat(imgPath,'S2 peak plot and DTA'),hgexport('factorystyle'),'Format','tiff')
hgsave(h2,strcat(imgPath,'S2 peak plot and DTA'),'-v7.3')
solTools = {sol_mats_regfine.tool, sol_mats_reg2fine.tool, sol_mats_1nfine.tool, sol_mats_1n2fine.tool};
regParams = [sol_mats_regfine.L sol_mats_reg2fine.L sol_mats_1nfine.L sol_mats_1n2fine.L];
save([dataPath '/all data.mat'],'DTMS','DTMSerr','DTA','DTAerr','pSR','pSRerr','flocMax','regParams','solTools','-v7.3');
