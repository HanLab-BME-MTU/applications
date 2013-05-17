%% compareFTTC_BEM_FastBEM
function [pSR,pSRerr,regParams] = calculateForceCG(imScale,forceType,percentNoise,savePath,forceMeshFastBEM,forceMeshFastBEMfine,M_FastBEM,M_FastBEMfine,basisClassTablePath,basisClassTablePathFine)
% calculateForceCG(forceType,percentNoise,savePath,forceMesh,M,M_FastBEM)
% This function uses conjugate gradient method to use sparse M.

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

meshPtsFwdSol=2^10;
% L=0;
csScale = imScale/5;
numPoints_u=50*imScale;       %(must be even number)
numPoints_f=50*csScale;       %(also must be even number - why???)
% numPoints_out=50;

xmin =1;
xmax =50*imScale;
ymin =1;
ymax =50*imScale;

xmaxf = 1+5*(xmax/5-1); % always 5 pixel grid for coarse grid
ymaxf = 1+5*(ymax/5-1);

imgPath = [savePath '/img/'];
dataPath = [savePath '/data/'];
if (~exist(imgPath,'dir') || ~exist(dataPath,'dir'))
    mkdir(imgPath);
    mkdir(dataPath);
end
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
        20*imScale,27*imScale; 
        21*imScale,20*imScale; % 2nd cluster
        
        27*imScale,39*imScale; 
        31*imScale,42*imScale; 
        37*imScale,43*imScale;
        31*imScale,36*imScale;
        38*imScale,40*imScale;
        35*imScale,35*imScale; % 3rd cluster

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
        850,-280;
        400,-10;

        500,-600;
        450,-700;
        250,-900;
        800,-1200;
        280,-800;
        700,-1000;

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
force_y7 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(7,1),fpos(7,2),fvec(7,1),fvec(7,2),forceType,'smallFA'); %group2
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
force_y14 = assumedForceShifted(2,x_mat_u,y_mat_u,fpos(14,1),fpos(14,2),fvec(14,1),fvec(14,2),forceType,'smallFA'); %group2
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
    assumedForceShifted(1,x,y,20*imScale,27*imScale,850,-280,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,21*imScale,20*imScale,400,-10,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,27*imScale,39*imScale,500,-600,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,31*imScale,42*imScale,450,-700,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,37*imScale,43*imScale,250,-900,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,31*imScale,36*imScale,800,-1200,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,38*imScale,40*imScale,280,-800,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,35*imScale,35*imScale,700,-1000,forceType,'smallFA')+...
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
    assumedForceShifted(2,x,y,20*imScale,27*imScale,850,-280,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,21*imScale,20*imScale,400,-10,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,27*imScale,39*imScale,500,-600,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,31*imScale,42*imScale,450,-700,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,37*imScale,43*imScale,250,-900,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,31*imScale,36*imScale,800,-1200,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,38*imScale,40*imScale,280,-800,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,35*imScale,35*imScale,700,-1000,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,33*imScale,9*imScale,-450,-800,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,39*imScale,18*imScale,-400,550,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,43*imScale,24*imScale,-350,400,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,38*imScale,12*imScale,-350,500,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,42*imScale,17*imScale,-150,600,forceType,'smallFA'),'fft',[],meshPtsFwdSol);
toc

display('forward solution for sparse solution...')
tic
% for sparse sampling
[x_mat_f, y_mat_f]=meshgrid(linspace(xmin,xmaxf,numPoints_f) , linspace(ymin,ymaxf,numPoints_f));

[ux_f, uy_f]=fwdSolution(x_mat_f,y_mat_f,E,[],[],[],[],...
    @(x,y) assumedForceShifted(1,x,y,10*imScale,10*imScale,1200,-100,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,12*imScale,11*imScale,1800,-160,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,15*imScale,10*imScale,1500,-130,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,11*imScale,14*imScale,1300,-110,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,15*imScale,26*imScale,1100,-380,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,17*imScale,31*imScale,1300,-480,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,20*imScale,27*imScale,850,-280,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,21*imScale,20*imScale,400,-10,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,27*imScale,39*imScale,500,-600,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,31*imScale,42*imScale,450,-700,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,37*imScale,43*imScale,250,-900,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,31*imScale,36*imScale,800,-1200,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,38*imScale,40*imScale,280,-800,forceType,'smallFA')+...
    assumedForceShifted(1,x,y,35*imScale,35*imScale,700,-1000,forceType,'smallFA')+...
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
    assumedForceShifted(2,x,y,20*imScale,27*imScale,850,-280,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,21*imScale,20*imScale,400,-10,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,27*imScale,39*imScale,500,-600,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,31*imScale,42*imScale,450,-700,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,37*imScale,43*imScale,250,-900,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,31*imScale,36*imScale,800,-1200,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,38*imScale,40*imScale,280,-800,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,35*imScale,35*imScale,700,-1000,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,33*imScale,9*imScale,-450,-800,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,39*imScale,18*imScale,-400,550,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,43*imScale,24*imScale,-350,400,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,38*imScale,12*imScale,-350,500,forceType,'smallFA')+...
    assumedForceShifted(2,x,y,42*imScale,17*imScale,-150,600,forceType,'smallFA'),'fft',[],meshPtsFwdSol);
% [ux uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,@(x,y) assumedForce(2,x,y),@(x,y) assumedForce(2,x,y),'fft',[],meshPtsFwdSol);
toc
display(['adding noise (' num2str(percentNoise*100) '%) ...'])
if addNoise==1
    max_f=max([max(max(ux_f)) max(max(uy_f))]); 
    noise_x=normrnd(0,percentNoise*max_f,numPoints_f,numPoints_f);
    noise_y=normrnd(0,percentNoise*max_f,numPoints_f,numPoints_f);
    ux_f=ux_f+noise_x;
    uy_f=uy_f+noise_y;
    
    max_u=max([max(max(ux)) max(max(uy))]); 
    noise_xu=normrnd(0,percentNoise*max_u,numPoints_u,numPoints_u);
    noise_yu=normrnd(0,percentNoise*max_u,numPoints_u,numPoints_u);
    ux=ux+noise_xu;
    uy=uy+noise_yu;
end

u_mat_f(:,:,1)=ux_f;
u_mat_f(:,:,2)=uy_f;
ux_f_vec=reshape(ux_f,[],1);
uy_f_vec=reshape(uy_f,[],1);

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
%% Mesh creation
x_vec_f=reshape(x_mat_f,[],1);
y_vec_f=reshape(y_mat_f,[],1);
if isempty(forceMeshFastBEM)
    display('creating mesh and basis function')
    tic
    forceMeshFastBEM=createMeshAndBasisFastBEM(x_vec_f,y_vec_f,true,[],0);
    toc
end
%% FTTC-reconstruction with regularization 
L=2e-7;
grid_mat_f(:,:,1)=x_mat_f;
grid_mat_f(:,:,2)=y_mat_f;
i_max = size(grid_mat_f,1);
j_max = size(grid_mat_f,2);
cluster_size = grid_mat_f(1,1,1) - grid_mat_f(2,2,1);
[~,~,force_FTTC, ~] = reg_fourier_TFM(grid_mat_f,u_mat_f,E,s, pix_durch_my, cluster_size, i_max, j_max, L);  

rec_force_FTTCreg(:,:,1)=reshape(force_FTTC(:,1),i_max,j_max);
rec_force_FTTCreg(:,:,2)=reshape(force_FTTC(:,2),i_max,j_max);
%% FTTC-reconstruction with regularization with fine mesh
L=2e-7;
grid_mat_u(:,:,1)=x_mat_u; %dense
grid_mat_u(:,:,2)=y_mat_u; %dense
i_max = size(grid_mat_u,1);
j_max = size(grid_mat_u,2);
cluster_size = grid_mat_u(1,1,1) - grid_mat_u(2,2,1);

[~,~,force_FTTCfine, ~] = reg_fourier_TFM(grid_mat_u,u_mat,E,s, pix_durch_my, cluster_size, i_max, j_max, L);  

rec_force_FTTCregfine(:,:,1)=reshape(force_FTTCfine(:,1),i_max,j_max);
rec_force_FTTCregfine(:,:,2)=reshape(force_FTTCfine(:,2),i_max,j_max);
%% FastBEM
% FastBEM with regularization
% basisClassTablePath = '/home/sh268/orchestra/home/Documents/TFM-simulation/basisFunction5050.mat';
L = 1e-7;
display(['FastBEM (' num2str(L) ')...'])
[fx_FastBEMreg,fy_FastBEMreg,~,~,M_FastBEM,~,~,~,sol_mats_reg]...
    = BEM_force_reconstruction(x_mat_f,y_mat_f,ux_f,uy_f,forceMeshFastBEM,...
    E,L,[],[],'fast',meshPtsFwdSol,...
    'QR','basisClassTblPath',basisClassTablePath,'imgRows',ymax,'imgCols',xmax,'fwdMap',M_FastBEM);
rec_force_FastBEMreg(:,:,1)=fx_FastBEMreg;
rec_force_FastBEMreg(:,:,2)=fy_FastBEMreg;
% save(strcat(dataPath,'FastBEMreg.mat'),'rec_force_FastBEMreg','sol_mats_reg'); % sol_mats_reg.L: regularization parameter
%% FastBEM with fine mesh (For figure 3) 
% Mesh creation
if isempty(forceMeshFastBEMfine)
    display('creating fine mesh and basis function (1x1)')
    tic
    forceMeshFastBEMfine=createMeshAndBasisFastBEM(x_vec_u,y_vec_u,true,[],0);
    toc
end
L = 1e-7;
display(['FastBEM with fine mesh (' num2str(L) ')...'])
[fx_FastBEMregfine,fy_FastBEMregfine,~,~,M_FastBEMfine,~,~,~,sol_mats_regfine]...
    = iterFastBEM(x_mat_u,y_mat_u,ux,uy,forceMeshFastBEMfine,...
    E,L,[],[],'fast',meshPtsFwdSol,...
    'CG','basisClassTblPath',basisClassTablePathFine,'imgRows',ymax,'imgCols',xmax,'fwdMap',M_FastBEMfine);
rec_force_FastBEMregfine(:,:,1)=fx_FastBEMregfine;
rec_force_FastBEMregfine(:,:,2)=fy_FastBEMregfine;
%% FastBEM-reconstruction with laplacian regularization
L = 1e-8; %sol_mats_reg.L;
display(['FastBEM with 2nd order regularization (' num2str(L) ')...'])
tic
[fx_FastBEMreg2,fy_FastBEMreg2,~,~,~,~,~,~,sol_mats_reg2]...
    = BEM_force_reconstruction(x_mat_f,y_mat_f,ux_f,uy_f,forceMeshFastBEM,...
    E,L,[],[],'fast',meshPtsFwdSol,...
    'LaplacianRegularization','basisClassTblPath',basisClassTablePath,...
    'imgRows',ymax,'imgCols',xmax,'fwdMap',M_FastBEM);
toc
rec_force_FastBEMreg2(:,:,1)=fx_FastBEMreg2;
rec_force_FastBEMreg2(:,:,2)=fy_FastBEMreg2;
%% FastBEM-reconstruction with laplacian regularization with fine mesh
L = 1e-8; %sol_mats_reg.L;
display(['FastBEM with 2nd order regularization (' num2str(L) ')...'])
tic
[fx_FastBEMreg2fine,fy_FastBEMreg2fine,~,~,~,~,~,~,sol_mats_reg2fine]...
    = iterFastBEM(x_mat_u,y_mat_u,ux,uy,forceMeshFastBEMfine,...
    E,L,[],[],'fast',meshPtsFwdSol,...
    'LaplacianRegularization','basisClassTblPath',basisClassTablePathFine,...
    'imgRows',ymax,'imgCols',xmax,'fwdMap',M_FastBEMfine);
toc
rec_force_FastBEMreg2fine(:,:,1)=fx_FastBEMreg2fine;
rec_force_FastBEMreg2fine(:,:,2)=fy_FastBEMreg2fine;


%% FastBEM-reconstruction with 1Norm-0th order regularization with fine mesh
L = 1e-4; %sol_mats_reg.L;
display(['FastBEM with 1Norm-0th order regularization with fine mesh (' num2str(L) ')...'])
tic
[fx_FastBEM1nfine,fy_FastBEM1nfine,~,~,~,~,~,~,sol_mats_1nfine]...
    = iterFastBEM(x_mat_u,y_mat_u,ux,uy,forceMeshFastBEMfine,...
    E,L,[],[],'fast',meshPtsFwdSol,...
    '1NormReg','basisClassTblPath',basisClassTablePathFine,...
    'imgRows',ymax,'imgCols',xmax,'fwdMap',M_FastBEMfine);
toc
rec_force_FastBEM1nfine(:,:,1)=fx_FastBEM1nfine;
rec_force_FastBEM1nfine(:,:,2)=fy_FastBEM1nfine;

% %% quick saving of M and forcemesh
% save([dataPath '/MandforceMesh.mat'],'M_FastBEM','M_FastBEMfine','forceMeshFastBEM','forceMeshFastBEMfine','-v7.3');
%% Quiver plots of BEM results compared to original, FTTC, regularized FTTC 
% h5 = figure;
% forceScale=.1*sqrt(max(max(force_x1.^2+force_y1.^2)));
% hold on
% quiver(x_mat_u,y_mat_u,force_x/forceScale,force_y/forceScale,0,'k');
% quiver(grid_mat_f(:,:,1),grid_mat_f(:,:,2),rec_force_FTTCreg(:,:,1)/forceScale,rec_force_FTTCreg(:,:,2)/forceScale,0,'g');
% quiver(x_out,y_out,rec_force_FastBEMreg(:,:,1)/forceScale,rec_force_FastBEMreg(:,:,2)/forceScale,0,'b');
% quiver(x_out,y_out,rec_force_FastBEMreg2(:,:,1)/forceScale,rec_force_FastBEMreg2(:,:,2)/forceScale,0,'r');
% 
% hold off
% xlim([xmin-1 xmax+1])
% ylim([ymin-1 ymax+1])
% % Make the text of the legend italic and color it brown
% hleg = legend('Org Force','Regularized FTTC',...
%     'FastBEMreg','FastBEM laplacian',...
%               'Location','NorthWest');
% set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
% title('reconstructed forces')
% % We see that non-regularized force has 1) incorrect directions on some
% % points near high forces 2) forces at the boundary of the image. In case
% % of regularized reconstructed forces, the orientation looks better, but
% % the it underestimates the high forces compared original given forces.
% % BEM is appears to be correct even without the regularization. It has no boundary effect.
% hgexport(h5,strcat(imgPath,'force vector field'),hgexport('factorystyle'),'Format','tiff')
% hgsave(h5,strcat(imgPath,'force vector field'),'-v7.3')
% delete(h5)

%% Heat map presentation for each force magnitude
pos = [x_vec_u y_vec_u]; %dense
pos_f = [x_vec_f y_vec_f]; %sparse
disp_vec = [ux_vec uy_vec]; %dense
disp_f_vec = [ux_f_vec uy_f_vec]; %sparse

rec_force_org_vec = [reshape(force_x,[],1) reshape(force_y,[],1)];
[grid_matfine,tmat_org, ~, ~] = interp_vec2grid(pos+disp_vec, rec_force_org_vec,[],grid_mat_u); %1:cluster size
tnorm_org = (tmat_org(:,:,1).^2 + tmat_org(:,:,2).^2).^0.5; %this should be fine mesh
fmat_org(:,:,1) = force_x;
fmat_org(:,:,2) = force_y;
fnorm_org = (force_x.^2 + force_y.^2).^0.5; %this should be fine mesh

rec_force_FTTCreg_vec = [reshape(rec_force_FTTCreg(:,:,1),[],1) reshape(rec_force_FTTCreg(:,:,2),[],1)];
[~,tmat_fttcreg, ~, ~] = interp_vec2grid(pos_f+disp_f_vec, rec_force_FTTCreg_vec,[],grid_mat_f); %1:cluster size
tnorm_fttcreg = (tmat_fttcreg(:,:,1).^2 + tmat_fttcreg(:,:,2).^2).^0.5;
fnorm_fttcreg = (rec_force_FTTCreg(:,:,1).^2 + rec_force_FTTCreg(:,:,2).^2).^0.5;

rec_force_FastBEMreg_vec = [reshape(rec_force_FastBEMreg(:,:,1),[],1) reshape(rec_force_FastBEMreg(:,:,2),[],1)];
[grid_mat,tmat_FastBEMreg, ~, ~] = interp_vec2grid(pos_f+disp_f_vec, rec_force_FastBEMreg_vec,[],grid_mat_f); %1:cluster size
tnorm_FastBEMreg = (tmat_FastBEMreg(:,:,1).^2 + tmat_FastBEMreg(:,:,2).^2).^0.5;
fmat_FastBEMreg(:,:,1) = reshape(rec_force_FastBEMreg(:,:,1),size(x_mat_f));
fmat_FastBEMreg(:,:,2) = reshape(rec_force_FastBEMreg(:,:,2),size(y_mat_f));
fnorm_FastBEMreg = (fmat_FastBEMreg(:,:,1).^2 + fmat_FastBEMreg(:,:,2).^2).^0.5;

rec_force_FastBEMreg2_vec = [reshape(rec_force_FastBEMreg2(:,:,1),[],1) reshape(rec_force_FastBEMreg2(:,:,2),[],1)];
[~,tmat_FastBEMreg2, ~, ~] = interp_vec2grid(pos_f+disp_f_vec, rec_force_FastBEMreg2_vec,[],grid_mat_f); %1:cluster size
tnorm_FastBEMreg2 = (tmat_FastBEMreg2(:,:,1).^2 + tmat_FastBEMreg2(:,:,2).^2).^0.5;
fmat_FastBEMreg2(:,:,1) = reshape(rec_force_FastBEMreg2(:,:,1),size(x_mat_f));
fmat_FastBEMreg2(:,:,2) = reshape(rec_force_FastBEMreg2(:,:,2),size(y_mat_f));
fnorm_FastBEMreg2 = (fmat_FastBEMreg2(:,:,1).^2 + fmat_FastBEMreg2(:,:,2).^2).^0.5;

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
fnorm_FastBEMn1fine = (fmat_FastBEM1nfine(:,:,1).^2 + fmat_FastBEM1nfine(:,:,2).^2).^0.5;

% %% Residuals
% % magnitude residual
% res_tnorm_fttc = tnorm_org-tnorm_fttc;
% res_tnorm_fttcreg = tnorm_org-tnorm_fttcreg;
% res_tnorm_BEM = tnorm_org-tnorm_BEM;
% res_tnorm_BEMreg = tnorm_org-tnorm_BEMreg;
% res_tnorm_FastBEM = tnorm_org-tnorm_FastBEM;
% res_tnorm_FastBEMreg = tnorm_org-tnorm_FastBEMreg;
% res_tnorm_FastBEMregsc = tnorm_org-tnorm_FastBEMregsc;
% res_tnorm_FastBEMreg2 = tnorm_org-tnorm_FastBEMreg2;
% 
% h52=figure(52); 
% set(h52, 'Position', [100 300 900 900])
% tmin = min(min(res_tnorm_fttcreg));
% tmax = max(max(res_tnorm_fttcreg));
% tmin = min(tmin,min(min(res_tnorm_fttc)));
% tmax = max(tmax,max(max(res_tnorm_fttc)));
% colormap('jet');
% subplot(3,3,1), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_fttc), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FTTC residual')
% caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
% xlim([xmin xmax])
% ylim([ymin ymax])
% 
% subplot(3,3,2), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_fttcreg), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FTTC regularization res')
% caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
% xlim([xmin xmax])
% ylim([ymin ymax])
% 
% subplot(3,3,4), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_BEM), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('BEM residual')
% caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
% xlim([xmin xmax])
% ylim([ymin ymax])
% 
% subplot(3,3,5), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_BEMreg), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('BEM regularization res')
% caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
% xlim([xmin xmax])
% ylim([ymin ymax])
% 
% subplot(3,3,7), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_FastBEM), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FastBEM residual')
% caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
% xlim([xmin xmax])
% ylim([ymin ymax])
% 
% subplot(3,3,8), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_FastBEMreg), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FastBEM regularization res')
% caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
% xlim([xmin xmax])
% ylim([ymin ymax])
% 
% subplot(3,3,9), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_FastBEMregsc), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FastBEM regularization res')
% caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
% xlim([xmin xmax])
% ylim([ymin ymax])
% 
% subplot(3,3,10), surf(grid_mat(:,:,1), grid_mat(:,:,2), res_tnorm_FastBEMreg2), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FastBEM 2nd order regularization')
% caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
% xlim([xmin xmax])
% ylim([ymin ymax])
% 
% subplot(3,3,6), caxis([tmin-0.01 tmax+0.01]), axis off; colorbar('West')
%% Total residual
% tot_res(1) = sum((rec_force_org_vec(:,1)-rec_force_FTTCreg_vec(:,1)).^2+...
%     (rec_force_org_vec(:,2)-rec_force_FTTCreg_vec(:,2)).^2)^.5;
% tot_res(2) = sum((rec_force_org_vec(:,1)-rec_force_FastBEMreg_vec(:,1)).^2+...
%     (rec_force_org_vec(:,2)-rec_force_FastBEMreg_vec(:,2)).^2)^.5;
% tot_res(3) = sum((rec_force_org_vec(:,1)-rec_force_FastBEMreg2_vec(:,1)).^2+...
%     (rec_force_org_vec(:,2)-rec_force_FastBEMreg2_vec(:,2)).^2)^.5;
% first one is fttc without regularization, second bar is total residual
% in regularized fttc. Total residual is lower in regularized fttc.

%% Force at cross-sectional line at y=10 for large FA, high peak
% For force 1
tnorm1D_org_g1 = fnorm_org(11*imScale,12*imScale-5*imScale:12*imScale+5*imScale);
tnorm1D_FTTCreg_g1 = fnorm_fttcreg(round(11*csScale)+1,round(12*csScale)+1-round(5*csScale):round(12*csScale)+1+round(5*csScale));
tnorm1D_FastBEMreg_g1 = fnorm_FastBEMreg(round(11*csScale)+1,round(12*csScale)+1-round(5*csScale):round(12*csScale)+1+round(5*csScale));
tnorm1D_FastBEMreg2_g1 = fnorm_FastBEMreg2(round(11*csScale)+1,round(12*csScale)+1-round(5*csScale):round(12*csScale)+1+round(5*csScale));
tnorm1D_FTTCregfine_g1 = fnorm_fttcregfine(11*imScale,12*imScale-5*imScale:12*imScale+5*imScale);
tnorm1D_FastBEMregfine_g1 = fnorm_FastBEMregfine(11*imScale,12*imScale-5*imScale:12*imScale+5*imScale);
tnorm1D_FastBEMreg2fine_g1 = fnorm_FastBEMreg2fine(11*imScale,12*imScale-5*imScale:12*imScale+5*imScale);
tnorm1D_FastBEMn1fine_g1 = fnorm_FastBEMn1fine(11*imScale,12*imScale-5*imScale:12*imScale+5*imScale);
x1D_g1 = grid_mat_f(round(11*csScale)+1,round(12*csScale)+1-round(5*csScale):round(12*csScale)+1+round(5*csScale),1);
x1Du_g1 = grid_mat_u(11*imScale,12*imScale-5*imScale:12*imScale+5*imScale,1);
%% Force at cross-sectional line at force_4 for large FA, low peak
% force_x4 = assumedForceShifted(1,x_mat_u,y_mat_u,37*10,45*10,300,-800,forceType,'largeFA');
% force_y4 = assumedForceShifted(2,x_mat_u,y_mat_u,37*10,45*10,300,-800,forceType,'largeFA'); %group2
tnorm1D_org_g2 = fnorm_org(36*imScale,31*imScale-5*imScale:31*imScale+5*imScale);
tnorm1D_FTTCreg_g2 = fnorm_fttcreg(round(36*csScale)+1,round(31*csScale)+1-round(5*csScale):round(31*csScale)+1+round(5*csScale));
tnorm1D_FastBEMreg_g2 = fnorm_FastBEMreg(round(36*csScale)+1,round(31*csScale)+1-round(5*csScale):round(31*csScale)+1+round(5*csScale));
tnorm1D_FastBEMreg2_g2 = fnorm_FastBEMreg2(round(36*csScale)+1,round(31*csScale)+1-round(5*csScale):round(31*csScale)+1+round(5*csScale));
tnorm1D_FTTCregfine_g2 = fnorm_fttcregfine(36*imScale,31*imScale-5*imScale:31*imScale+5*imScale);
tnorm1D_FastBEMregfine_g2 = fnorm_FastBEMregfine(36*imScale,31*imScale-5*imScale:31*imScale+5*imScale);
tnorm1D_FastBEMreg2fine_g2 = fnorm_FastBEMreg2fine(36*imScale,31*imScale-5*imScale:31*imScale+5*imScale);
tnorm1D_FastBEMn1fine_g2 = fnorm_FastBEMn1fine(36*imScale,31*imScale-5*imScale:31*imScale+5*imScale);
x1D_g2 = grid_mat_f(round(36*csScale)+1,round(31*csScale)+1-round(5*csScale):round(31*csScale)+1+round(5*csScale),1);
x1Du_g2 = grid_mat_u(36*imScale,31*imScale-5*imScale:31*imScale+5*imScale,1);
%% Force at cross-sectional line at force_7 for small FA, high peak
tnorm1D_org_g3 = fnorm_org(31*imScale,17*imScale-5*imScale:17*imScale+5*imScale);
tnorm1D_FTTCreg_g3 = fnorm_fttcreg(round(31*csScale)+1,round(17*csScale)+1-round(5*csScale):round(17*csScale)+1+round(5*csScale));
tnorm1D_FastBEMreg_g3 = fnorm_FastBEMreg(round(31*csScale)+1,round(17*csScale)+1-round(5*csScale):round(17*csScale)+1+round(5*csScale));
tnorm1D_FastBEMreg2_g3 = fnorm_FastBEMreg2(round(31*csScale)+1,round(17*csScale)+1-round(5*csScale):round(17*csScale)+1+round(5*csScale));
tnorm1D_FTTCregfine_g3 = fnorm_fttcregfine(31*imScale,17*imScale-5*imScale:17*imScale+5*imScale);
tnorm1D_FastBEMregfine_g3 = fnorm_FastBEMregfine(31*imScale,17*imScale-5*imScale:17*imScale+5*imScale);
tnorm1D_FastBEMreg2fine_g3 = fnorm_FastBEMreg2fine(31*imScale,17*imScale-5*imScale:17*imScale+5*imScale);
tnorm1D_FastBEMn1fine_g3 = fnorm_FastBEMn1fine(31*imScale,17*imScale-5*imScale:17*imScale+5*imScale);
x1D_g3 = grid_mat_f(round(31*csScale)+1,round(17*csScale)+1-round(5*csScale):round(17*csScale)+1+round(5*csScale),1);
x1Du_g3 = grid_mat_u(31*imScale,17*imScale-5*imScale:17*imScale+5*imScale,1);
%% Force at cross-sectional line at force_10 for small FA, low peak
% force_x10 = assumedForceShifted(1,x_mat_u,y_mat_u,12*10,20*10,600,-480,forceType,'smallFA');
% force_y10 = assumedForceShifted(2,x_mat_u,y_mat_u,12*10,20*10,600,-480,forceType,'smallFA'); %group4
tnorm1D_org_g4 = fnorm_org(20*imScale,21*imScale-5*imScale:21*imScale+5*imScale);
tnorm1D_FTTCreg_g4 = fnorm_fttcreg(round(20*csScale)+1,round(21*csScale)+1-round(5*csScale):round(21*csScale)+1+round(5*csScale));
tnorm1D_FastBEMreg_g4 = fnorm_FastBEMreg(round(20*csScale)+1,round(21*csScale)+1-round(5*csScale):round(21*csScale)+1+round(5*csScale));
tnorm1D_FastBEMreg2_g4 = fnorm_FastBEMreg2(round(20*csScale)+1,round(21*csScale)+1-round(5*csScale):round(21*csScale)+1+round(5*csScale));
tnorm1D_FTTCregfine_g4 = fnorm_fttcregfine(20*imScale,21*imScale-5*imScale:21*imScale+5*imScale);
tnorm1D_FastBEMregfine_g4 = fnorm_FastBEMregfine(20*imScale,21*imScale-5*imScale:21*imScale+5*imScale);
tnorm1D_FastBEMreg2fine_g4 = fnorm_FastBEMreg2fine(20*imScale,21*imScale-5*imScale:21*imScale+5*imScale);
tnorm1D_FastBEMn1fine_g4 = fnorm_FastBEMn1fine(20*imScale,21*imScale-5*imScale:21*imScale+5*imScale);
x1D_g4 = grid_mat_f(round(20*csScale)+1,round(21*csScale)+1-round(5*csScale):round(21*csScale)+1+round(5*csScale),1);
x1Du_g4 = grid_mat_u(20*imScale,21*imScale-5*imScale:21*imScale+5*imScale,1);
% %% Boundary cutting for FastBEM
% tnorm_FastBEMregsc(1,:)=0;
% tnorm_FastBEMregsc(:,1)=0;
% tnorm_FastBEMregsc(end,:)=0;
% tnorm_FastBEMregsc(:,end)=0;
% tnorm_FastBEMreg2(1,:)=0;
% tnorm_FastBEMreg2(:,1)=0;
% tnorm_FastBEMreg2(end,:)=0;
% tnorm_FastBEMreg2(:,end)=0;

%% Precision in locating peak force
tnorm(:,:,1) = fnorm_fttcreg;
tnorm(:,:,2) = fnorm_FastBEMreg;
tnorm(:,:,3) = fnorm_FastBEMreg2;
tnormd(:,:,1) = fnorm_fttcregfine;
tnormd(:,:,2) = fnorm_FastBEMregfine;
tnormd(:,:,3) = fnorm_FastBEMreg2fine;
tnormd(:,:,4) = fnorm_FastBEMn1fine;
nMethods = 7;

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
backgroundMask = sqrt((x_mat_f-(locMaxI(1,1)/5+1)).^2+(y_mat_f-(locMaxI(1,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(2,1)/5+1)).^2+(y_mat_f-(locMaxI(2,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(3,1)/5+1)).^2+(y_mat_f-(locMaxI(3,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(4,1)/5+1)).^2+(y_mat_f-(locMaxI(4,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(5,1)/5+1)).^2+(y_mat_f-(locMaxI(5,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(6,1)/5+1)).^2+(y_mat_f-(locMaxI(6,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(7,1)/5+1)).^2+(y_mat_f-(locMaxI(7,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(8,1)/5+1)).^2+(y_mat_f-(locMaxI(8,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(9,1)/5+1)).^2+(y_mat_f-(locMaxI(9,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(10,1)/5+1)).^2+(y_mat_f-(locMaxI(10,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(11,1)/5+1)).^2+(y_mat_f-(locMaxI(11,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(12,1)/5+1)).^2+(y_mat_f-(locMaxI(12,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(13,1)/5+1)).^2+(y_mat_f-(locMaxI(13,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(14,1)/5+1)).^2+(y_mat_f-(locMaxI(14,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(15,1)/5+1)).^2+(y_mat_f-(locMaxI(15,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(16,1)/5+1)).^2+(y_mat_f-(locMaxI(16,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(17,1)/5+1)).^2+(y_mat_f-(locMaxI(17,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(18,1)/5+1)).^2+(y_mat_f-(locMaxI(18,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_f-(locMaxI(19,1)/5+1)).^2+(y_mat_f-(locMaxI(19,2)/5+1)).^2)>6;
backgroundMaskfine = sqrt((x_mat_u-(locMaxI(1,1))).^2+(y_mat_u-(locMaxI(1,2))).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(2,1))).^2+(y_mat_u-(locMaxI(2,2))).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(3,1))).^2+(y_mat_u-(locMaxI(3,2))).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(4,1))).^2+(y_mat_u-(locMaxI(4,2))).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(5,1))).^2+(y_mat_u-(locMaxI(5,2))).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(6,1)/5+1)).^2+(y_mat_u-(locMaxI(6,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(7,1)/5+1)).^2+(y_mat_u-(locMaxI(7,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(8,1)/5+1)).^2+(y_mat_u-(locMaxI(8,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(9,1)/5+1)).^2+(y_mat_u-(locMaxI(9,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(10,1)/5+1)).^2+(y_mat_u-(locMaxI(10,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(11,1)/5+1)).^2+(y_mat_u-(locMaxI(11,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(12,1)/5+1)).^2+(y_mat_u-(locMaxI(12,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(13,1)/5+1)).^2+(y_mat_u-(locMaxI(13,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(14,1)/5+1)).^2+(y_mat_u-(locMaxI(14,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(15,1)/5+1)).^2+(y_mat_u-(locMaxI(15,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(16,1)/5+1)).^2+(y_mat_u-(locMaxI(16,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(17,1)/5+1)).^2+(y_mat_u-(locMaxI(17,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(18,1)/5+1)).^2+(y_mat_u-(locMaxI(18,2)/5+1)).^2)>6 ...
                    & sqrt((x_mat_u-(locMaxI(19,1)/5+1)).^2+(y_mat_u-(locMaxI(19,2)/5+1)).^2)>6;

for jj=1:3 
    for k=1:nPeaks
        flocMax(k,jj) = tnorm(round(locMaxI(k,1)/5)+1,round(locMaxI(k,2)/5)+1,jj);
        flocMax(k,jj+3) = tnormd(locMaxI(k,1),locMaxI(k,2),jj);
        flocMaxRatio(k,jj) = flocMax(k,jj)/fnorm_org(locMaxI(k,1),locMaxI(k,2)); % ratio based on original force norm (0~1)
        flocMaxRatio(k,jj+3) = flocMax(k,jj+3)/fnorm_org(locMaxI(k,1),locMaxI(k,2)); % ratio based on original force norm (0~1)
        if jj==1
            recF_vec = [rec_force_FTTCreg(round(locMaxI(k,1)/5)+1,round(locMaxI(k,2)/5)+1,1),...
                rec_force_FTTCreg(round(locMaxI(k,1)/5)+1,round(locMaxI(k,2)/5)+1,2)];
            recF_vecfine = [rec_force_FTTCregfine(locMaxI(k,1),locMaxI(k,2),1),...
                rec_force_FTTCregfine(locMaxI(k,1),locMaxI(k,2),2)];
        elseif jj==2
            recF_vec = [fmat_FastBEMreg(round(locMaxI(k,1)/5)+1,round(locMaxI(k,2)/5)+1,1),...
                fmat_FastBEMreg(round(locMaxI(k,1)/5)+1,round(locMaxI(k,2)/5)+1,2)];
            recF_vecfine = [fmat_FastBEMregfine(locMaxI(k,1)+1,locMaxI(k,2),1),...
                fmat_FastBEMregfine(locMaxI(k,1),locMaxI(k,2),2)];
        elseif jj==3
            recF_vec = [fmat_FastBEMreg2(round(locMaxI(k,1)/5)+1,round(locMaxI(k,2)/5)+1,1),...
                fmat_FastBEMreg2(round(locMaxI(k,1)/5)+1,round(locMaxI(k,2)/5)+1,2)];
            recF_vecfine = [fmat_FastBEMreg2fine(locMaxI(k,1),locMaxI(k,2),1),...
                fmat_FastBEMreg2fine(locMaxI(k,1),locMaxI(k,2),2)];
        end
        orgF_vec = [fmat_org(locMaxI(k,1),locMaxI(k,2),1),...
                    fmat_org(locMaxI(k,1),locMaxI(k,2),2)];

        flocMaxAngle(k,jj) = acosd(recF_vec*orgF_vec'/...
            (norm(recF_vec)*norm(orgF_vec))); % in degree
        flocMaxAngle(k,jj+3) = acosd(recF_vecfine*orgF_vec'/...
            (norm(recF_vecfine)*norm(orgF_vec))); % in degree
        %indices for surroundings (2 um radius)
        if k==1 || k==2 || k==5
            r = 3;
        else
            r=3;
        end
        ringMask = sqrt((x_mat_f-(locMaxI(k,1)/5+1)).^2+(y_mat_f-(locMaxI(k,2)/5+1)).^2)>r ...
                    & sqrt((x_mat_f-(locMaxI(k,1)/5+1)).^2+(y_mat_f-(locMaxI(k,2)/5+1)).^2)<=r+5 ...
                    & backgroundMask;
        ringMaskfine = sqrt((x_mat_u-locMaxI(k,1)).^2+(y_mat_u-locMaxI(k,2)).^2)>r ...
                    & sqrt((x_mat_u-(locMaxI(k,1))).^2+(y_mat_u-(locMaxI(k,2))).^2)<=r+1 ...
                    & backgroundMaskfine;
        neighborTnorm = ringMask.*tnorm(:,:,jj);
        neighborTnormfine = ringMaskfine.*tnormd(:,:,jj);
        meanNeiTnorm = sum(sum(neighborTnorm))/sum(sum(ringMask));
        meanNeiTnormfine = sum(sum(neighborTnormfine))/sum(sum(ringMaskfine));
        flocMaxDTMS(k,jj) = meanNeiTnorm/flocMax(k,jj);
        flocMaxDTMS(k,jj+3) = meanNeiTnormfine/flocMax(k,jj+3);
        % Peak force localization(position) match - skipped for now.
    end
    % statistic of flocMaxRatio
    pSR(jj) = mean(flocMaxRatio(:,jj));
    pSRerr(jj) = sqrt(std(flocMaxRatio(:,jj))/nPeaks);
    DTA(jj) = mean(flocMaxAngle(:,jj));
    DTAerr(jj) = sqrt(std(flocMaxAngle(:,jj))/nPeaks);
    DTMS(jj) = mean(flocMaxDTMS(:,jj));
    DTMSerr(jj) = sqrt(std(flocMaxDTMS(:,jj))/nPeaks);
    pSR(jj+3) = mean(flocMaxRatio(:,jj+3));
    pSRerr(jj+3) = sqrt(std(flocMaxRatio(:,jj+3))/nPeaks);
    DTA(jj+3) = mean(flocMaxAngle(:,jj+3));
    DTAerr(jj+3) = sqrt(std(flocMaxAngle(:,jj+3))/nPeaks);
    DTMS(jj+3) = mean(flocMaxDTMS(:,jj+3));
    DTMSerr(jj+3) = sqrt(std(flocMaxDTMS(:,jj+3))/nPeaks);
end

% for jj=7
jj=4;
for k=1:nPeaks
    flocMax(k,jj+3) = tnormd(locMaxI(k,1),locMaxI(k,2),jj);
    flocMaxRatio(k,jj+3) = flocMax(k,jj+3)/fnorm_org(locMaxI(k,1),locMaxI(k,2)); % ratio based on original force norm (0~1)
    recF_vecfine = [fmat_FastBEMreg2fine(locMaxI(k,1),locMaxI(k,2),1),...
        fmat_FastBEMreg2fine(locMaxI(k,1),round(locMaxI(k,2))+1,2)];
    orgF_vec = [fmat_org(locMaxI(k,1),locMaxI(k,2),1),...
                fmat_org(locMaxI(k,1),locMaxI(k,2),2)];
    flocMaxAngle(k,jj+3) = acosd(recF_vecfine*orgF_vec'/...
        (norm(recF_vecfine)*norm(orgF_vec))); % in degree
    %indices for surroundings (2 um radius)
    r=3;
    ringMaskfine = sqrt((x_mat_u-locMaxI(k,1)).^2+(y_mat_u-locMaxI(k,2)).^2)>r ...
                & sqrt((x_mat_u-(locMaxI(k,1))).^2+(y_mat_u-(locMaxI(k,2))).^2)<=r+1 ...
                & backgroundMaskfine;
    neighborTnormfine = ringMaskfine.*tnormd(:,:,jj);
    meanNeiTnormfine = sum(sum(neighborTnormfine))/sum(sum(ringMaskfine));
    flocMaxDTMS(k,jj+3) = meanNeiTnormfine/flocMax(k,jj+3);
end
% statistic of flocMaxRatio
pSR(jj+3) = mean(flocMaxRatio(:,jj+3));
pSRerr(jj+3) = sqrt(std(flocMaxRatio(:,jj+3))/nPeaks);
DTA(jj+3) = mean(flocMaxAngle(:,jj+3));
DTAerr(jj+3) = sqrt(std(flocMaxAngle(:,jj+3))/nPeaks);
DTMS(jj+3) = mean(flocMaxDTMS(:,jj+3));
DTMSerr(jj+3) = sqrt(std(flocMaxDTMS(:,jj+3))/nPeaks);

%% Figures
%% Quiver plots of BEM results compared to original, FTTC, regularized FTTC 
h5 = figure;
forceScale=.5*sqrt(max(max(force_x1.^2+force_y1.^2)));
hold on
quiver(x_mat_u,y_mat_u,force_x/forceScale,force_y/forceScale,0,'k');
quiver(grid_mat_f(:,:,1),grid_mat_f(:,:,2),rec_force_FTTCreg(:,:,1)/forceScale,rec_force_FTTCreg(:,:,2)/forceScale,0,'c');
quiver(grid_mat_f(:,:,1),grid_mat_f(:,:,2),fmat_FastBEMreg(:,:,1)/forceScale,fmat_FastBEMreg(:,:,2)/forceScale,0,'y');
quiver(grid_mat_f(:,:,1),grid_mat_f(:,:,2),fmat_FastBEMreg2(:,:,1)/forceScale,fmat_FastBEMreg2(:,:,2)/forceScale,0,'m');
quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),rec_force_FTTCregfine(:,:,1)/forceScale,rec_force_FTTCregfine(:,:,2)/forceScale,0,'g');
quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),fmat_FastBEMregfine(:,:,1)/forceScale,fmat_FastBEMregfine(:,:,2)/forceScale,0,'b');
quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),fmat_FastBEMreg2fine(:,:,1)/forceScale,fmat_FastBEMreg2fine(:,:,2)/forceScale,0,'r');
quiver(grid_mat_u(:,:,1),grid_mat_u(:,:,2),fmat_FastBEM1nfine(:,:,1)/forceScale,fmat_FastBEM1nfine(:,:,2)/forceScale,0,'r');

hold off
xlim([xmin-1 xmax+1])
ylim([ymin-1 ymax+1])
% Make the text of the legend italic and color it brown
hleg = legend('Org Force','FTTC',...
    'FastBEMreg','FastBEM laplacian',...
    'FTTC fine','FastBEM fine','Lap fine','L1Norm',...
              'Location','NorthEast');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1],'Fontsize',6)
title('reconstructed forces')
hgexport(h5,strcat(imgPath,'force vector field'),hgexport('factorystyle'),'Format','tiff')
hgsave(h5,strcat(imgPath,'force vector field'),'-v7.3')
% delete(h5)

%% Heatmaps
h1 = figure;
hold off
set(h1, 'Position', [100 100 1500 1200])
%% Heatmaps
tmin = min(min(tnorm_org))-500;
tmax = max(max(tnorm_org))*1.3;
colormap('jet');
subplot(4,5,1), surf(grid_matfine(:,:,1), grid_matfine(:,:,2), tnorm_org,'FaceColor','interp',...
	'EdgeColor','none', 'FaceLighting','phong'), zlim([tmin-1 tmax+1]), ...
    view(0,90), title('Original force')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
% % hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
xlim([xmin xmax])
ylim([ymin ymax])

subplot(4,5,2), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_fttcreg), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FTTC regularization')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
xlim([xmin xmax])
ylim([ymin ymax])

subplot(4,5,3), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_FastBEMreg), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FastBEM regularization')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
xlim([xmin xmax])
ylim([ymin ymax])

subplot(4,5,4), surf(grid_mat(:,:,1), grid_mat(:,:,2), tnorm_FastBEMreg2), zlim([tmin-1 tmax+1])
view(0,90), shading interp; title('FastBEM 2nd order regularization')
set(gca, 'DataAspectRatio', [1,1,1]);
hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
xlim([xmin xmax])
ylim([ymin ymax])
caxis([tmin-0.01 tmax+0.01]);


subplot(4,5,6), surf(grid_matfine(:,:,1), grid_matfine(:,:,2), tnorm_fttcregfine), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FTTC regularization')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
% hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
xlim([xmin xmax])
ylim([ymin ymax])

subplot(4,5,7), surf(grid_matfine(:,:,1), grid_matfine(:,:,2), tnorm_FastBEMregfine), zlim([tmin-1 tmax+1]), view(0,90), shading interp; title('FastBEM regularization')
caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
% hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
xlim([xmin xmax])
ylim([ymin ymax])

subplot(4,5,8), surf(grid_matfine(:,:,1), grid_matfine(:,:,2), tnorm_FastBEMreg2fine), zlim([tmin-1 tmax+1])
view(0,90), shading interp; title('FastBEM 2nd order regularization')
set(gca, 'DataAspectRatio', [1,1,1]);
% hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
xlim([xmin xmax])
ylim([ymin ymax])
caxis([tmin-0.01 tmax+0.01]);

subplot(4,5,9), surf(grid_matfine(:,:,1), grid_matfine(:,:,2), tnorm_FastBEM1nfine), zlim([tmin-1 tmax+1])
view(0,90), shading interp; title('FastBEM 1norm regularization')
set(gca, 'DataAspectRatio', [1,1,1]);
% hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
xlim([xmin xmax])
ylim([ymin ymax])
caxis([tmin-0.01 tmax+0.01]);

subplot(4,5,10), caxis([tmin-0.01 tmax+0.01]), axis off; colorbar('West')
%% profile plot
subplot(4,5,11), %for high peak
plot(x1Du_g1,tnorm1D_org_g1,'k','Linewidth',2)
hold on
plot(x1D_g1,tnorm1D_FTTCreg_g1,'g')
plot(x1D_g1,tnorm1D_FastBEMreg_g1,'b')
plot(x1D_g1,tnorm1D_FastBEMreg2_g1,'r')
plot(x1Du_g1,tnorm1D_FTTCregfine_g1,'g','Linewidth',2)
plot(x1Du_g1,tnorm1D_FastBEMregfine_g1,'b','Linewidth',2)
plot(x1Du_g1,tnorm1D_FastBEMreg2fine_g1,'r','Linewidth',2)
plot(x1Du_g1,tnorm1D_FastBEMn1fine_g1,'m','Linewidth',2)
xlabel('x')
ylabel('Traction Stress (Pa)')
ymax1 = max(tnorm1D_org_g1)*1.5;
ylim([0 ymax1]);
title({'1D stress profile for','large FA, high peak'})

subplot(4,5,12), %for high peak
plot(x1Du_g2,tnorm1D_org_g2,'k','Linewidth',2)
hold on
plot(x1D_g2,tnorm1D_FTTCreg_g2,'g')
plot(x1D_g2,tnorm1D_FastBEMreg_g2,'b')
plot(x1D_g2,tnorm1D_FastBEMreg2_g2,'r')
plot(x1Du_g2,tnorm1D_FTTCregfine_g2,'g','Linewidth',2)
plot(x1Du_g2,tnorm1D_FastBEMregfine_g2,'b','Linewidth',2)
plot(x1Du_g2,tnorm1D_FastBEMreg2fine_g2,'r','Linewidth',2)
plot(x1Du_g2,tnorm1D_FastBEMn1fine_g2,'m','Linewidth',2)
xlabel('x')
ylabel('Traction Stress (Pa)')
ylim([0 ymax1]);
title({'1D stress profile for','large FA, low peak'})

subplot(4,5,13), %for high peak
plot(x1Du_g3,tnorm1D_org_g3,'k','Linewidth',2)
hold on
plot(x1D_g3,tnorm1D_FTTCreg_g3,'g')
plot(x1D_g3,tnorm1D_FastBEMreg_g3,'b')
plot(x1D_g3,tnorm1D_FastBEMreg2_g3,'r')
plot(x1Du_g3,tnorm1D_FTTCregfine_g3,'g','Linewidth',2)
plot(x1Du_g3,tnorm1D_FastBEMregfine_g3,'b','Linewidth',2)
plot(x1Du_g3,tnorm1D_FastBEMreg2fine_g3,'r','Linewidth',2)
plot(x1Du_g3,tnorm1D_FastBEMn1fine_g3,'m','Linewidth',2)
xlabel('x')
ylabel('Traction Stress (Pa)')
ylim([0 ymax1]);
title({'1D stress profile for','small FA, high peak'})

subplot(4,5,14), %for low peak
plot(x1Du_g4,tnorm1D_org_g4,'k','Linewidth',2)
hold on
plot(x1D_g4,tnorm1D_FTTCreg_g4,'g')
plot(x1D_g4,tnorm1D_FastBEMreg_g4,'b')
plot(x1D_g4,tnorm1D_FastBEMreg2_g4,'r')
plot(x1Du_g4,tnorm1D_FTTCregfine_g4,'g','Linewidth',2)
plot(x1Du_g4,tnorm1D_FastBEMregfine_g4,'b','Linewidth',2)
plot(x1Du_g4,tnorm1D_FastBEMreg2fine_g4,'r','Linewidth',2)
plot(x1Du_g4,tnorm1D_FastBEMn1fine_g4,'r','Linewidth',2)
xlabel('x')
ylabel('Traction Stress (Pa)')
ylim([0 ymax1]);
title({'1D stress profile for','small FA, low peak'})
hleg = legend('Original','FTTC','FastBEM0th','FastBEM2nd','FTTCfine','BEM0thfine','2ndfine','Location','North');
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1],'FontSize',7)

% subplot(4,4,8),
% h1D_org = plot(x1Du_g4,tnorm1D_org_g4,'k','Linewidth',2);
% hold on
% h1D_fttc = plot(x1D_g4,tnorm1D_FTTCreg_g4,'g');
% h1D_reg = plot(x1D_g4,tnorm1D_FastBEMreg_g4,'b');
% h1D_reg2 = plot(x1D_g4,tnorm1D_FastBEMreg2_g4,'r');
% axis off;
% hleg = legend('Original Stress','FTTC reg','FastBEM reg','FastBEM laplacian','Location','South');
% set([h1D_org h1D_reg h1D_reg2 h1D_fttc],'Visible','off')
% set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1],'FontSize',7)

% %% Mean Squared Deviation
% subplot(4,4,12), bar(tot_res)
% ylim([4000 max(tot_res(1))*1.1]);
% 
% title('mean squared deviation')
% xlabel({'FTTC reg vs.','vs. FastBEM reg ',' vs. FastBEM 2nd reg'});
%% Peak plot
subplot(4,5,16)
plot(flocMaxOrg/1000,flocMax(:,1)/1000,'g.'), hold on
plot(flocMaxOrg/1000,flocMax(:,2)/1000,'b.'), 
plot(flocMaxOrg/1000,flocMax(:,3)/1000,'r.'), 
plot(flocMaxOrg/1000,flocMax(:,4)/1000,'go'), 
plot(flocMaxOrg/1000,flocMax(:,5)/1000,'bo'), 
plot(flocMaxOrg/1000,flocMax(:,6)/1000,'ro'), 
plot(flocMaxOrg/1000,flocMax(:,7)/1000,'mo'), 
plot(0:.1:3.2,0:.1:3.2,'k--'),hold off
xlim([0 2.7]),ylim([0 3.3]);
xlabel('Original Peak Stress (kPa)')
ylabel('Reconstructed Peak Stress (kPa)')
title({'Reconstructed Peak Stress','w.r.t. Original Peak Stress'})
hleg2=legend('FTTC','0th','2nd','FTTCf','0thf','2ndf','1norm','1:1','Location','NorthWest');
set(hleg2,'FontAngle','italic','TextColor',[.3,.2,.1],'FontSize',7)

%% Peak Magnitude Ratio
subplot(4,5,17),
bar(pSR), hold on
errorbar(pSR,pSRerr)
title('Peak Stress ratio')
xlabel({'FTTC reg vs. ',' FastBEM reg vs. ',' FastBEM 2nd reg',' FastBEM 1norm'});
%% Deviation of traction angle (DTA)
subplot(4,5,18),
bar(DTA), hold on
errorbar(DTA,DTAerr)
title('Deviation of Traction Angle')
xlabel({'FTTC reg vs. ',' FastBEM reg vs. ',' FastBEM 2nd reg',' FastBEM 1norm'});
%% Deviation of traction magnitude surroundings (DTMS)
subplot(4,5,19),
bar(DTMS), hold on
errorbar(DTMS,DTMSerr)
title('Deviation of Traction Magnitude from Surroundings')
xlabel({'FTTC reg vs. ',' FastBEM reg vs. ',' FastBEM 2nd reg',' FastBEM 1norm'});

%% save
regParams = [sol_mats_reg.L sol_mats_reg2.L sol_mats_regfine.L sol_mats_reg2fine.L];

hgexport(h1,strcat(imgPath,'fig3 TFM accuracy with finer mesh'),hgexport('factorystyle'),'Format','tiff')
hgsave(h1,strcat(imgPath,'fig3 TFM accuracy with finer mesh'),'-v7.3')
% delete(h1)

save([dataPath '/all data.mat'],'DTMS','DTMSerr','DTA','DTAerr','pSR','pSRerr','flocMax','regParams','-v7.3');
