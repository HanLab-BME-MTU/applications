%% calculateForceFineMeshImgScaled_fn
function [] = testForceReconstruction(forceType,percentNoise,savePath,forceMeshFastBEM,M_FastBEM,basisClassTablePathFine)
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
meshPtsFwdSol=2^10;
% L=0;
% numPoints_out=50;


imgPath = [savePath '/img/'];
dataPath = [savePath '/data/'];
if (~exist(imgPath,'dir') || ~exist(dataPath,'dir'))
    mkdir(imgPath);
    mkdir(dataPath);
end
cd(imgPath);
%% Mesh generation and artificial force generation
gridSpacing = 1;
xmin = gridSpacing;
ymin = gridSpacing;
gridSpacingf = 4;
xmax=gridSpacingf*2*(floor(500/gridSpacingf/2));
ymax=gridSpacingf*2*(floor(500/gridSpacingf/2));
xminf = round(gridSpacingf/2);
yminf = round(gridSpacingf/2);

[x_mat_f, y_mat_f]=meshgrid(xminf:gridSpacingf:xmax,yminf:gridSpacingf:ymax);

[x_mat_u, y_mat_u]=meshgrid(xmin:gridSpacing:xmax,ymin:gridSpacing:ymax);
x_vec_u = reshape(x_mat_u,[],1);
y_vec_u = reshape(y_mat_u,[],1);
% % Hexagonal force mesh
% displField.pos = [x_vec_u y_vec_u];
% [xvec,yvec]=createHexGridFromDisplField(displField,1);

% assumes um/pixel = 0.1, i.e. 10 pixel = 1 um.
% I need force distribution with multiple sources, unit: Pa -Sangyoon 013113
% I need to create 4 types of forces (large FA (20 pix), large force (2kPa) - group1
%                                     large FA, small force (500-1000 Pa) - group2
%                                     small FA (5 pix), large force - group3
%                                     small FA, small force - group4)
%% Forward solution
display('forward solution for fine solution...')
tic
% % temporary - get the coordinates
% figure, imshow(ones(500))
% hold on
% fposxy = [];
% n = 0;
% % Loop, picking up the points.
% disp('Left mouse button picks points.')
% disp('Right mouse button picks last point.')
% but = 1;
% while but == 1
%     [xi,yi,but] = ginput(1);
%     plot(xi,yi,'ro')
%     n = n+1;
%     text(xi,yi-8,num2str(n));
%     fposxy(n,:) = [xi yi];
% end
% hold off
%% positions for NA for 2nd geometry
posNA = [  192.0000  472.0000
  201.0000  450.0000
  212.0000  431.0000
  223.0000  410.0000
  233.0000  391.0000
  245.0000  378.0000
  259.0000  360.0000
  271.0000  342.0000
  285.0000  327.0000
  300.0000  309.0000
  323.0000  292.0000
  342.0000  283.0000
  359.0000  274.0000
  376.0000  260.0000
  388.0000  249.0000
  403.0000  245.0000
  415.0000  236.0000
  428.0000  236.0000
  444.0000  224.0000
  466.0000  220.0000];
%% positions for FA for 2nd geometry
% focal adhesion
% % temporary - get the coordinates
% figure, imshow(ones(500))
% hold on
% plot(posNA(:,1),posNA(:,2),'bo')
% fposxy = [];
% n = 0;
% % Loop, picking up the points.
% disp('Left mouse button picks points.')
% disp('Right mouse button picks last point.')
% but = 1;
% while but == 1
%     [xi,yi,but] = ginput(1);
%     plot(xi,yi,'ro')
%     n = n+1;
%     text(xi,yi-8,num2str(n));
%     fposxy(n,:) = [xi yi];
% end
% hold off
posFA = [265.0000  430.0000
  276.0000  411.0000
  329.0000  338.0000
  342.0000  335.0000
  413.0000  275.0000
  426.0000  272.0000
  330.0000  403.0000
  404.0000  340.0000];
%% forcetype1;
[ux, uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,...
    @(x,y) assumedForceAniso2D(1,x,y,139,267,150,420,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x,y,156,232,110*2,450*2,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x,y,184,212,60,500,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x,y,217,200,20,520,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x,y,246,195,0,550*2,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x,y,272,199,-30*2,510*2,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x,y,297,211,-75,495,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x,y,323,231,-100*2,440*2,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x,y,346,266,-140,400,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x,y,42,323,480,1600,500/108,2000/108,forceType)+...
    assumedForceAniso2D(1,x,y,54,319,400,1660,500/108,2100/108,forceType)+...
    assumedForceAniso2D(1,x,y,96,298,320,1720,500/108,2000/108,forceType)+...
    assumedForceAniso2D(1,x,y,139,288,240,1780,500/108,2300/108,forceType)+...
    assumedForceAniso2D(1,x,y,180,280,160,1840,500/108,2000/108,forceType)+...
    assumedForceAniso2D(1,x,y,225,276,80,1900,500/108,2400/108,forceType)+...
    assumedForceAniso2D(1,x,y,263,273,0,2000,550/108,2300/108,forceType)+...
    assumedForceAniso2D(1,x,y,301,275,-80,1900,600/108,2600/108,forceType)+...
    assumedForceAniso2D(1,x,y,331,281,-160,1840,500/108,2200/108,forceType)+...
    assumedForceAniso2D(1,x,y,347,279,-240,1780,400/108,2000/108,forceType)+...
    assumedForceAniso2D(1,x,y,381,290,-320,1720,600/108,2700/108,forceType)+...
    assumedForceAniso2D(1,x,y,417,303,-400,1660,500/108,2300/108,forceType)+...
    assumedForceAniso2D(1,x,y,455,317,-480,1600,450/108,2100/108,forceType)+...
    assumedForceAniso2D(1,x,y,69,359,600,2500,700/108,3000/108,forceType)+...
    assumedForceAniso2D(1,x,y,91,349,500,2600,600/108,2600/108,forceType)+...
    assumedForceAniso2D(1,x,y,115,337,400,2700,700/108,3500/108,forceType)+...
    assumedForceAniso2D(1,x,y,129,338,300,2800,700/108,3000/108,forceType)+...
    assumedForceAniso2D(1,x,y,168,327,200,2900,600/108,2600/108,forceType)+...
    assumedForceAniso2D(1,x,y,186,321,100,3000,600/108,3000/108,forceType)+...
    assumedForceAniso2D(1,x,y,270,318,0,3100,650/108,3500/108,forceType)+...
    assumedForceAniso2D(1,x,y,288,320,-100,3000,600/108,2600/108,forceType)+...
    assumedForceAniso2D(1,x,y,312,323,-200,2900,600/108,3000/108,forceType)+...
    assumedForceAniso2D(1,x,y,361,330,-300,2800,600/108,2600/108,forceType)+...
    assumedForceAniso2D(1,x,y,383,337,-400,2700,600/108,3500/108,forceType)+...
    assumedForceAniso2D(1,x,y,423,352,-500,2600,600/108,3000/108,forceType)+...
    assumedForceAniso2D(1,x,y,432,363,-600,2500,600/108,2600/108,forceType),...
    @(x,y) assumedForceAniso2D(2,x,y,139,267,150,420,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x,y,156,232,110*2,450*2,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x,y,184,212,60,500,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x,y,217,200,20,520,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x,y,246,195,0,550*2,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x,y,272,199,-30*2,510*2,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x,y,297,211,-75,495,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x,y,323,231,-100*2,440*2,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x,y,346,266,-140,400,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x,y,42,323,480,1600,500/108,2000/108,forceType)+...
    assumedForceAniso2D(2,x,y,54,319,400,1660,500/108,2100/108,forceType)+...
    assumedForceAniso2D(2,x,y,96,298,320,1720,500/108,2000/108,forceType)+...
    assumedForceAniso2D(2,x,y,139,288,240,1780,500/108,2300/108,forceType)+...
    assumedForceAniso2D(2,x,y,180,280,160,1840,500/108,2000/108,forceType)+...
    assumedForceAniso2D(2,x,y,225,276,80,1900,500/108,2400/108,forceType)+...
    assumedForceAniso2D(2,x,y,263,273,0,2000,550/108,2300/108,forceType)+...
    assumedForceAniso2D(2,x,y,301,275,-80,1900,600/108,2600/108,forceType)+...
    assumedForceAniso2D(2,x,y,331,281,-160,1840,500/108,2200/108,forceType)+...
    assumedForceAniso2D(2,x,y,347,279,-240,1780,400/108,2000/108,forceType)+...
    assumedForceAniso2D(2,x,y,381,290,-320,1720,600/108,2700/108,forceType)+...
    assumedForceAniso2D(2,x,y,417,303,-400,1660,500/108,2300/108,forceType)+...
    assumedForceAniso2D(2,x,y,455,317,-480,1600,450/108,2100/108,forceType)+...
    assumedForceAniso2D(2,x,y,69,359,600,2500,700/108,3000/108,forceType)+...
    assumedForceAniso2D(2,x,y,91,349,500,2600,600/108,2600/108,forceType)+...
    assumedForceAniso2D(2,x,y,115,337,400,2700,700/108,3500/108,forceType)+...
    assumedForceAniso2D(2,x,y,129,338,300,2800,700/108,3000/108,forceType)+...
    assumedForceAniso2D(2,x,y,168,327,200,2900,600/108,2600/108,forceType)+...
    assumedForceAniso2D(2,x,y,186,321,100,3000,600/108,3000/108,forceType)+...
    assumedForceAniso2D(2,x,y,270,318,0,3100,650/108,3500/108,forceType)+...
    assumedForceAniso2D(2,x,y,288,320,-100,3000,600/108,2600/108,forceType)+...
    assumedForceAniso2D(2,x,y,312,323,-200,2900,600/108,3000/108,forceType)+...
    assumedForceAniso2D(2,x,y,361,330,-300,2800,600/108,2600/108,forceType)+...
    assumedForceAniso2D(2,x,y,383,337,-400,2700,600/108,3500/108,forceType)+...
    assumedForceAniso2D(2,x,y,423,352,-500,2600,600/108,3000/108,forceType)+...
    assumedForceAniso2D(2,x,y,432,363,-600,2500,600/108,2600/108,forceType),'fft',[],meshPtsFwdSol); %,'conv',[],meshPtsFwdSol);
force_x = assumedForceAniso2D(1,x_mat_u,y_mat_u,139,267,150,420,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,156,232,110*2,450*2,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,184,212,60,500,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,217,200,20,520,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,246,195,0,550*2,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,272,199,-30*2,510*2,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,297,211,-75,495,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,323,231,-100*2,440*2,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,346,266,-140,400,400/72,500/72,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,42,323,480,1600,500/108,2000/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,54,319,400,1660,500/108,2100/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,96,298,320,1720,500/108,2000/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,139,288,240,1780,500/108,2300/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,180,280,160,1840,500/108,2000/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,225,276,80,1900,500/108,2400/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,263,273,0,2000,550/108,2300/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,301,275,-80,1900,600/108,2600/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,331,281,-160,1840,500/108,2200/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,347,279,-240,1780,400/108,2000/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,381,290,-320,1720,600/108,2700/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,417,303,-400,1660,500/108,2300/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,455,317,-480,1600,450/108,2100/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,69,359,600,2500,700/108,3000/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,91,349,500,2600,600/108,2600/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,115,337,400,2700,700/108,3500/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,129,338,300,2800,700/108,3000/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,168,327,200,2900,600/108,2600/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,186,321,100,3000,600/108,3000/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,270,318,0,3100,650/108,3500/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,288,320,-100,3000,600/108,2600/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,312,323,-200,2900,600/108,3000/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,361,330,-300,2800,600/108,2600/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,383,337,-400,2700,600/108,3500/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,423,352,-500,2600,600/108,3000/108,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,432,363,-600,2500,600/108,2600/108,forceType);
force_y = assumedForceAniso2D(2,x_mat_u,y_mat_u,139,267,150,420,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,156,232,110*2,450*2,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,184,212,60,500,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,217,200,20,520,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,246,195,0,550*2,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,272,199,-30*2,510*2,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,297,211,-75,495,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,323,231,-100*2,440*2,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,346,266,-140,400,400/72,500/72,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,42,323,480,1600,500/108,2000/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,54,319,400,1660,500/108,2100/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,96,298,320,1720,500/108,2000/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,139,288,240,1780,500/108,2300/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,180,280,160,1840,500/108,2000/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,225,276,80,1900,500/108,2400/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,263,273,0,2000,550/108,2300/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,301,275,-80,1900,600/108,2600/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,331,281,-160,1840,500/108,2200/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,347,279,-240,1780,400/108,2000/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,381,290,-320,1720,600/108,2700/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,417,303,-400,1660,500/108,2300/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,455,317,-480,1600,450/108,2100/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,69,359,600,2500,700/108,3000/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,91,349,500,2600,600/108,2600/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,115,337,400,2700,700/108,3500/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,129,338,300,2800,700/108,3000/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,168,327,200,2900,600/108,2600/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,186,321,100,3000,600/108,3000/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,270,318,0,3100,650/108,3500/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,288,320,-100,3000,600/108,2600/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,312,323,-200,2900,600/108,3000/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,361,330,-300,2800,600/108,2600/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,383,337,-400,2700,600/108,3500/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,423,352,-500,2600,600/108,3000/108,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,432,363,-600,2500,600/108,2600/108,forceType);

% % 2nd force composition
% [ux, uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,...
%     @(x,y) assumedForceAniso2D(1,x,y,192,472,450,10,500/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x,y,201,450,760,60,400/72,500/72,forceType)+...
%     assumedForceAniso2D(1,x,y,212,431,600,70,500/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x,y,223,410,540,90,400/72,500/72,forceType)+...
%     assumedForceAniso2D(1,x,y,233,391,930,220,500/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x,y,245,378,670,510,500/72,500/72,forceType)+...
%     assumedForceAniso2D(1,x,y,259,360,620,495,500/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x,y,271,342,560,440,400/72,500/72,forceType)+...
%     assumedForceAniso2D(1,x,y,285,327,490,400,400/72,500/72,forceType)+...
%     assumedForceAniso2D(1,x,y,300,309,480,600,500/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x,y,323,292,450,660,500/72,700/72,forceType)+...
%     assumedForceAniso2D(1,x,y,342,283,410,520,500/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x,y,359,274,390,480,500/72,700/72,forceType)+...
%     assumedForceAniso2D(1,x,y,376,260,340,540,500/72,500/72,forceType)+...
%     assumedForceAniso2D(1,x,y,388,249,300,760,500/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x,y,403,245,280,800,550/72,700/72,forceType)+...
%     assumedForceAniso2D(1,x,y,415,236,240,800,600/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x,y,428,236,170,840,500/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x,y,444,224,90,980,400/72,500/72,forceType)+...
%     assumedForceAniso2D(1,x,y,466,220,10,920,600/72,700/72,forceType)+...
%     assumedForceAniso2D(1,x,y,265,430,1660,160,500/72,2300/72,forceType)+...
%     assumedForceAniso2D(1,x,y,276,411,1580,200,450/72,2100/72,forceType)+...
%     assumedForceAniso2D(1,x,y,329,338,1300,920,700/72,3000/72,forceType)+...
%     assumedForceAniso2D(1,x,y,342,335,1200,1000,600/72,2600/72,forceType)+...
%     assumedForceAniso2D(1,x,y,413,275,300,2700,700/72,3500/72,forceType)+...
%     assumedForceAniso2D(1,x,y,426,272,200,2800,700/72,3000/72,forceType)+...
%     assumedForceAniso2D(1,x,y,330,403,2100,200,600/72,2600/72,forceType)+...
%     assumedForceAniso2D(1,x,y,404,340,1500,1300,600/72,3000/72,forceType),...
%     @(x,y) assumedForceAniso2D(2,x,y,192,472,450,10,500/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x,y,201,450,760,60,400/72,500/72,forceType)+...
%     assumedForceAniso2D(2,x,y,212,431,600,70,500/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x,y,223,410,540,90,400/72,500/72,forceType)+...
%     assumedForceAniso2D(2,x,y,233,391,930,220,500/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x,y,245,378,670,510,500/72,500/72,forceType)+...
%     assumedForceAniso2D(2,x,y,259,360,620,495,500/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x,y,271,342,560,440,400/72,500/72,forceType)+...
%     assumedForceAniso2D(2,x,y,285,327,490,400,400/72,500/72,forceType)+...
%     assumedForceAniso2D(2,x,y,300,309,480,600,500/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x,y,323,292,450,660,500/72,700/72,forceType)+...
%     assumedForceAniso2D(2,x,y,342,283,410,520,500/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x,y,359,274,390,480,500/72,700/72,forceType)+...
%     assumedForceAniso2D(2,x,y,376,260,340,540,500/72,500/72,forceType)+...
%     assumedForceAniso2D(2,x,y,388,249,300,760,500/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x,y,403,245,280,800,550/72,700/72,forceType)+...
%     assumedForceAniso2D(2,x,y,415,236,240,800,600/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x,y,428,236,170,840,500/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x,y,444,224,90,980,400/72,500/72,forceType)+...
%     assumedForceAniso2D(2,x,y,466,220,10,920,600/72,700/72,forceType)+...
%     assumedForceAniso2D(2,x,y,265,430,1660,160,500/72,2300/72,forceType)+...
%     assumedForceAniso2D(2,x,y,276,411,1580,200,450/72,2100/72,forceType)+...
%     assumedForceAniso2D(2,x,y,329,338,1300,920,700/72,3000/72,forceType)+...
%     assumedForceAniso2D(2,x,y,342,335,1200,1000,600/72,2600/72,forceType)+...
%     assumedForceAniso2D(2,x,y,413,275,300,2700,700/72,3500/72,forceType)+...
%     assumedForceAniso2D(2,x,y,426,272,200,2800,700/72,3000/72,forceType)+...
%     assumedForceAniso2D(2,x,y,330,403,2100,200,600/72,2600/72,forceType)+...
%     assumedForceAniso2D(2,x,y,404,340,1500,1300,600/72,3000/72,forceType),'fft',[],meshPtsFwdSol); %,'conv',[],meshPtsFwdSol);
% force_x = assumedForceAniso2D(1,x_mat_u,y_mat_u,192,472,450,10,500/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,201,450,760,60,400/72,500/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,212,431,600,70,500/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,223,410,540,90,400/72,500/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,233,391,930,220,500/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,245,378,670,510,500/72,500/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,259,360,620,495,500/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,271,342,560,440,400/72,500/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,285,327,490,400,400/72,500/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,300,309,480,600,500/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,323,292,450,660,500/72,700/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,342,283,410,520,500/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,359,274,390,480,500/72,700/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,376,260,340,540,500/72,500/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,388,249,300,760,500/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,403,245,280,800,550/72,700/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,415,236,240,800,600/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,428,236,170,840,500/72,600/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,444,224,90,980,400/72,500/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,466,220,10,920,600/72,700/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,265,430,1660,160,500/72,2300/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,276,411,1580,200,450/72,2100/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,329,338,1300,920,700/72,3000/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,342,335,1200,1000,600/72,2600/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,413,275,300,2700,700/72,3500/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,426,272,200,2800,700/72,3000/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,330,403,2100,200,600/72,2600/72,forceType)+...
%     assumedForceAniso2D(1,x_mat_u,y_mat_u,404,340,1500,1300,600/72,3000/72,forceType);
% force_y = assumedForceAniso2D(2,x_mat_u,y_mat_u,192,472,450,10,500/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,201,450,760,60,400/72,500/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,212,431,600,70,500/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,223,410,540,90,400/72,500/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,233,391,930,220,500/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,245,378,670,510,500/72,500/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,259,360,620,495,500/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,271,342,560,440,400/72,500/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,285,327,490,400,400/72,500/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,300,309,480,600,500/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,323,292,450,660,500/72,700/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,342,283,410,520,500/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,359,274,390,480,500/72,700/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,376,260,340,540,500/72,500/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,388,249,300,760,500/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,403,245,280,800,550/72,700/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,415,236,240,800,600/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,428,236,170,840,500/72,600/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,444,224,90,980,400/72,500/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,466,220,10,920,600/72,700/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,265,430,1660,160,500/72,2300/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,276,411,1580,200,450/72,2100/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,329,338,1300,920,700/72,3000/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,342,335,1200,1000,600/72,2600/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,413,275,300,2700,700/72,3500/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,426,272,200,2800,700/72,3000/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,330,403,2100,200,600/72,2600/72,forceType)+...
%     assumedForceAniso2D(2,x_mat_u,y_mat_u,404,340,1500,1300,600/72,3000/72,forceType);

toc

display(['adding noise (' num2str(percentNoise*100) '%) ...'])
if addNoise==1
    max_u=max([max(max(ux)) max(max(uy))]); 
    noise_xu=normrnd(0,percentNoise*max_u,floor((xmax-xmin+1)/gridSpacing),floor((ymax-ymin+1)/gridSpacing));
    noise_yu=normrnd(0,percentNoise*max_u,floor((xmax-xmin+1)/gridSpacing),floor((ymax-ymin+1)/gridSpacing));
    ux=ux+noise_xu;
    uy=uy+noise_yu;
end
% u_mat(:,:,1)=ux;
% u_mat(:,:,2)=uy;
ux_vec=reshape(ux,[],1);
uy_vec=reshape(uy,[],1);
% u=vertcat(ux_vec,uy_vec);

%% grid setting
pix_durch_my=1; %seem to be important only for the energy
grid_mat_f(:,:,1)=x_mat_f; %dense
grid_mat_f(:,:,2)=y_mat_f; %dense
[grid_mat_f,u_mat, i_max,j_max] = interp_vec2grid([x_vec_u y_vec_u], [ux_vec uy_vec],[],grid_mat_f);
x_vec_f = reshape(x_mat_f,[],1);
y_vec_f = reshape(y_mat_f,[],1);
ux_vec_f = reshape(u_mat(:,:,1),[],1);
uy_vec_f = reshape(u_mat(:,:,2),[],1);
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
L= 5e-6;
% i_max = size(grid_mat_f,1);
% j_max = size(grid_mat_f,2);
cluster_size = grid_mat_f(1,1,1) - grid_mat_f(2,2,1);

[~,~,force_FTTC, ~] = reg_fourier_TFM(grid_mat_f,u_mat,E,s, pix_durch_my, cluster_size, i_max, j_max, L);  
fttcL = L;
rec_force_FTTC(:,:,1)=reshape(force_FTTC(:,1),i_max,j_max);
rec_force_FTTC(:,:,2)=reshape(force_FTTC(:,2),i_max,j_max);
%% FastBEM with fine mesh (For figure 3) 
% Mesh creation
if isempty(forceMeshFastBEM)
    display('creating fine mesh and basis function (1x1)')
    tic
    forceMeshFastBEM=createMeshAndBasisFastBEM(x_vec_f,y_vec_f,false,[],0);
    toc
end
L = 1e-6;
display(['FastBEM (' num2str(L) ')...'])
[fx_FastBEM,fy_FastBEM,~,~,M_FastBEM,~,~,~,sol_mats_BEM]...
    = BEM_force_reconstruction(x_mat_f,y_mat_f,u_mat(:,:,1),u_mat(:,:,2),forceMeshFastBEM,...
    E,L,[],[],'fast',meshPtsFwdSol,...
    'QR','basisClassTblPath',basisClassTablePathFine,'imgRows',ymax,'imgCols',xmax,'fwdMap',M_FastBEM);
rec_force_FastBEMreg(:,:,1)=fx_FastBEM;
rec_force_FastBEMreg(:,:,2)=fy_FastBEM;
%% FastBEM-reconstruction with laplacian regularization with fine mesh
L = 1e-6; %sol_mats_reg.L;
display(['FastBEM with 2nd order regularization (' num2str(L) ')...'])
tic
[fx_FastBEMreg2,fy_FastBEMreg2,~,~,~,~,~,~,sol_mats_reg2]...
    = BEM_force_reconstruction(x_mat_f,y_mat_f,u_mat(:,:,1),u_mat(:,:,2),forceMeshFastBEM,...
    E,L,[],[],'fast',meshPtsFwdSol,...
    'LaplacianReg','basisClassTblPath',basisClassTablePathFine,...
    'imgRows',ymax,'imgCols',xmax,'fwdMap',M_FastBEM);
toc
rec_force_FastBEMreg2(:,:,1)=fx_FastBEMreg2;
rec_force_FastBEMreg2(:,:,2)=fy_FastBEMreg2;

%% FastBEM-reconstruction with 1Norm-0th order regularization with fine mesh
L = 1e-3; %sol_mats_reg.L;
display(['FastBEM with 1Norm-0th order regularization with fine mesh (' num2str(L) ')...'])
tic
[fx_FastBEM1n,fy_FastBEM1n,~,~,~,~,~,~,sol_mats_1n]...
    = BEM_force_reconstruction(x_mat_f,y_mat_f,u_mat(:,:,1),u_mat(:,:,2),forceMeshFastBEM,...
    E,L,[],[],'fast',meshPtsFwdSol,...
    '1NormReg','basisClassTblPath',basisClassTablePathFine,...
    'imgRows',ymax,'imgCols',xmax,'fwdMap',M_FastBEM);
toc
rec_force_FastBEM1n(:,:,1)=fx_FastBEM1n;
rec_force_FastBEM1n(:,:,2)=fy_FastBEM1n;

%% FastBEM-reconstruction with 1Norm-2th order regularization with fine mesh
L = 1e-3; %sol_mats_reg.L;
display(['FastBEM with 1Norm-2nd order regularization with fine mesh (' num2str(L) ')...'])
tic
[fx_FastBEM1n2,fy_FastBEM1n2,~,~,~,~,~,~,sol_mats_1n2]...
    = BEM_force_reconstruction(x_mat_f,y_mat_f,u_mat(:,:,1),u_mat(:,:,2),forceMeshFastBEM,...
    E,L,[],[],'fast',meshPtsFwdSol,...
    '1NormRegLaplacian','basisClassTblPath',basisClassTablePathFine,...
    'imgRows',ymax,'imgCols',xmax,'fwdMap',M_FastBEM);
toc
rec_force_FastBEM1n2(:,:,1)=fx_FastBEM1n2;
rec_force_FastBEM1n2(:,:,2)=fy_FastBEM1n2;

%% quick saving of M and forcemesh
save([dataPath '/MandforceMesh.mat'],'M_FastBEM','forceMeshFastBEM',...
    'rec_force_FTTC','rec_force_FastBEMreg','rec_force_FastBEMreg2','fy_FastBEMreg2','rec_force_FastBEM1n',...
    'rec_force_FastBEM1n2','-v7.3');
%% Heat map presentation for each force magnitude
pos = [x_vec_f y_vec_f]; %sparse
disp_vec = [ux_vec_f uy_vec_f]; %sparse
posu = [x_vec_u y_vec_u]; %dense
disp_vecu = [ux_vec uy_vec]; %dense
grid_mat_u(:,:,1) = x_mat_u;
grid_mat_u(:,:,2) = y_mat_u;

rec_force_org_vec = [reshape(force_x,[],1) reshape(force_y,[],1)];
[~,tmat_org, ~, ~] = interp_vec2grid(posu+disp_vecu, rec_force_org_vec,[],grid_mat_u); %1:cluster size
tnorm_org = (tmat_org(:,:,1).^2 + tmat_org(:,:,2).^2).^0.5; %this should be fine mesh
fmat_org(:,:,1) = force_x;
fmat_org(:,:,2) = force_y;
fnorm_org = (force_x.^2 + force_y.^2).^0.5; %this should be fine mesh

rec_force_FTTCregfine_vec = [reshape(rec_force_FTTC(:,:,1),[],1) reshape(rec_force_FTTC(:,:,2),[],1)];
rec_force_FTTCinterp(:,:,1) = griddata(x_mat_f,y_mat_f,rec_force_FTTC(:,:,1),x_mat_u,y_mat_u,'cubic');
rec_force_FTTCinterp(:,:,2) = griddata(x_mat_f,y_mat_f,rec_force_FTTC(:,:,2),x_mat_u,y_mat_u,'cubic');
[~,tmat_fttcreg, ~, ~] = interp_vec2grid(pos+disp_vec, rec_force_FTTCregfine_vec,[],grid_mat_u); %1:cluster size
tnorm_fttcreg = (tmat_fttcreg(:,:,1).^2 + tmat_fttcreg(:,:,2).^2).^0.5;
fnorm_fttcregc = (rec_force_FTTC(:,:,1).^2 + rec_force_FTTC(:,:,2).^2).^0.5;
fnorm_fttcreg = griddata(x_mat_f,y_mat_f,fnorm_fttcregc,x_mat_u,y_mat_u,'cubic');

rec_force_FastBEMreg_vec = [reshape(rec_force_FastBEMreg(:,:,1),[],1) reshape(rec_force_FastBEMreg(:,:,2),[],1)];
[~,tmat_FastBEMregfine, ~, ~] = interp_vec2grid(pos+disp_vec, rec_force_FastBEMreg_vec,[],grid_mat_u); %1:cluster size
tnorm_FastBEMregfine = (tmat_FastBEMregfine(:,:,1).^2 + tmat_FastBEMregfine(:,:,2).^2).^0.5;
fmat_FastBEMregfine(:,:,1) = reshape(rec_force_FastBEMreg(:,:,1),size(x_mat_f));
fmat_FastBEMregfine(:,:,2) = reshape(rec_force_FastBEMreg(:,:,2),size(y_mat_f));
fnorm_FastBEMregfinec = (fmat_FastBEMregfine(:,:,1).^2 + fmat_FastBEMregfine(:,:,2).^2).^0.5;
fnorm_FastBEMreg = griddata(x_mat_f,y_mat_f,fnorm_FastBEMregfinec,x_mat_u,y_mat_u,'cubic');
fmat_FastBEMreginterp(:,:,1) = griddata(x_mat_f,y_mat_f,fmat_FastBEMregfine(:,:,1),x_mat_u,y_mat_u,'cubic');
fmat_FastBEMreginterp(:,:,2) = griddata(x_mat_f,y_mat_f,fmat_FastBEMregfine(:,:,2),x_mat_u,y_mat_u,'cubic');

rec_force_FastBEMreg2_vec = [reshape(rec_force_FastBEMreg2(:,:,1),[],1) reshape(rec_force_FastBEMreg2(:,:,2),[],1)];
[~,tmat_FastBEMreg2, ~, ~] = interp_vec2grid(pos+disp_vec, rec_force_FastBEMreg2_vec,[],grid_mat_u); %1:cluster size
tnorm_FastBEMreg2fine = (tmat_FastBEMreg2(:,:,1).^2 + tmat_FastBEMreg2(:,:,2).^2).^0.5;
fmat_FastBEMreg2fine(:,:,1) = reshape(rec_force_FastBEMreg2(:,:,1),size(x_mat_f));
fmat_FastBEMreg2fine(:,:,2) = reshape(rec_force_FastBEMreg2(:,:,2),size(y_mat_f));
fnorm_FastBEMreg2finec = (fmat_FastBEMreg2fine(:,:,1).^2 + fmat_FastBEMreg2fine(:,:,2).^2).^0.5;
fnorm_FastBEMreg2 = griddata(x_mat_f,y_mat_f,fnorm_FastBEMreg2finec,x_mat_u,y_mat_u,'cubic');
fmat_FastBEMreg2interp(:,:,1) = griddata(x_mat_f,y_mat_f,fmat_FastBEMreg2fine(:,:,1),x_mat_u,y_mat_u,'cubic');
fmat_FastBEMreg2interp(:,:,2) = griddata(x_mat_f,y_mat_f,fmat_FastBEMreg2fine(:,:,2),x_mat_u,y_mat_u,'cubic');

rec_force_FastBEM1nfine_vec = [reshape(rec_force_FastBEM1n(:,:,1),[],1) reshape(rec_force_FastBEM1n(:,:,2),[],1)];
[~,tmat_FastBEM1nfine, ~, ~] = interp_vec2grid(pos+disp_vec, rec_force_FastBEM1nfine_vec,[],grid_mat_u); %1:cluster size
tnorm_FastBEM1nfine = (tmat_FastBEM1nfine(:,:,1).^2 + tmat_FastBEM1nfine(:,:,2).^2).^0.5;
fmat_FastBEM1nfine(:,:,1) = reshape(rec_force_FastBEM1n(:,:,1),size(x_mat_f));
fmat_FastBEM1nfine(:,:,2) = reshape(rec_force_FastBEM1n(:,:,2),size(y_mat_f));
fnorm_FastBEM1nfinec = (fmat_FastBEM1nfine(:,:,1).^2 + fmat_FastBEM1nfine(:,:,2).^2).^0.5;
fnorm_FastBEM1n = griddata(x_mat_f,y_mat_f,fnorm_FastBEM1nfinec,x_mat_u,y_mat_u,'cubic');
fmat_FastBEMreg1ninterp(:,:,1) = griddata(x_mat_f,y_mat_f,fmat_FastBEM1nfine(:,:,1),x_mat_u,y_mat_u,'cubic');
fmat_FastBEMreg1ninterp(:,:,2) = griddata(x_mat_f,y_mat_f,fmat_FastBEM1nfine(:,:,2),x_mat_u,y_mat_u,'cubic');

rec_force_FastBEM1n2fine_vec = [reshape(rec_force_FastBEM1n2(:,:,1),[],1) reshape(rec_force_FastBEM1n2(:,:,2),[],1)];
[~,tmat_FastBEM1n2fine, ~, ~] = interp_vec2grid(pos+disp_vec, rec_force_FastBEM1n2fine_vec,[],grid_mat_u); %1:cluster size
tnorm_FastBEM1n2fine = (tmat_FastBEM1n2fine(:,:,1).^2 + tmat_FastBEM1n2fine(:,:,2).^2).^0.5;
fmat_FastBEM1n2fine(:,:,1) = reshape(rec_force_FastBEM1n2(:,:,1),size(x_mat_f));
fmat_FastBEM1n2fine(:,:,2) = reshape(rec_force_FastBEM1n2(:,:,2),size(y_mat_f));
fnorm_FastBEM1n2finec = (fmat_FastBEM1n2fine(:,:,1).^2 + fmat_FastBEM1n2fine(:,:,2).^2).^0.5;
fnorm_FastBEM1n2 = griddata(x_mat_f,y_mat_f,fnorm_FastBEM1n2finec,x_mat_u,y_mat_u,'cubic');
fmat_FastBEMreg1n2interp(:,:,1) = griddata(x_mat_f,y_mat_f,fmat_FastBEM1n2fine(:,:,1),x_mat_u,y_mat_u,'cubic');
fmat_FastBEMreg1n2interp(:,:,2) = griddata(x_mat_f,y_mat_f,fmat_FastBEM1n2fine(:,:,2),x_mat_u,y_mat_u,'cubic');

%% Force at cross-sectional line for large FA, high peak
% For force 1
rangeP = 20;
tnorm1D_org_g1 = fnorm_org(273-rangeP:273+rangeP,263);
tnorm1D_FTTCregfine_g1 = fnorm_fttcreg(273-rangeP:273+rangeP,263);
tnorm1D_FastBEMregfine_g1 = fnorm_FastBEMreg(273-rangeP:273+rangeP,263);
tnorm1D_FastBEMreg2fine_g1 = fnorm_FastBEMreg2(273-rangeP:273+rangeP,263);
tnorm1D_FastBEM1nfine_g1 = fnorm_FastBEM1n(273-rangeP:273+rangeP,263);
tnorm1D_FastBEM1n2fine_g1 = fnorm_FastBEM1n2(273-rangeP:273+rangeP,263);
y1Du_g1 = grid_mat_u(273-rangeP:273+rangeP,263,2);
%% Force at cross-sectional line for NA, low peak
% For force 1
rangeP = 20;
tnorm1D_org_g3 = fnorm_org(195-rangeP:195+rangeP,246);
tnorm1D_FTTCregfine_g3 = fnorm_fttcreg(195-rangeP:195+rangeP,246);
tnorm1D_FastBEMregfine_g3 = fnorm_FastBEMreg(195-rangeP:195+rangeP,246);
tnorm1D_FastBEMreg2fine_g3 = fnorm_FastBEMreg2(195-rangeP:195+rangeP,246);
tnorm1D_FastBEM1nfine_g3 = fnorm_FastBEM1n(195-rangeP:195+rangeP,246);
tnorm1D_FastBEM1n2fine_g3 = fnorm_FastBEM1n2(195-rangeP:195+rangeP,246);
y1Du_g3 = grid_mat_u(195-rangeP:195+rangeP,246,2);
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
tnormd(:,:,1) = fnorm_fttcreg;
tnormd(:,:,2) = fnorm_FastBEMreg;
tnormd(:,:,3) = fnorm_FastBEMreg2;
tnormd(:,:,4) = fnorm_FastBEM1n;
tnormd(:,:,5) = fnorm_FastBEM1n2;
nMethods = 5;

% zeroI = [1 1]; %index of zero coordinates
% [~,locMaxI,~] = findMaxScoreI(tnorm_org,zeroI,1,0); %from dense data
% posFA = [   42.0000  323.0000
%    54.0000  319.0000
%    96.0000  298.0000
%   139.0000  288.0000
%   180.0000  280.0000
%   225.0000  276.0000
%   263.0000  273.0000
%   301.0000  275.0000
%   331.0000  281.0000
%   347.0000  279.0000
%   381.0000  290.0000
%   417.0000  303.0000
%   455.0000  317.0000
%    69.0000  359.0000
%    91.0000  349.0000
%   115.0000  337.0000
%   129.0000  338.0000
%   168.0000  327.0000
%   186.0000  321.0000
%   270.0000  318.0000
%   288.0000  320.0000
%   312.0000  323.0000
%   361.0000  330.0000
%   383.0000  337.0000
%   423.0000  352.0000
%   432.0000  363.0000];
posFA = [139.0000  267.0000
  156.0000  232.0000
  184.0000  212.0000
  217.0000  200.0000
  246.0000  195.0000
  272.0000  199.0000
  297.0000  211.0000
  323.0000  231.0000
  346.0000  266.0000];
locMaxI = posFA(:,2:-1:1);
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
% r = 1.1; %
backgroundMaskfine = fnorm_org==0;

for jj=1:nMethods 
    for k=1:nPeaks
        flocMax(k,jj) = tnormd(locMaxI(k,1),locMaxI(k,2),jj);
        flocMaxRatio(k,jj) = flocMax(k,jj)/fnorm_org(locMaxI(k,1),locMaxI(k,2)); % ratio based on original force norm (0~1)
        if jj==1
            recF_vecfine = [rec_force_FTTCinterp(locMaxI(k,1),locMaxI(k,2),1),...
                rec_force_FTTCinterp(locMaxI(k,1),locMaxI(k,2),2)];
        elseif jj==2
            recF_vecfine = [fmat_FastBEMreginterp(locMaxI(k,1)+1,locMaxI(k,2),1),...
                fmat_FastBEMreginterp(locMaxI(k,1),locMaxI(k,2),2)];
        elseif jj==3
            recF_vecfine = [fmat_FastBEMreg2interp(locMaxI(k,1),locMaxI(k,2),1),...
                fmat_FastBEMreg2interp(locMaxI(k,1),locMaxI(k,2),2)];
        elseif jj==4
            recF_vecfine = [fmat_FastBEMreg1ninterp(locMaxI(k,1),locMaxI(k,2),1),...
                fmat_FastBEMreg1ninterp(locMaxI(k,1),locMaxI(k,2),2)];
        elseif jj==5
            recF_vecfine = [fmat_FastBEMreg1n2interp(locMaxI(k,1),locMaxI(k,2),1),...
                fmat_FastBEMreg1n2interp(locMaxI(k,1),locMaxI(k,2),2)];
        end
        orgF_vec = [fmat_org(locMaxI(k,1),locMaxI(k,2),1),...
                    fmat_org(locMaxI(k,1),locMaxI(k,2),2)];

        flocMaxAngle(k,jj) = acosd(recF_vecfine*orgF_vec'/...
            (norm(recF_vecfine)*norm(orgF_vec))); % in degree
        %indices for surroundings (2 um radius)
        % find a label of the according force
        backgroundMaskLabel = bwlabel(~backgroundMaskfine);
        iLabel = backgroundMaskLabel(locMaxI(k,1),locMaxI(k,2));
        currentMaxMask = backgroundMaskLabel==iLabel;
        ringMaskfine =  bwdist(currentMaxMask)>0 ...
                    & bwdist(currentMaxMask)<2 ...
                    & backgroundMaskfine;
        neighborTnormfine = ringMaskfine.*tnormd(:,:,jj);
        meanNeiTnormfine = nansum(neighborTnormfine(:))/sum(ringMaskfine(:));
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
forceScale=.5*sqrt(max(max(force_x.^2+force_y.^2)));
hold on
quiver(x_mat_u,y_mat_u,force_x/forceScale,force_y/forceScale,0,'k');
quiver(grid_mat_f(:,:,1),grid_mat_f(:,:,2),rec_force_FTTC(:,:,1)/forceScale,rec_force_FTTC(:,:,2)/forceScale,0,'g');
quiver(grid_mat_f(:,:,1),grid_mat_f(:,:,2),fmat_FastBEMregfine(:,:,1)/forceScale,fmat_FastBEMregfine(:,:,2)/forceScale,0,'b');
quiver(grid_mat_f(:,:,1),grid_mat_f(:,:,2),fmat_FastBEMreg2fine(:,:,1)/forceScale,fmat_FastBEMreg2fine(:,:,2)/forceScale,0,'r');
quiver(grid_mat_f(:,:,1),grid_mat_f(:,:,2),fmat_FastBEM1nfine(:,:,1)/forceScale,fmat_FastBEM1nfine(:,:,2)/forceScale,0,'c');
quiver(grid_mat_f(:,:,1),grid_mat_f(:,:,2),fmat_FastBEM1n2fine(:,:,1)/forceScale,fmat_FastBEM1n2fine(:,:,2)/forceScale,0,'m');

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
tmin = min(min(tnorm_org));
tmax = max(max(tnorm_org))*1.1;
colormap('jet');
subplot(3,4,1), imshow(tnorm_org,[tmin-1 tmax+1]); colormap jet;
% surf(grid_matfine(:,:,1), grid_matfine(:,:,2), tnorm_org,'FaceColor','interp',...
% 	'EdgeColor','none', 'FaceLighting','phong'), zlim([tmin-1 tmax+1]), ...
%     view(0,90), 
% caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
title('Original force')
% % hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
% xlim([xmin xmax])
% ylim([ymin ymax])

subplot(3,4,2), imshow(tnorm_fttcreg,[tmin tmax]); 
title(['FTTC regularization L=' num2str(fttcL)])
% caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
% hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
% xlim([xmin xmax])
% ylim([ymin ymax])

subplot(3,4,3), imshow(tnorm_FastBEMregfine,[tmin tmax]);  
title(['FastBEM regularization L=' num2str(sol_mats_BEM.L)])
% caxis([tmin-0.01 tmax+0.01]),set(gca, 'DataAspectRatio', [1,1,1]);
% hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
% xlim([xmin xmax])
% ylim([ymin ymax])

subplot(3,4,5),  imshow(tnorm_FastBEMreg2fine,[tmin tmax]); 
% view(0,90), shading interp; 
title(['FastBEM 2nd order regularization L=' num2str(sol_mats_reg2.L)] )
% set(gca, 'DataAspectRatio', [1,1,1]);
% hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
% xlim([xmin xmax])
% ylim([ymin ymax])
% caxis([tmin-0.01 tmax+0.01]);

subplot(3,4,6), imshow(tnorm_FastBEM1nfine,[tmin tmax]);
% view(0,90), shading interp; 
title(['FastBEM 1norm regularization L=' num2str(sol_mats_1n.L)])
% set(gca, 'DataAspectRatio', [1,1,1]);
% hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
% xlim([xmin xmax])
% ylim([ymin ymax])
% caxis([tmin-0.01 tmax+0.01]);

subplot(3,4,7), imshow(tnorm_FastBEM1n2fine,[tmin tmax]);colormap jet;
% view(0,90), shading interp; 
title(['FastBEM 1norm Laplacian L=' num2str(sol_mats_1n2.L)])
% set(gca, 'DataAspectRatio', [1,1,1]);
% hold on, line([4 16],[10 10],[tmax+1 tmax+1],'Color','y'), line([16 28],[15 15],[tmax+1 tmax+1],'Color','c')
% xlim([xmin xmax])
% ylim([ymin ymax])
% caxis([tmin-0.01 tmax+0.01]);

subplot(3,4,8), caxis([tmin-0.01 tmax+0.01]), axis off; colorbar('West')
%% profile plot
subplot(3,4,9), %for high peak
plot(y1Du_g1,tnorm1D_org_g1,'k','Linewidth',2)
hold on
plot(y1Du_g1,tnorm1D_FTTCregfine_g1,'c','Linewidth',2)
plot(y1Du_g1,tnorm1D_FastBEMregfine_g1,'g','Linewidth',2)
plot(y1Du_g1,tnorm1D_FastBEMreg2fine_g1,'r','Linewidth',2)
plot(y1Du_g1,tnorm1D_FastBEM1nfine_g1,'b','Linewidth',2)
plot(y1Du_g1,tnorm1D_FastBEM1n2fine_g1,'m','Linewidth',2)
xlabel('x')
ylabel('Traction Stress (Pa)')
% ymax1 = max(tnorm1D_org_g1)*3;
% ylim([0 ymax1]);
xlim([y1Du_g1(1) y1Du_g1(end)]);
title({'1D stress profile for','large FA, high peak'})
%% profile plot for peak from small area
subplot(3,4,10), %for high peak
plot(y1Du_g3,tnorm1D_org_g3,'k','Linewidth',2)
hold on
plot(y1Du_g3,tnorm1D_FTTCregfine_g3,'c','Linewidth',2)
plot(y1Du_g3,tnorm1D_FastBEMregfine_g3,'g','Linewidth',2)
plot(y1Du_g3,tnorm1D_FastBEMreg2fine_g3,'r','Linewidth',2)
plot(y1Du_g3,tnorm1D_FastBEM1nfine_g3,'b','Linewidth',2)
plot(y1Du_g3,tnorm1D_FastBEM1n2fine_g3,'m','Linewidth',2)
xlabel('x')
ylabel('Traction Stress (Pa)')
% ymax1 = max(tnorm1D_org_g3)*3;
% ylim([0 ymax1]);
xlim([y1Du_g3(1) y1Du_g3(end)]);
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
subplot(3,4,11),
bar(pSR), hold on
errorbar(pSR,pSRerr)
title('Peak Stress ratio')
xlabel({'FTTC reg vs. ',' FastBEM reg vs. ',' FastBEM 2nd reg',' FastBEM 1norm',' 1norm 2nd'});
%% Deviation of traction magnitude surroundings (DTMS)
subplot(3,4,12),
bar(DTMS), hold on
errorbar(DTMS,DTMSerr)
title('Deviation of Traction Magnitude from Surroundings')
xlabel({'FTTC reg vs. ',' FastBEM reg vs. ',' FastBEM 2nd reg',' FastBEM 1norm',' 1norm 2nd'});

%% save
hgexport(h1,strcat(imgPath,'fig1 TFM accuracy with finer mesh'),hgexport('factorystyle'),'Format','tiff')
hgsave(h1,strcat(imgPath,'fig1 TFM accuracy with finer mesh'),'-v7.3')
print(h1,[imgPath '/fig1.eps'],'-depsc')
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
solTools = {sol_mats_BEM.tool, sol_mats_reg2.tool, sol_mats_1n.tool, sol_mats_1n2.tool};
regParams = [sol_mats_BEM.L sol_mats_reg2.L sol_mats_1n.L sol_mats_1n2.L];
save([dataPath '/all data.mat'],'DTMS','DTMSerr','DTA','DTAerr','pSR','pSRerr','flocMax','regParams','solTools','-v7.3');
%% L-curve analysis with L2 0th
disp('calculating L-curve with L2 0th...')
M = M_FastBEM;
clear M_FastBEM
MpM = M'*M;
u=vertcat(ux_vec_f,uy_vec_f);
[eyeWeights,~] =getGramMatrix(forceMeshFastBEM);
% original force
force_x_f = force_x(xminf:gridSpacingf:xmax,yminf:gridSpacingf:ymax);
force_y_f = force_y(xminf:gridSpacingf:xmax,yminf:gridSpacingf:ymax);
force_x_f(:,end)=[];
force_x_f(end,:)=[];
force_x_f(1,:)=[];
force_x_f(:,1)=[];
force_y_f(:,end)=[];
force_y_f(end,:)=[];
force_y_f(1,:)=[];
force_y_f(:,1)=[];

force_x_vec_f=reshape(force_x_f,[],1);
force_y_vec_f=reshape(force_y_f,[],1);
force_0=vertcat(force_x_vec_f,force_y_vec_f);

%% examine a logarithmically spaced range of regularization parameters (L2)
lambda=10.^(-9:0.125:-3);
rho=zeros(length(lambda),1);
eta=zeros(length(lambda),1);
fErr=zeros(length(lambda),1);
fCoeff=zeros(size(M,2),length(lambda));
for i=1:length(lambda);
  fCoeff(:,i)=(MpM+lambda(i)*eyeWeights)\(M'*u);
  rho(i)=norm(M*fCoeff(:,i)-u); %residual norm
  eta(i)=norm(eyeWeights*fCoeff(:,i),1); % semi norm
  % force error
  fErr(i)=norm(fCoeff(:,i)-force_0);
end
%% Find the corner of the Tikhonov L-curve
% [reg_corner,ireg_corner,~]=l_curve_corner(rho,eta,lambda);
[reg_corner,ireg_corner,~]=regParamSelecetionLcurve(rho,eta,lambda);
[~,fminIdx]=min(fErr);

save([dataPath '/LcurveL2-0th.mat'],'rho','eta','fErr','reg_cornerL1','ireg_cornerL1','lambda','fCoeff','fminIdx','-v7.3');

%% Plot the sparse deconvolution L-curve.
hLcurve = figure;
set(hLcurve, 'Position', [50 300 800 400])
subplot(1,2,1)
loglog(rho,eta,'k-');
xlabel('Residual Norm ||Gm-d||_{2}');
ylabel('Simi-Norm ||Lm||_{2}');
hold on
% mark and label the corner
H=loglog(rho(ireg_corner),eta(ireg_corner),'r*');
loglog(rho(fminIdx),eta(fminIdx),'ko');
set(H,'markersize',6)
H=text(rho(ireg_corner),1.1*eta(ireg_corner),...
    ['    ',num2str(lambda(ireg_corner),'%5.1e')]);
set(H,'Fontsize',7);
subplot(1,2,2)
loglog(lambda,rho/max(rho(:)),'g')
hold on
loglog(lambda,eta/max(eta(:)),'b')
loglog(lambda,fErr/max(fErr(:)),'r')
% [~,fminIdx]=min(fErr);
loglog(lambda(fminIdx),fErr(fminIdx)/max(fErr(:)),'ko')
loglog(reg_corner,rho(ireg_corner)/max(rho(:)),'g*')
hold on
loglog(reg_corner,eta(ireg_corner)/max(eta(:)),'b*')
loglog(reg_corner,fErr(ireg_corner)/max(fErr(:)),'r*')
% axis([1e-2 100 0.001 1e8])
disp('Printing L-curve...')
% print -deps2 nameSave
print(hLcurve,[imgPath '/Lcurve.eps'],'-depsc')
hgexport(hLcurve,strcat(imgPath,'Lcurve'),hgexport('factorystyle'),'Format','tiff')
hgsave(hLcurve,strcat(imgPath,'Lcurve'),'-v7.3')
%% showing force for L2
[fx,fy,x_out,y_out]=calcForcesFromCoef(forceMeshFastBEM,fCoeff(:,37),grid_mat_f(:,:,1),grid_mat_f(:,:,2),'new');
generateHeatmapFromGridData(x_out,y_out,fx,fy,[dataPath filesep 'L2forcemap at Lcorner'])
[fx,fy,x_out,y_out]=calcForcesFromCoef(forceMeshFastBEM,fCoeff(:,31),grid_mat_f(:,:,1),grid_mat_f(:,:,2),'new');
generateHeatmapFromGridData(x_out,y_out,fx,fy,[dataPath filesep 'L1forcemap at fErr minimum'])

%% examine a logarithmically spaced range of regularization parameters (L1)
MpM = M_FastBEM'*M_FastBEM;
maxIter=4;
tolx=1e-2;
tolr=1e-7;
lambda=10.^(-7:0.125:-1);
rho=zeros(length(lambda),1);
eta=zeros(length(lambda),1);
fErr=zeros(length(lambda),1);
fCoeff=zeros(size(M,2),length(lambda));
for i=1:length(lambda)
  fCoeff(:,i)=iterativeL1Regularization(M,MpM,u,eyeWeights,lambda(i),maxIter,tolx,tolr);
  rho(i)=norm(M*fCoeff(:,i)-u); %residual norm
  eta(i)=norm(eyeWeights*fCoeff(:,i),1); % semi norm
  % force error
  fErr(i)=norm(fCoeff(:,i)-force_0);
end
%% Find the corner of the Tikhonov L-curve
[reg_cornerL1,ireg_cornerL1,~]=regParamSelecetionLcurve(rho,eta,lambda);
[~,fminIdx]=min(fErr);

save([dataPath '/LcurveL1-0th.mat'],'rho','eta','fErr','reg_cornerL1','ireg_cornerL1','lambda','fCoeff','fminIdx','-v7.3');
%% Plot the sparse deconvolution L-curve.
hLcurve = figure;
set(hLcurve, 'Position', [50 300 800 400])
subplot(1,2,1)
loglog(rho,eta,'k-');
xlabel('Residual Norm ||Gm-d||_{2}');
ylabel('Simi-Norm ||Lm||_{2}');
hold on
% mark and label the corner
H=loglog(rho(ireg_cornerL1),eta(ireg_cornerL1),'r*');
loglog(rho(fminIdx),eta(fminIdx),'ko');
set(H,'markersize',6)
H=text(rho(ireg_cornerL1),1.1*eta(ireg_cornerL1),...
    ['    ',num2str(lambda(ireg_cornerL1),'%5.1e')]);
set(H,'Fontsize',7);
subplot(1,2,2)
loglog(lambda,rho/max(rho(:)),'g')
hold on
loglog(lambda,eta/max(eta(:)),'b')
loglog(lambda,fErr/max(fErr(:)),'r')
[~,fminIdx]=min(fErr);
loglog(lambda(fminIdx),fErr(fminIdx)/max(fErr(:)),'ko')
loglog(reg_cornerL1,rho(ireg_cornerL1)/max(rho(:)),'g*')
hold on
loglog(reg_cornerL1,eta(ireg_cornerL1)/max(eta(:)),'b*')
loglog(reg_cornerL1,fErr(ireg_cornerL1)/max(fErr(:)),'r*')
% axis([1e-2 100 0.001 1e8])
disp('Printing L-curve...')
% print -deps2 nameSave
print(hLcurve,[imgPath '/Lcurve L1 0th.eps'],'-depsc')
hgexport(hLcurve,strcat(imgPath,'Lcurve L1 0th'),hgexport('factorystyle'),'Format','tiff')
hgsave(hLcurve,strcat(imgPath,'Lcurve L1 0th'),'-v7.3')

%% showing force
[fx,fy,x_out,y_out]=calcForcesFromCoef(forceMeshFastBEM,fCoeff(:,35),grid_mat_f(:,:,1),grid_mat_f(:,:,2),'new');
generateHeatmapFromGridData(x_out,y_out,fx,fy,[dataPath filesep 'L1forcemap'])
[fx,fy,x_out,y_out]=calcForcesFromCoef(forceMeshFastBEM,fCoeff(:,1),grid_mat_f(:,:,1),grid_mat_f(:,:,2),'new');
generateHeatmapFromGridData(x_out,y_out,fx,fy,[dataPath filesep 'L1forcemap1e-7'])
[fx,fy,x_out,y_out]=calcForcesFromCoef(forceMeshFastBEM,fCoeff(:,33),grid_mat_f(:,:,1),grid_mat_f(:,:,2),'new');
generateHeatmapFromGridData(x_out,y_out,fx,fy,[dataPath filesep 'L1forcemap1e-3'])