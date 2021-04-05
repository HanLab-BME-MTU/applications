%% set up
clc
clear
%% Preparing synthetic bead images
meshPtsFwdSol=2^9;
xmax=meshPtsFwdSol;
ymax=meshPtsFwdSol;
nPoints = 8000; % was 10000
bead_r = 40; % nm
pixSize = 72; % nm/pix 90x
sigma = 1.73; % was 1.6 before
Aorg = 300+1000*rand(1,nPoints);
[refimg,bead_x, bead_y, ~, Av] = simGaussianBeads(xmax,ymax, sigma, ...
        'npoints', nPoints, 'Border', 'truncated','A',Aorg);
    
%% Displacement Field
E=8000;  %Young's modulus, unit: Pa %changed from 8000 to 5000
forceType = 'groupForce';
gridSpacing = 1;
xmin = gridSpacing;
ymin = gridSpacing;
[x_mat_u, y_mat_u]=meshgrid(xmin:gridSpacing:xmax,ymin:gridSpacing:ymax);

%% posNA - get the coordinates
% figure, imshow(ones(500))
% hold on
% posNA = [];
% n = 0;
% % Loop, picking up the points.
% disp('Left mouse button picks points.')
% disp('Right mouse button picks last point.')
% 
% but = 1;
% while but == 1
%     [xi,yi,but] = ginput(1);
%     plot(xi,yi,'ro')
%     n = n+1;
%     text(xi,yi-8,num2str(n));
%     posNA(n,:) = [xi yi];
% end
% hold off
% %% posFA- get the coordinates
% figure, imshow(ones(500))
% 
% hold on
% posFA = [];
% n = 0;
% % Loop, picking up the points.
% disp('Left mouse button picks points.')
% disp('Right mouse button picks last point.')
% 
% but = 1;
% while but == 1
%     [xi,yi,but] = ginput(1);
%     plot(xi,yi,'ro')
%     n = n+1;
%     text(xi,yi-8,num2str(n));
%     posFA(n,:) = [xi yi];
% end
% hold off

%% display selected positions
% disp([posFA])
% disp([posNA])
% copy and paste coordinates of postFA and posNA from command window

posFA = [347.0000  135.0000
  239.0000  109.0000
  161.0000  132.0000
  158.0000  220.0000
  176.0000  332.0000
  265.0000  209.0000
  233.0000  377.0000
  440.0000  206.0000];

posNA = [74.0000  318.0000
  112.0000  149.0000
  234.0000   68.0000
  396.0000  102.0000
  125.0000   48.0000];
%%  ------------------------ original force
% for ii = numel(posNA)
% force_x = assumedForceAniso2D(1,x_mat_u,y_mat_u,posNA(ii,1),posNA(ii,2),
% force = assumedForceAniso2D(j,x,y,xshift,yshift,wx,wy,d1,d2,forceType) takes grid x and y and make anisotropic gaussian
%   distributed force field of which sources are xshift and yshift.
%   input     :   x           grid of x coordinates
%                 y           grid of y coordinates
%                 xshift      x value of point source of force
%                 yshift      y value of point source of force
%                 wx          x component of force orientation
%                 wy          y component of force orientation
%                 d1          diameter of adhesion in the direction of minor
%                 axis
%                 d2          diameter of adhesion in the direction of major
%                 axis
%                 forceType   type of force
%                 ('pointForce','groupForce', or 'smoothForce')

force_x = assumedForceAniso2D(1,x_mat_u,y_mat_u,posNA(1,1),posNA(1,2),300,620,400/72,500/72,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posNA(2,1),posNA(2,2),200,650,400/72,500/72,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posNA(3,1),posNA(3,2),150,700,400/72,500/72,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posNA(4,1),posNA(4,2),125,420,400/72,500/72,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posNA(5,1),posNA(5,2),175,520,400/72,500/72,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(1,1),posFA(1,2),500,1850,700/108,700/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(2,1),posFA(2,2),500,1500,1000/108,1000/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(3,1),posFA(3,2),800,1800,2000/108,2000/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(4,1),posFA(4,2),500,1500,1000/108,1000/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(5,1),posFA(5,2),500,1840,700/108,700/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(6,1),posFA(6,2),8000,10000,1500/108,10000/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(7,1),posFA(7,2),500,2000,1000/108,2500/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(8,1),posFA(8,2),500,2000,1000/108,2500/108,forceType);
      
force_y = assumedForceAniso2D(2,x_mat_u,y_mat_u,posNA(1,1),posNA(1,2),300,620,400/72,500/72,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posNA(2,1),posNA(2,2),200,650,400/72,500/72,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posNA(3,1),posNA(3,2),150,700,400/72,500/72,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posNA(4,1),posNA(4,2),125,420,400/72,500/72,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posNA(5,1),posNA(5,2),175,520,400/72,500/72,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posFA(1,1),posFA(1,2),500,1850,700/108,700/108,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posFA(2,1),posFA(2,2),500,1500,1000/108,1000/108,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posFA(3,1),posFA(3,2),800,1800,2000/108,2000/108,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posFA(4,1),posFA(4,2),500,1500,1000/108,1000/108,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posFA(5,1),posFA(5,2),500,1840,700/108,700/108,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posFA(6,1),posFA(6,2),8000,10000,1500/108,10000/108,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posFA(7,1),posFA(7,2),500,2000,1000/108,2500/108,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posFA(8,1),posFA(8,2),500,2000,1000/108,2500/108,forceType);
      
%% force noise
% figure 1
maxFnoise = 100;
fnorm_org = (force_x.^2 + force_y.^2).^0.5; %this should be fine mesh
foreground = fnorm_org>120;
background = ~foreground;
force_x = (maxFnoise*rand(ymax,xmax)).*background + force_x;
force_y = (maxFnoise*rand(ymax,xmax)).*background + force_y;
fnorm_org = (force_x.^2 + force_y.^2).^0.5; %this should be fine mesh
figure, imshow(fnorm_org,[0 12000]), colormap jet

%% displacement field
% figure 2
[ux, uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,...
    force_x,force_y,'fft',[],meshPtsFwdSol,50000,0.5,false); %,'conv',[],meshPtsFwdSol);
figure, imshow(uy,[]), colormap jet
%%
% % figure (3); surf(x_mat_u,y_mat_u,ux)
% % [x_mat_u,y_mat_u,uy] = peaks(512);
% % uy(:,:,1) = zeros(512); % red
% % uy(:,:,2) = ones(512).*linspace(0,256,512); % green
% % uy(:,:,3) = ones(512).*linspace(0,1,512); % blue
% %[X,Y] = meshgrid(x_mat_u,y_mat_u);
% Z = uy;Crop The rectangular crop area must not be outside the image
% 
% figure (3); surf(x_mat_u,y_mat_u,uy)
% grid on
% iy = y_mat_u(x_mat_u == 1);
% iz = Z(y_mat_u == 1);
% 
% figure(4); 
% plot(iy,iz)
% grid on
% 
% figure(5);
% ix = x_mat_u(y_mat_u == 256);
% iz = uy(y_mat_u == 256);
% plot(ix,iz)
% grid on

% % xslice = 300;
% % yslice = [];
% % zslice = [];
% 
% % grid on
% % iy = y_mat_u(x_mat_u == 1);
% % iz = uy(x_mat_u == 1);
%% finding displacement at bead location
nPoints = length(bead_x);
bead_ux = zeros(size(bead_x));
bead_uy = zeros(size(bead_y));
for k=1:nPoints
    [~,indcol_closest_x] = min(abs(x_mat_u(1,:)-bead_x(k)),[],2);
    [~,indrow_closest_y] = min(abs(y_mat_u(:,1)-bead_y(k)),[],1);
    row_bottom = max(1,indrow_closest_y-2);
    row_top = min(size(x_mat_u,2),indrow_closest_y+2);
    col_bottom = max(1,indcol_closest_x-2);
    col_top = min(size(y_mat_u,1),indcol_closest_x+2);
    loc_xmat = x_mat_u(row_bottom:row_top,col_bottom:col_top);
    loc_ymat = y_mat_u(row_bottom:row_top,col_bottom:col_top);
    loc_ux = ux(row_bottom:row_top,col_bottom:col_top);
    loc_uy = uy(row_bottom:row_top,col_bottom:col_top);
    bead_ux(k) = interp2(loc_xmat,loc_ymat,loc_ux,bead_x(k),bead_y(k));
    if isnan(bead_ux(k))
        bead_ux(k) = ux(indrow_closest_y,indcol_closest_x);
    end
    bead_uy(k) = interp2(loc_xmat,loc_ymat,loc_uy,bead_x(k),bead_y(k));
    if isnan(bead_uy(k))
        bead_uy(k) = uy(indrow_closest_y,indcol_closest_x);
    end
end

% pixelSize = 0.108; % assuming 60x objective um/pixel
beadimg = simGaussianBeads(xmax,ymax, sigma,'x',bead_x+bead_ux,'y',bead_y+bead_uy,'A',Av, 'Border', 'truncated');

%% saving
% dataPath='/storage/network/TFM_Development/TFM2D/PIVimprovement/simulation';
dataPath=pwd;
%PIVImprovement - simulation folder
imgPath=[dataPath filesep 'Beads'];
refPath=[dataPath filesep 'Reference'];
if ~exist(refPath,'dir')
    mkdir(imgPath);
    mkdir(refPath);
end
imwrite(uint16(refimg*2^16/max(max(refimg))),[refPath filesep 'img1ref.tif'],'tif')
imwrite(uint16(beadimg*2^16/max(max(beadimg))),[imgPath filesep 'img2bead.tif'],'tif')
%% display original forcefield
generateHeatmapFromGridData(x_mat_u,y_mat_u,force_x,force_y,[dataPath '/Original forcefield'],0,0,6000,false,430,430);
%% display measured displacement field
% load displacement field
    displPath = [dataPath filesep 'TFMPackage/displacementField'];
    displFile = [displPath filesep 'displField.mat'];  %/displField.mat'
    load(displFile)
    generateHeatmapFromField(displField,displPath,0,40,'jet',430,430,false);
    save('displField.mat')
    %% display corrected displacement field
    % load displacement field
    displPath = [dataPath filesep 'TFMPackage/correctedDisplacementField'];
    displFile = [displPath filesep 'displField.mat'];  %/displField.mat'
    load(displFile)
    generateHeatmapFromField(displField,displPath,0,16,[],430,430,false);
%% view refimg and beadimg
figure, imshow(refimg,[])
figure, imshow(beadimg,[])
