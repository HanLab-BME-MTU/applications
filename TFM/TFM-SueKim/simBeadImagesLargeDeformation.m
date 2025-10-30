% This was not chosen to be used for the study. - Sangyoon Mar 2021
% We chose to use simBeadImagesLargerDeformation.m
%% Preparing synthetic bead images
clc
clear all
meshPtsFwdSol=2^9;
xmax=meshPtsFwdSol;
ymax=meshPtsFwdSol;
nPoints = 5000; % was 25000
bead_r = 40; % nm
pixSize = 72; % nm/pix 90x
sigma = 1.73; % was 1.6 before
Aorg = 300+1000*rand(1,nPoints);
[refimg,bead_x, bead_y, ~, Av] = simGaussianBeads(xmax,ymax, sigma, ...
        'npoints', nPoints, 'Border', 'truncated','A',Aorg);
%% may not add noise
%% Displacement Field

E=8000;  %Young's modulus, unit: Pa
forceType = 'groupForce';
gridSpacing = 1;
xmin = gridSpacing;
ymin = gridSpacing;
[x_mat_u, y_mat_u]=meshgrid(xmin:gridSpacing:xmax,ymin:gridSpacing:ymax);

%% posNA - get the coordinates
figure, imshow(ones(500))
hold on
posNA = [];
n = 0;
% Loop, picking up the points.
disp('Left mouse button picks points.')
disp('Right mouse button picks last point.')

but = 1;
while but == 1
    [xi,yi,but] = ginput(1);
    plot(xi,yi,'ro')
    n = n+1;
    text(xi,yi-8,num2str(n));
    posNA(n,:) = [xi yi];
end
hold off
%% posFA - get the coordinates
figure, imshow(ones(500))

hold on
posFA = [];
n = 0
% Loop, picking up the points.
disp('Left mouse button picks points.')
disp('Right mouse button picks last point.')

but = 1;
while but == 1
    [xi,yi,but] = ginput(1);
    plot(xi,yi,'ro')
    n = n+1;
    text(xi,yi-8,num2str(n));
    posFA(n,:) = [xi yi];
end
hold off

% %% posBigFA - get the coordinates
% figure, imshow(ones(500))
% hold on
% posBigFA = [];
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
%     posBigFA(n,:) = [xi yi];
% end
% hold off

%% display selected positions
disp([posFA])
disp([posNA])
% disp([posBigFA])

posFA = [ 103.0000  352.0000
  174.0000  277.0000
  232.0000  211.0000
  312.0000  173.0000
  391.0000  146.0000
  368.0000  217.0000
  313.0000  262.0000
  278.0000  311.0000
  208.0000  337.0000
  163.0000  366.0000
  181.0000  222.0000
  267.0000  154.0000 ];

posNA = [ 116.0000  232.0000
  198.0000  152.0000
  314.0000  121.0000 ];

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
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(1,1),posFA(1,2),800,1600,500/108,2000/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(2,1),posFA(2,2),600,2000,500/500,2100/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(3,1),posFA(3,2),2500,10000,2000/108,5200/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(4,1),posFA(4,2),440,2800,1000/108,2500/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(5,1),posFA(5,2),160,1840,500/108,2000/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(6,1),posFA(6,2),80,1900,500/108,2300/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(7,1),posFA(7,2),0,3000,800/108,2500/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(8,1),posFA(8,2),100,1900,1000/108,3200/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(9,1),posFA(9,2),200,1840,500/108,2200/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(10,1),posFA(10,2),300,2500,400/108,2000/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(11,1),posFA(11,2),450,2800,1000/108,2500/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(12,1),posFA(12,2),500,3000,1000/108,3000/108,forceType);
      
force_y = assumedForceAniso2D(2,x_mat_u,y_mat_u,posNA(1,1),posNA(1,2),300,620,400/72,500/72,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posNA(2,1),posNA(2,2),200,650,400/72,500/72,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posNA(3,1),posNA(3,2),150,700,400/72,500/72,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posFA(1,1),posFA(1,2),800,1600,500/108,2000/108,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posFA(2,1),posFA(2,2),600,2000,500/500,2100/108,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posFA(3,1),posFA(3,2),2500,10000,2000/108,5200/108,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posFA(4,1),posFA(4,2),440,2800,1000/108,2500/108,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posFA(5,1),posFA(5,2),160,1840,500/108,2000/108,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posFA(6,1),posFA(6,2),80,1900,500/108,2300/108,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posFA(7,1),posFA(7,2),0,3000,800/108,2500/108,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posFA(8,1),posFA(8,2),100,1900,1000/108,3200/108,forceType)+...
          assumedForceAniso2D(2,x_mat_u,y_mat_u,posFA(9,1),posFA(9,2),200,1840,500/108,2200/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(10,1),posFA(10,2),300,2500,400/108,2000/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(11,1),posFA(11,2),450,2800,1000/108,2500/108,forceType)+...
          assumedForceAniso2D(1,x_mat_u,y_mat_u,posFA(12,1),posFA(12,2),500,3000,1000/108,3000/108,forceType);
 
%% force noise
maxFnoise = 100;
fnorm_org = (force_x.^2 + force_y.^2).^0.5; %this should be fine mesh
foreground = fnorm_org>120;
background = ~foreground;
force_x = (maxFnoise*rand(ymax,xmax)).*background + force_x;
force_y = (maxFnoise*rand(ymax,xmax)).*background + force_y;
fnorm_org = (force_x.^2 + force_y.^2).^0.5; %this should be fine mesh
figure, imshow(fnorm_org,[0 500]), colormap jet

%% displacement field
[ux, uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,...
    force_x,force_y,'fft',[],meshPtsFwdSol,50000,0.5,false); %,'conv',[],meshPtsFwdSol);
figure, imshow(uy,[]), colormap jet

%% finding displacement at bead location
nPoints = length(bead_x);
bead_ux = zeros(size(bead_x));
bead_uy = zeros(size(bead_y));
for k=1:nPoints
    % identify closest integer location in mat_u
    [~,indcol_closest_x] = min(abs(x_mat_u(1,:)-bead_x(k)),[],2);
    [~,indrow_closest_y] = min(abs(y_mat_u(:,1)-bead_y(k)),[],1);
    % get the neighborhood with the x,y location
    row_bottom = max(1,indrow_closest_y-2);
    row_top = min(size(x_mat_u,2),indrow_closest_y+2);
    col_bottom = max(1,indcol_closest_x-2);
    col_top = min(size(y_mat_u,1),indcol_closest_x+2);
    loc_xmat = x_mat_u(row_bottom:row_top,col_bottom:col_top);
    loc_ymat = y_mat_u(row_bottom:row_top,col_bottom:col_top);
    loc_ux = ux(row_bottom:row_top,col_bottom:col_top);
    loc_uy = uy(row_bottom:row_top,col_bottom:col_top);
    % interpolate for subpixel location
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

%% datapath
dataPath='/storage/network/TFM_Development/TFM2D/PIVimprovement/simulation';
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
generateHeatmapFromGridData(x_mat_u,y_mat_u,force_x,force_y,[dataPath '/Original forcefield'],0,0,3100,false,430,430);
%% display original forcefield with 200 max
generateHeatmapFromGridData(x_mat_u,y_mat_u,force_x,force_y,[dataPath '/Original forcefield 200'],0,0,200,false,430,430);
%% display original forcefield with 500 max
generateHeatmapFromGridData(x_mat_u,y_mat_u,force_x,force_y,[dataPath '/Original forcefield 500'],0,0,500,false,430,430);

%% display original displacement field
    generateHeatmapFromGridData(x_mat_u,y_mat_u,ux,uy,[dataPath '/Original displacementField'],0,0,16,false,430,430);
    %% display measured displacement field
    % load displacement field
    displPath = [dataPath filesep 'TFMPackage/displacementField'];
    displFile =[dataPath filesep 'TFMPackage/displacementField/displField.mat'];
    load(displFile)
    generateHeatmapFromField(displField,displPath,0,16,'jet',430,430,false);
    %% display corrected displacement field
    % load displacement field
    displPath = [dataPath filesep 'TFMPackage/correctedDisplacementField'];
    displFile =[dataPath filesep 'TFMPackage/correctedDisplacementField/displField.mat'];
    load(displFile)
    generateHeatmapFromField(displField,displPath,0,16,[],430,430,false);
%% view refimg and beadimg
figure, imshow(refimg,[])
figure, imshow(beadimg,[])
