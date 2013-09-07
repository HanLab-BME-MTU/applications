%% Preparing synthetic bead images

% reference image (500x500)
% reference image (300x200)
xmax=500;
ymax=500;
nPoints = 25000;
bead_r = 50; % nm
pixSize = 72; % nm/pix 90x
sigma = 1.6;
[refimg,bead_x, bead_y, ~, Av] = simGaussianBeads(xmax,ymax, sigma, ...
        'npoints', nPoints, 'Border', 'truncated','A',0.3+rand(1,nPoints));
%% Noise addition (5%)
refimg = refimg+0.05*rand(ymax,xmax)*max(refimg(:));

% bead images
%% displacement field
E=8000;  %Young's modulus, unit: Pa
meshPtsFwdSol=2^10;
forceType = 'groupForce';

gridSpacing = 1;
xmin = gridSpacing;
ymin = gridSpacing;

[x_mat_u, y_mat_u]=meshgrid(xmin:gridSpacing:xmax,ymin:gridSpacing:ymax);

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

posNA = [139.0000  267.0000
  156.0000  232.0000
  184.0000  212.0000
  217.0000  200.0000
  246.0000  195.0000
  272.0000  199.0000
  297.0000  211.0000
  323.0000  231.0000
  346.0000  266.0000];
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
posFA = [   42.0000  323.0000
   54.0000  319.0000
   96.0000  298.0000
  139.0000  288.0000
  180.0000  280.0000
  225.0000  276.0000
  263.0000  273.0000
  301.0000  275.0000
  331.0000  281.0000
  347.0000  279.0000
  381.0000  290.0000
  417.0000  303.0000
  455.0000  317.0000
   69.0000  359.0000
   91.0000  349.0000
  115.0000  337.0000
  129.0000  338.0000
  168.0000  327.0000
  186.0000  321.0000
  270.0000  318.0000
  288.0000  320.0000
  312.0000  323.0000
  361.0000  330.0000
  383.0000  337.0000
  423.0000  352.0000
  432.0000  363.0000];

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
%% Noise addition (5%)
beadimg = beadimg+0.05*rand(ymax,xmax)*max(beadimg(:));
%% saving
dataPath='/hms/scratch1/sh268/multiForceTesting2';
imgPath=[dataPath filesep 'Beads'];
refPath=[dataPath filesep 'Reference'];
if ~exist(refPath,'dir')
    mkdir(imgPath);
    mkdir(refPath);
end

imwrite(uint16(refimg*2^16/max(max(refimg))),[refPath filesep 'img1ref.tif'],'tif')
imwrite(uint16(beadimg*2^16/max(max(beadimg))),[imgPath filesep 'img2bead.tif'],'tif')

%%  ------------------------ original force

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
%-------------------------------------------

% run the TFM Package and obtain displacement field and force field

%--------------------------------------------
%% display original forcefield
generateHeatmapFromGridData(x_mat_u,y_mat_u,force_x,force_y,[dataPath '/Original forcefield'],3100,460,460)

%% display original displacement field
    generateHeatmapFromGridData(x_mat_u,y_mat_u,ux,uy,dataPath)
    %% display measured displacement field
    % load displacement field
    displPath = [dataPath filesep 'TFMPackage/displacementField'];
    displFile =[dataPath filesep 'TFMPackage/displacementField/displField.mat'];
    load(displFile)
    generateHeatmapFromField(displField,displPath,9.6);
    %% display corrected displacement field
    % load displacement field
    displPath = [dataPath filesep 'TFMPackage/correctedDisplacementField'];
    displFile =[dataPath filesep 'TFMPackage/correctedDisplacementField/displField.mat'];
    load(displFile)
    generateHeatmapFromField(displField,displPath,9.6);

%     %% display measured displacement field
%     % load displacement field
%     displPath = '/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/multiForceTesting/TFMPackage/displacementField cL=11';
%     displFile ='/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/multiForceTesting/TFMPackage/displacementField cL=11/displField.mat';
%     load(displFile)
%     generateHeatmapFromField(displField,displPath,9.6)
% 
%% display measured corrected displacement field
% load displacement field
load('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/multiForceTesting/TFMPackage/correctedDisplacementField/displField.mat')
[reg_grid,~,~,spacing]=createRegGridFromDisplField(displField,2); %2=2 times fine interpolation
[grid_mat,iu_mat,~,~] = interp_vec2grid(displField.pos, displField.vec,[],reg_grid);
pos = [reshape(grid_mat(:,:,1),[],1) reshape(grid_mat(:,:,2),[],1)]; %dense
disp_vec = [reshape(iu_mat(:,:,1),[],1) reshape(iu_mat(:,:,2),[],1)]; 

umnorm = (iu_mat(:,:,1).^2 + iu_mat(:,:,2).^2).^0.5;
umMap = griddata(grid_mat(:,:,1),grid_mat(:,:,2),umnorm,XI,YI,'cubic');
ummin = min(min(umnorm));
ummax = max(max(umnorm));

h3=figure('color','w');
subplot('Position',[0 0 0.8 1])
imshow(umMap,[ummin ummax]), colormap jet;
%quiver plot
hold on
hqm = quiver(displField.pos(:,1),displField.pos(:,2), displField.vec(:,1)./dispScale,displField.vec(:,2)./dispScale,0,'Color',[75/255 0/255 130/255]);

subplot('Position',[0.8 0.1 0.1 0.8])
axis tight
caxis([ummin ummax]), axis off
hc = colorbar('West');
set(hc,'Fontsize',18)

% saving
% Set up the output file path
outputFilePath = [dataPath filesep 'displacementField corrected'];
tifPath = [outputFilePath filesep 'tifs'];
figPath = [outputFilePath filesep 'figs'];
epsPath = [outputFilePath filesep 'eps'];
if ~exist(tifPath,'dir') || ~exist(epsPath,'dir')
    mkdir(tifPath);
    mkdir(figPath);
    mkdir(epsPath);
end

I = getframe(h3);
imwrite(I.cdata, strcat(tifPath,'/displFieldMagTif','.tif'));
hgsave(h3,strcat(figPath,'/displFieldMagFig','-v7.3'))
print(h3,strcat(epsPath,'/displFieldMagEps.eps'),'-depsc2')
    %% display calculated force field
% load displacement field
forcePath = '/hms/scratch1/sh268/multiForceTesting2/TFMPackage/forceField_L2_2nd Lcurve lambda';
forceFile ='/hms/scratch1/sh268/multiForceTesting2/TFMPackage/forceField/forceField.mat';
load(forceFile)
generateHeatmapFromField(forceField,forcePath);

%% L-curve for L1 0th
MpM=M'*M;
maxIter = 7;
tolx = 0.2;
tolr = 1e-7;
disp('L-curve ...')
%         [~,L] = calculateLfromLcurveSparse(M,MpM,u,eyeWeights,maxIter,tolx,tolr,solMethodBEM);
alphas=10.^(-7:.125:-2);
rho=zeros(length(alphas),1);
eta=zeros(length(alphas),1);
msparse=zeros(size(M,2),length(alphas));
fErr=zeros(length(alphas),1);
[eyeWeights,~] =getGramMatrix(forceMesh);

% original force at force base nodes
xminf = forceMesh.basis(1).node(1);
yminf = forceMesh.basis(1).node(2);
gridSpacingf = forceMesh.basis(2).node(2)-forceMesh.basis(1).node(2);
xmaxf = forceMesh.basis(end).node(1);
ymaxf = forceMesh.basis(end).node(2);

force_x_f = force_x(xminf:gridSpacingf:xmaxf,yminf:gridSpacingf:ymaxf);
force_y_f = force_y(xminf:gridSpacingf:xmaxf,yminf:gridSpacingf:ymaxf);

force_x_vec_f=reshape(force_x_f,[],1);
force_y_vec_f=reshape(force_y_f,[],1);
force_0=vertcat(force_x_vec_f,force_y_vec_f);

%% L-curve L1 0th with conventianal force error criterion
for i=1:length(alphas);
  msparse(:,i)=iterativeL1Regularization(M,MpM,u,eyeWeights,alphas(i),maxIter,tolx,tolr);
  rho(i)=norm(M*msparse(:,i)-u);
  eta(i)=norm(eyeWeights*msparse(:,i),1);
  % force error
  fErr(i)=norm(msparse(:,i)-force_0);
  disp([num2str(i) ' out of ' num2str(length(alphas))]);
end
%% new force error criterion
% points at the BG
posBG = [   42.0000  423.0000
   54.0000  419.0000
   96.0000  398.0000
  139.0000  438.0000
  180.0000  430.0000
  225.0000  476.0000
  263.0000  433.0000
  301.0000  435.0000
  331.0000  481.0000
  347.0000  479.0000
  381.0000  430.0000
  417.0000  403.0000
  455.0000  417.0000
   69.0000  459.0000
   91.0000  449.0000
  115.0000  437.0000
  129.0000  438.0000
  168.0000  427.0000
  186.0000  421.0000
  270.0000  418.0000
  288.0000  420.0000
  312.0000  423.0000
  361.0000  430.0000
  383.0000  437.0000
  423.0000  452.0000
  432.0000  463.0000
  139+5.0000  267.0000
  156.0000  232+5.0000
  184+5.0000  212.0000
  217.0000  200+5.0000
  246+5.0000  195.0000
  272.0000  199+5.0000
  297.0000  211+5.0000
  323+5.0000  231.0000
  346.0000  266+5.0000
  42+10.0000  323.0000
   54.0000  319+10.0000
   96+10.0000  298.0000
  139.0000  288+10.0000
  180.0000  280+10.0000
  225+10.0000  276.0000
  263.0000  273+10.0000
  301+10.0000  275.0000
  331.0000  281+10.0000
  347+10.0000  279.0000
  381.0000  290+10.0000
  417+10.0000  303.0000
  455.0000  317+10.0000
   69+10.0000  359.0000
   91.0000  349+10.0000
  115+10.0000  337.0000
  129.0000  338+10.0000
  168+10.0000  327.0000
  186.0000  321+10.0000
  270+10.0000  318.0000
  288.0000  320+10.0000
  312+10.0000  323.0000
  361.0000  330+10.0000
  383+10.0000  337.0000
  423.0000  352+10.0000
  432+10.0000  363+10.0000];
locMaxI_FA = posFA(:,2:-1:1);
locMaxI_NA = posNA(:,2:-1:1);
locMaxI_BG = posBG(:,2:-1:1);
locMaxI_tot = [locMaxI_FA; locMaxI_NA; locMaxI_BG];
%% merge them
for i=1:length(alphas)
    [fx,fy,~,~]=calcForcesFromCoef(forceMesh,msparse(:,i),x_mat_u,y_mat_u,'new');
    fErr(i)=sum(sqrt((diag(force_x(locMaxI_tot(:,1),locMaxI_tot(:,2)))-...
        diag(fx(locMaxI_tot(:,1),locMaxI_tot(:,2)))).^2+...
        (diag(force_y(locMaxI_tot(:,1),locMaxI_tot(:,2)))-...
        diag(fy(locMaxI_tot(:,1),locMaxI_tot(:,2)))).^2));
    disp([num2str(i) '     ' num2str(alphas(i)) '      ' num2str(fErr(i))]);
end

%% Find the corner of the Tikhonov L-curve for L1 0th
% [reg_corner,ireg_corner,~]=l_curve_corner(rho,eta,alphas);
[reg_corner,ireg_corner,~]=regParamSelecetionLcurve(rho,eta,alphas);
% [reg_corner,ireg_corner,~]=regParamSelecetionLcurve(rho,eta,alphas);
[~,fminIdx]=min(fErr);
% [~,fminIdx]=min(fErr3);
% save([dataPath '/LcurveL1-0th.mat'],'rho','eta','rho2','eta2','rho3','eta3','fErr3','reg_corner','ireg_corner','alphas','msparse','fminIdx','-v7.3');
save([dataPath '/LcurveL1-0th.mat'],'forceMesh','rho','eta','fErr','reg_corner','ireg_corner','alphas','msparse','fminIdx','-v7.3');
%% showing force for L1 0th
load([dataPath '/LcurveL1-0th.mat']);

[fx,fy,x_out,y_out]=calcForcesFromCoef(forceMesh,msparse(:,ireg_corner),x_mat_u,y_mat_u,'new');
% msparse_Lcorner 
% forcePath = '/hms/scratch1/sh268/multiForceTesting2/TFMPackage/L1-0th forcemap at Lcorner';
% forceField.pos=[reshape(x_out,[],1) reshape(y_out,[],1)];
% forceField.vec=[reshape(fx,[],1) reshape(fy,[],1)];
% 
% generateHeatmapFromField(forceField,forcePath,3100);

generateHeatmapFromGridData(x_out,y_out,fx,fy,[dataPath filesep 'L1-0th forcemap at Lcorner'],3100,false,460,460)
[fx,fy,x_out,y_out]=calcForcesFromCoef(forceMesh,msparse(:,ceil((fminIdx+ireg_corner)/2)),x_mat_u,y_mat_u,'new');
generateHeatmapFromGridData(x_out,y_out,fx,fy,[dataPath filesep 'L1-0th forcemap at fErr minimum'],3100,false,460,460);

 %% L-curve for L1 2nd
MpM=M'*M;
maxIter = 10;
tolx = 2e-2;
tolr = 1e-7;
disp('L-curve ...')
%         [~,L] = calculateLfromLcurveSparse(M,MpM,u,eyeWeights,maxIter,tolx,tolr,solMethodBEM);
alphas=10.^(-7:.125:-2);
rho=zeros(length(alphas),1);
eta=zeros(length(alphas),1);
msparse=zeros(size(M,2),length(alphas));
fErr=zeros(length(alphas),1);
% Laplacian
nBasis = forceMesh.numBasis;
basisx=zeros(nBasis,1);
basisy=zeros(nBasis,1);
for k=1:nBasis
    basisx(k) = forceMesh.basis(k).node(1);
    basisy(k) = forceMesh.basis(k).node(2);
end
nBasisx = size(unique(basisx),1);
nBasisy = size(unique(basisy),1);

%% Laplacian
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
%% Lcurve for L1 2nd
for i=1:length(alphas);
  msparse(:,i)=iterativeL1Regularization(M,MpM,u,-Lap,alphas(i),maxIter,tolx,tolr);
  rho(i)=norm(M*msparse(:,i)-u);
  eta(i)=norm(-Lap*msparse(:,i),1);
  % force error
  fErr(i)=norm(msparse(:,i)-force_0);
  disp([num2str(i) 'out of ' num2str(length(alphas))]);
end

%% Find the corner of the Tikhonov L-curve
% [reg_corner,ireg_corner,~]=l_curve_corner(rho,eta,alphas);
[reg_corner,ireg_corner,~]=regParamSelecetionLcurve(rho,eta,alphas);
[~,fminIdx]=min(fErr);
[xgrid,ygrid]=meshgrid(xminf:gridSpacingf:xmaxf,yminf:gridSpacingf:ymaxf);

save([dataPath '/LcurveL1-2nd.mat'],'forceMesh','rho','eta','fErr','reg_corner','ireg_corner','alphas','msparse','fminIdx','xgrid','ygrid','-v7.3');
%% showing force for L1 2nd
load([dataPath '/LcurveL1-2nd.mat'])
[fx,fy,x_out,y_out]=calcForcesFromCoef(forceMesh,msparse(:,ireg_corner),xgrid,ygrid,'new');
generateHeatmapFromGridData(x_out,y_out,fx,fy,[dataPath filesep 'L1-2nd forcemap at Lcorner'],3100)
[fx,fy,x_out,y_out]=calcForcesFromCoef(forceMesh,msparse(:,fminIdx),xgrid,ygrid,'new');
generateHeatmapFromGridData(x_out,y_out,fx,fy,[dataPath filesep 'L1-2nd forcemap at fErr minimum'],3100)
%% L-curve analysis for L2 0th
disp('calculating L-curve with L2 0th...')
load('/hms/scratch1/sh268/multiForceTesting2/TFMPackage/forceField/BEMParams.mat');
MpM = M'*M;
[eyeWeights,~] =getGramMatrix(forceMesh);
% original force at force base nodes
xminf = forceMesh.basis(1).node(1);
yminf = forceMesh.basis(1).node(2);
gridSpacingf = forceMesh.basis(2).node(2)-forceMesh.basis(1).node(2);
xmaxf = forceMesh.basis(end).node(1);
ymaxf = forceMesh.basis(end).node(2);

force_x_f = force_x(xminf:gridSpacingf:xmaxf,yminf:gridSpacingf:ymaxf);
force_y_f = force_y(xminf:gridSpacingf:xmaxf,yminf:gridSpacingf:ymaxf);

force_x_vec_f=reshape(force_x_f,[],1);
force_y_vec_f=reshape(force_y_f,[],1);
force_0=vertcat(force_x_vec_f,force_y_vec_f);
%% L2 0th L curve calculation
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
  disp([num2str(i) ' out of ' num2str(length(lambda))]);
end

%% Find the corner of the Tikhonov L-curve
[reg_corner,ireg_corner,kappa]=regParamSelecetionLcurve(rho,eta,lambda);
[~,fminIdx]=min(fErr);

save([dataPath '/LcurveL2-0th.mat'],'forceMesh','rho','eta','fErr','reg_corner','ireg_corner','lambda','fCoeff','fminIdx','-v7.3');

%% showing force for L2
load([dataPath '/LcurveL2-0th.mat'],'rho','eta','fErr','reg_corner','ireg_corner','lambda','fCoeff','fminIdx');
[fx,fy,x_out,y_out]=calcForcesFromCoef(forceMesh,fCoeff(:,ireg_corner),xgrid,ygrid,'new');
% generateHeatmapFromField(forceField,forcePath,3100,'jet');
generateHeatmapFromGridData(x_out,y_out,fx,fy,[dataPath filesep 'L2-0th forcemap at Lcorner'],3100,460,460)
[fx,fy,x_out,y_out]=calcForcesFromCoef(forceMesh,fCoeff(:,fminIdx),xgrid,ygrid,'new');
generateHeatmapFromGridData(x_out,y_out,fx,fy,[dataPath filesep 'L2-0th forcemap at fErr minimum'],3100,460,460)
%% Lcurve for L2 2nd 
lambda=10.^(-10:0.125:-3);
rho=zeros(length(lambda),1);
eta=zeros(length(lambda),1);
fErr=zeros(length(lambda),1);
fCoeff=zeros(size(M,2),length(lambda));
for i=1:length(lambda);
  fCoeff(:,i)=(MpM+lambda(i)*(-Lap))\(M'*u);
  rho(i)=norm(M*fCoeff(:,i)-u); %residual norm
  eta(i)=norm((-Lap)*fCoeff(:,i),1); % semi norm
  % force error
  fErr(i)=norm(fCoeff(:,i)-force_0);
  disp([num2str(i) ' out of ' num2str(length(lambda))]);
end
[reg_corner,ireg_corner,kappa]=regParamSelecetionLcurve(rho,eta,lambda);
[~,fminIdx]=min(fErr);

save([dataPath '/LcurveL2-2nd.mat'],'forceMesh','rho','eta','fErr','reg_corner','ireg_corner','lambda','fCoeff','fminIdx','-v7.3');

%% showing force for L2 2nd
% load([dataPath '/LcurveL2-0th.mat'],'rho','eta','fErr','reg_corner','ireg_corner','lambda','fCoeff','fminIdx');
[fx,fy,x_out,y_out]=calcForcesFromCoef(forceMesh,fCoeff(:,ireg_corner),xgrid,ygrid,'new');
% generateHeatmapFromField(forceField,forcePath,3100,'jet');
generateHeatmapFromGridData(x_out,y_out,fx,fy,[dataPath filesep 'L2-2nd forcemap at Lcorner'],3100)
[fx,fy,x_out,y_out]=calcForcesFromCoef(forceMesh,fCoeff(:,fminIdx),xgrid,ygrid,'new');
generateHeatmapFromGridData(x_out,y_out,fx,fy,[dataPath filesep 'L2-2nd forcemap at fErr minimum'],3100)


%% Quantify forces at FAs and noise levle at background
nMethods =3;
locMaxI_FA = posFA(:,2:-1:1);
locMaxI_NA = posNA(:,2:-1:1);

nPeaks_FA = size(locMaxI_FA,1);
flocMaxFA = zeros(nPeaks_FA, nMethods);
flocMaxRatioFA = zeros(nPeaks_FA, nMethods);
pSR_FA = zeros(nMethods,1); %peak stress ratio
pSRerr_FA = zeros(nMethods,1);

nPeaks_NA = size(locMaxI_NA,1);
flocMaxNA = zeros(nPeaks_NA, nMethods);
flocMaxRatioNA = zeros(nPeaks_NA, nMethods);
pSR_NA = zeros(nMethods,1); %peak stress ratio
pSRerr_NA = zeros(nMethods,1);

fnorm_org = (force_x.^2 + force_y.^2).^0.5; %this should be fine mesh
flocMaxOrgFA = diag(fnorm_org(locMaxI_FA(:,1),locMaxI_FA(:,2)));
flocMaxOrgNA = diag(fnorm_org(locMaxI_NA(:,1),locMaxI_NA(:,2)));
% r = 1.1; %
forground = fnorm_org>0;
backgroundMaskfine = ~bwmorph(forground,'dilate',5);
[xgrid_fine,ygrid_fine] = meshgrid(xgrid(1,1):xgrid(end,end),ygrid(1,1):ygrid(end,end));
% for L2 Lcorner
load('LcurveL2-0th.mat')
[fx,fy]=calcForcesFromCoef(forceMesh,fCoeff(:,ireg_corner),xgrid,ygrid,'new');
fmag =  (fx.^2 + fy.^2).^0.5; 
tnormd(:,:,1) = griddata(xgrid,ygrid,fmag,xgrid_fine,ygrid_fine,'linear');
% for L2 optimal
[fx,fy]=calcForcesFromCoef(forceMesh,fCoeff(:,fminIdx),xgrid,ygrid,'new');
fmag =  (fx.^2 + fy.^2).^0.5; 
tnormd(:,:,2) = griddata(xgrid,ygrid,fmag,xgrid_fine,ygrid_fine,'linear');
% for L1 Lcorner
load('LcurveL1-0th.mat')
[fx,fy]=calcForcesFromCoef(forceMesh,msparse(:,ireg_corner),xgrid,ygrid,'new');
fmag =  (fx.^2 + fy.^2).^0.5; 
tnormd(:,:,3) = griddata(xgrid,ygrid,fmag,xgrid_fine,ygrid_fine,'linear');

%% actual calculation for forces at adhesions and noise
for jj=1:nMethods 
    for k=1:nPeaks_FA
        flocMaxFA(k,jj) = tnormd(locMaxI_FA(k,1)-ygrid(1,1)+1,locMaxI_FA(k,2)-xgrid(1,1)+1,jj);
        flocMaxRatioFA(k,jj) = flocMaxFA(k,jj)/fnorm_org(locMaxI_FA(k,1),locMaxI_FA(k,2)); % ratio based on original force norm (0~1)
    end
    for k=1:nPeaks_NA
        flocMaxNA(k,jj) = tnormd(locMaxI_NA(k,1)-ygrid(1,1)+1,locMaxI_NA(k,2)-xgrid(1,1)+1,jj);
        flocMaxRatioNA(k,jj) = flocMaxNA(k,jj)/fnorm_org(locMaxI_NA(k,1),locMaxI_NA(k,2)); % ratio based on original force norm (0~1)
    end
% statistic of flocMaxRatio
    pSR_FA(jj) = mean(flocMaxRatioFA(:,jj));
    pSRerr_FA(jj) = std(flocMaxRatioFA(:,jj))/sqrt(nPeaks_FA);
    pSR_NA(jj) = mean(flocMaxRatioNA(:,jj));
    pSRerr_NA(jj) = std(flocMaxRatioNA(:,jj))/sqrt(nPeaks_NA);
%Background noise
    tnormBG(:,:,jj) = tnormd(:,:,jj).*backgroundMaskfine(ygrid(1,1):ygrid(end,end),xgrid(1,1):xgrid(end,end));
    tempTNBG = tnormBG(:,:,jj);
    tempBGMask = backgroundMaskfine(ygrid(1,1):ygrid(end,end),xgrid(1,1):xgrid(end,end));
    noiseBG(jj) = sum(tempTNBG(:))/sum(tempBGMask(:)); %unit is still Pascal
    noiseBGerr(jj) = std(tempTNBG(:)/sum(tempBGMask(:)));%/sqrt(sum(tempBGMask(:))); %unit is still Pascal
end
    
save([dataPath '/FatFANA.mat'],'flocMaxFA','flocMaxNA','flocMaxRatioFA','flocMaxRatioNA',...
    'pSR_FA','pSR_NA','pSRerr_FA','pSRerr_NA' ,'flocMaxOrgFA','noiseBG','noiseBGerr','tnormBG','tnormd');
