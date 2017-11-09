function [ux, uy, uMag] = testDisplacementFromForce(forceType,f,E,neu)
% testDisplacementFromForce(forceType,forceMag,E,neu)
% calculate displacement field analytically for given forceMag and E and
% Poisson's ratio neu.

% input:
%       forceType:     'pointForce', 'groupForce'
%       f:    input force magnitude in Pa
%       E:       Young's modulus of the gel in Pa
%       neu:    Poisson's ratio (0 to 0.55)
% output:
%       uMag:      an array of displacement field
% Sangyoon Han August 2014

%% parameter setup
% E=8000;  %Young's modulus, unit: Pa

% neu=0.5;  %Poisson's ratio, only needed for FTTC
% L=0;
% numPoints_out=50;
%% reference image (50x50)
xmax=50;
ymax=50;

% forceType = 'groupForce';

gridSpacing = 1;
xmin = gridSpacing;
ymin = gridSpacing;

[x_mat_u, y_mat_u]=meshgrid(xmin:gridSpacing:xmax,ymin:gridSpacing:ymax);
meshPtsFwdSol = 1024;
% d = 7.8; %corresponding to .25 um2 area with 0.072 um/pix resolutionm
d = 4.7; %corresponding to 0.1 um2 area with 0.072 um/pix resolutionm
[ux, uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,...
    @(x,y) assumedForceAniso2D(1,x,y,25,25,0,f,d,d,forceType),...
    @(x,y) assumedForceAniso2D(2,x,y,25,25,0,f,d,d,forceType),'fft',[],meshPtsFwdSol,500,neu,false);

% figure, quiver(x_mat_u,y_mat_u,ux,uy)
uMag = (ux.^2+uy.^2).^0.5;
return
%% evaluation
E=1000;
neu=0.5;
forceType = 'groupForce';
f = 100;
uMag = testDisplacementFromForce(forceType,f,E,neu);