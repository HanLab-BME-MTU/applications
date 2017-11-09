%% single force experiment
% input parameters to be replaced with function inputs
f=500; %Pa
d=8;
%% reference image (300x200)
xmax=200;
ymax=300;

%% Now displacement field from given force
E=8000;  %Young's modulus, unit: Pa
meshPtsFwdSol=2^10;
forceType = 'groupForce';

gridSpacing = 1;
xmin = gridSpacing;
ymin = gridSpacing;

[x_mat_u, y_mat_u]=meshgrid(xmin:gridSpacing:xmax,ymin:gridSpacing:ymax);

[ux, uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,...
    @(x,y) assumedForceAniso2D(1,x,y,(100),150,0,f,d,d,forceType),...
    @(x,y) assumedForceAniso2D(2,x,y,(100),150,0,f,d,d,forceType),'fft',[],meshPtsFwdSol); %,'conv',[],meshPtsFwdSol);
force_x = assumedForceAniso2D(1,x_mat_u,y_mat_u,(100),150,0,f,d,d,forceType);
force_y = assumedForceAniso2D(2,x_mat_u,y_mat_u,(100),150,0,f,d,d,forceType);
%% finite thickness solution - 34 um
h = 34/0.072;
[ux_fin34, uy_fin34]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,...
    @(x,y) assumedForceAniso2D(1,x,y,(100),150,0,f,d,d,forceType),...
    @(x,y) assumedForceAniso2D(2,x,y,(100),150,0,f,d,d,forceType),'fft_finite',[],meshPtsFwdSol,h); %,'conv',[],meshPtsFwdSol);
%% finite thickness solution - 34 um
h = 10/0.072;
[ux_fin10, uy_fin10]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,...
    @(x,y) assumedForceAniso2D(1,x,y,(100),150,0,f,d,d,forceType),...
    @(x,y) assumedForceAniso2D(2,x,y,(100),150,0,f,d,d,forceType),'fft_finite',[],meshPtsFwdSol,h); %,'conv',[],meshPtsFwdSol);
%% finite thickness solution - 34 um
h = 1/0.072;
[ux_fin1, uy_fin1]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,...
    @(x,y) assumedForceAniso2D(1,x,y,(100),150,0,f,d,d,forceType),...
    @(x,y) assumedForceAniso2D(2,x,y,(100),150,0,f,d,d,forceType),'fft_finite',[],meshPtsFwdSol,h); %,'conv',[],meshPtsFwdSol);
%% plot
figure, subplot(3,2,1); imshow(force_y(120:180,70:130),[0 500]), colormap jet
subplot(3,2,3); imshow(uy(120:180,70:130),[0 0.59]), colormap jet
subplot(3,2,4); imshow(uy_fin34(120:180,70:130),[0 0.59]), colormap jet
subplot(3,2,5); imshow(uy_fin10(120:180,70:130),[0 0.59]), colormap jet
subplot(3,2,6); imshow(uy_fin1(120:180,70:130),[0 0.59]), colormap jet
%% profile
figure, plot(1:300,uy(:,100),'k'), hold on
plot(1:300,uy_fin34(:,100),'--', 'Color',[1 .5 0]),
plot(1:300,uy_fin10(:,100),'-.', 'Color',[0 153/255 76/255]),
plot(1:300,uy_fin1(:,100),'b:'),
legend('Boussinesq','34 um','10 um','1 um')
xlabel('x (pixel)')
ylabel('displacement  uy (pixel)')
%% Only G11
x=1:100; y=0; h = 34/0.072;
G11 = boussinesqGreens(1,1,x,y,E);
G11finite34 = finiteThicknessGreens(1,1,x,y,E,h);
h = 10/0.072;
G11finite10 = finiteThicknessGreens(1,1,x,y,E,h);
h = 1/0.072;
G11finite1 = finiteThicknessGreens(1,1,x,y,E,h);
figure, plot(x, G11,'k'), hold on, 
plot(x, G11finite34,'--', 'Color',[1 .5 0])
plot(x, G11finite10,'-.', 'Color',[0 153/255 76/255])
plot(x, G11finite1,'b:')
legend('Boussinesq','34 um','10 um','1 um')
xlabel('x (pixel)')
ylabel('G_{11}')
%% save
save('allData.mat')