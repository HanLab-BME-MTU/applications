function [vx,vy]=optical_flow_for_vif(im1, im2, mask1, mask2)
% optical flow analysis for vif structure

% set optical flow parameters (see Coarse2FineTwoFrames.m for the definition of the parameters)
alpha = 0.012;
ratio = 0.75;
minWidth = 20;
nOuterFPIterations = 7;
nInnerFPIterations = 1;
nSORIterations = 30;

para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];

% this is the core part of calling the mexed dll file for computing optical flow
% it also returns the time that is needed for two-frame estimation
tic;
[vx,vy,warpI2] = Coarse2FineTwoFrames(im1,im2,para);
toc

mask12 = (mask1+mask2)>0;

vx(mask12==0)=0;
vy(mask12==0)=0;

h101 = figure(101); hold off; imagesc(im1); colormap(gray);axis image; axis off;
hold on;

vx_desample = vx*0;
vx_desample(1:10:end,1:10:end)=vx(1:10:end,1:10:end);

vy_desample = vy*0;
vy_desample(1:10:end,1:10:end)=vy(1:10:end,1:10:end);

ind = find(abs(vx_desample)>mean2(abs(vx))*1/3 & abs(vy_desample)>mean2(abs(vy))*1/3 ...
    & mask12>0);

[x,y] = meshgrid(1:size(vx,2),1:size(vx,1));

CC  = [  0         0    0.5625
         0         0    0.6250
         0         0    0.6875
         0         0    0.7500
         0         0    0.8125
         0         0    0.8750
         0         0    0.9375
         0         0    1.0000
         0    0.0625    1.0000
         0    0.1250    1.0000
         0    0.1875    1.0000
         0    0.2500    1.0000
         0    0.3125    1.0000
         0    0.3750    1.0000
         0    0.4375    1.0000
         0    0.5000    1.0000
         0    0.5625    1.0000
         0    0.6250    1.0000
         0    0.6875    1.0000
         0    0.7500    1.0000
         0    0.8125    1.0000
         0    0.8750    1.0000
         0    0.9375    1.0000
         0    1.0000    1.0000
    0.0625    1.0000    0.9375
    0.1250    1.0000    0.8750
    0.1875    1.0000    0.8125
    0.2500    1.0000    0.7500
    0.3125    1.0000    0.6875
    0.3750    1.0000    0.6250
    0.4375    1.0000    0.5625
    0.5000    1.0000    0.5000
    0.5625    1.0000    0.4375
    0.6250    1.0000    0.3750
    0.6875    1.0000    0.3125
    0.7500    1.0000    0.2500
    0.8125    1.0000    0.1875
    0.8750    1.0000    0.1250
    0.9375    1.0000    0.0625
    1.0000    1.0000         0
    1.0000    0.9375         0
    1.0000    0.8750         0
    1.0000    0.8125         0
    1.0000    0.7500         0
    1.0000    0.6875         0
    1.0000    0.6250         0
    1.0000    0.5625         0
    1.0000    0.5000         0
    1.0000    0.4375         0
    1.0000    0.3750         0
    1.0000    0.3125         0
    1.0000    0.2500         0
    1.0000    0.1875         0
    1.0000    0.1250         0
    1.0000    0.0625         0
    1.0000         0         0
    0.9375         0         0
    0.8750         0         0
    0.8125         0         0
    0.7500         0         0
    0.6875         0         0
    0.6250         0         0
    0.5625         0         0
    0.5000         0         0];

angle_ind_matrix = nan(size(vx,1),size(vx,2));

angle = atan2(vx./sqrt(vx.^2+vy.^2),vy./sqrt(vx.^2+vy.^2));
angle_ind_matrix = round((angle + pi)/(2*pi)*64);
angle_ind_matrix(angle_ind_matrix<1)=1;
angle_ind_matrix(angle_ind_matrix>64)=64;

angle_ind_matrix_desample = angle_ind_matrix*0;
angle_ind_matrix_desample(1:10:end,1:10:end)=angle_ind_matrix(1:10:end,1:10:end);


for i_c = 1 : 64
    ind_c = find(angle_ind_matrix_desample==i_c);
    quiver(x(ind_c),y(ind_c),vx(ind_c),vy(ind_c),'color',CC(i_c,:));
end

