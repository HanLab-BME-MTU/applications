%This script is designed to extract the traction maps from TFM movieData
%files and plot a cross section of said traction maps through the
%centerline of the map.
%Also generates ground truth profile plot if ground truth data is provided.

close all; clear; 

%Experimental Data
load('movieData.mat')
forceField_FTTC = load('TFMPackage/forceField (FTTC)/forceField.mat','forceField');
forceField_FTTC = forceField_FTTC.forceField;
forceField_FastBEM = load('TFMPackage/forceField (FastBEM)/forceField.mat','forceField');
forceField_FastBEM = forceField_FastBEM.forceField;
forceField_FEM = load('TFMPackage/forceField (FEM)/forceField.mat','forceField');
forceField_FEM = forceField_FEM.forceField;
forceField_nonLinFEM_low = load('TFMPackage/forceField (nonlin FEM)/forceField.mat','forceField');
forceField_nonLinFEM_low = forceField_nonLinFEM_low.forceField;
forceField_nonLinFEM = load('TFMPackage/forceField (nonlin FEM highres)/forceField.mat','forceField');
forceField_nonLinFEM = forceField_nonLinFEM.forceField;
displField = load('TFMPackage/displacementField/displField.mat','displField');
displField = displField.displField;

%Ground Truth
%load('groundTruth.mat')

%Reassemble Matrix
%FTTC
forceVecX_FTTC = forceField_FTTC(1).vec(:,1) + forceField_FTTC(1).vec(:,2);
forceVecY_FTTC = forceField_FTTC(2).vec(:,1) + forceField_FTTC(2).vec(:,2);
%forceMatX_FTTC = zeros(length(force_x));
%forceMatY_FTTC = zeros(length(force_y));
for k = 1:length(forceField_FTTC(1).pos(:,1))
forceMatX_FTTC(forceField_FTTC(1).pos(k,1),forceField_FTTC(1).pos(k,2))...
    = forceVecX_FTTC(k);
forceMatY_FTTC(forceField_FTTC(2).pos(k,1),forceField_FTTC(2).pos(k,2))...
    = forceVecY_FTTC(k);
end
%FastBEM
forceVecX_FastBEM = forceField_FastBEM(1).vec(:,1) + forceField_FastBEM(1).vec(:,2);
forceVecY_FastBEM = forceField_FastBEM(2).vec(:,1) + forceField_FastBEM(2).vec(:,2);
%forceMatX_FastBEM = zeros(length(force_x));
%forceMatY_FastBEM = zeros(length(force_y));
for k = 1:length(forceField_FastBEM(1).pos(:,1))
forceMatX_FastBEM(forceField_FastBEM(1).pos(k,1),forceField_FastBEM(1).pos(k,2))...
    = forceVecX_FastBEM(k);
forceMatY_FastBEM(forceField_FastBEM(2).pos(k,1),forceField_FastBEM(2).pos(k,2))...
    = forceVecY_FastBEM(k);
end
%FEM
forceVecX_FEM = forceField_FEM(1).vec(:,1) + forceField_FEM(1).vec(:,2);
forceVecY_FEM = forceField_FEM(2).vec(:,1) + forceField_FEM(2).vec(:,2);
%forceMatX_FEM = zeros(length(force_x));
%forceMatY_FEM = zeros(length(force_y));
for k = 1:length(forceField_FEM(1).pos(:,1))
forceMatX_FEM(forceField_FEM(1).pos(k,1),forceField_FEM(1).pos(k,2))...
    = forceVecX_FEM(k);
forceMatY_FEM(forceField_FEM(2).pos(k,1),forceField_FEM(2).pos(k,2))...
    = forceVecY_FEM(k);
end
%non-linear FEM
forceVecX_nonLinFEM = forceField_nonLinFEM(1).vec(:,1) + forceField_nonLinFEM(1).vec(:,2);
forceVecY_nonLinFEM = forceField_nonLinFEM(2).vec(:,1) + forceField_nonLinFEM(2).vec(:,2);
%forceMatX_nonLinFEM = zeros(length(force_x));
%forceMatY_nonLinFEM = zeros(length(force_y));
for k = 1:length(forceField_nonLinFEM(1).pos(:,1))
forceMatX_nonLinFEM(forceField_nonLinFEM(1).pos(k,1),forceField_nonLinFEM(1).pos(k,2))...
    = forceVecX_nonLinFEM(k);
forceMatY_nonLinFEM(forceField_nonLinFEM(2).pos(k,1),forceField_nonLinFEM(2).pos(k,2))...
    = forceVecY_nonLinFEM(k);
end
%non-linear FEM low res
forceVecX_nonLinFEM_low = forceField_nonLinFEM_low(1).vec(:,1) + forceField_nonLinFEM_low(1).vec(:,2);
forceVecY_nonLinFEM_low = forceField_nonLinFEM_low(2).vec(:,1) + forceField_nonLinFEM_low(2).vec(:,2);
%forceMatX_nonLinFEM = zeros(length(force_x));
%forceMatY_nonLinFEM = zeros(length(force_y));
for k = 1:length(forceField_nonLinFEM_low(1).pos(:,1))
forceMatX_nonLinFEM_low(forceField_nonLinFEM_low(1).pos(k,1),forceField_nonLinFEM_low(1).pos(k,2))...
    = forceVecX_nonLinFEM_low(k);
forceMatY_nonLinFEM_low(forceField_nonLinFEM_low(2).pos(k,1),forceField_nonLinFEM_low(2).pos(k,2))...
    = forceVecY_nonLinFEM_low(k);
end

%Plotting Results
heatMap_FTTC = generateHeatmapFromField(forceField_FTTC);
heatMap_FastBEM = generateHeatmapFromField(forceField_FastBEM);
heatMap_FEM = generateHeatmapFromField(forceField_FEM);
heatMap_nonLinFEM = generateHeatmapFromField(forceField_nonLinFEM);
heatMap_nonLinFEM_low = generateHeatmapFromField(forceField_nonLinFEM_low);

%Heatmaps
figure
clims = [0 500];
colormap('jet')
subplot(2,3,1)
imagesc(heatMap_FTTC,clims); colorbar; 
title('FTTC Force Magnitude','FontSize',8)
subplot(2,3,2)
imagesc(heatMap_FastBEM,clims); colorbar;
title('FastBEM Force Magnitude','FontSize',8)
subplot(2,3,3)
imagesc(heatMap_FEM,clims); colorbar;
title('FEM Force Magnitude','FontSize',8)
subplot(2,3,4)
imagesc(heatMap_nonLinFEM,clims); colorbar;
title('Non-Linear FEM Force Magnitude','FontSize',8)
subplot(2,3,5)
imagesc(heatMap_nonLinFEM_low,clims); colorbar;
title('Non-Linear FEM Force Magnitude (Low Res)','FontSize',8)


%Profile Plot
figure
title('Comparing Traction Profiles of Different Force Reconstruction Methods')
xlabel('Pixel Number')
ylabel('Traction Magnitude')
midPoint = [length(heatMap_FTTC)/2 length(heatMap_FastBEM)/2 length(heatMap_FEM)/2 ...
    length(heatMap_nonLinFEM)/2 length(heatMap_nonLinFEM_low)/2]; %1,2,3: FTTC,FastBEM,FEM
plot(heatMap_FTTC(midPoint(1),:),'LineWidth',1.5)
hold on
plot(heatMap_FastBEM(midPoint(2),:),'LineWidth',1.5)
hold on
plot(heatMap_FEM(midPoint(3),:),'LineWidth',1.5)
hold on
plot(heatMap_nonLinFEM(midPoint(4),:),'LineWidth',1.5)
hold on
plot(heatMap_nonLinFEM_low(midPoint(5),:),'LineWidth',1.5)
legend('FTTC','FastBEM','FEM','Non-Linear FEM','Non-Linear FEM Low-Res')
grid on