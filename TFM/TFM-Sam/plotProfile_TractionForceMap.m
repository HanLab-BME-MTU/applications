%This script is designed to extract the traction maps from TFM movieData
%files and plot a cross section of said traction maps through the
%centerline of the map.
%Also generates ground truth profile plot if ground truth data is provided.

%Experimental Data
load('movieData.mat')
forceField_FTTC = load('TFMPackage/forceField (FTTC)/forceField.mat','forceField');
forceField_FTTC = forceField_FTTC.forceField;
forceField_FastBEM = load('TFMPackage/forceField (FastBEM)/forceField.mat','forceField');
forceField_FastBEM = forceField_FastBEM.forceField;
forceField_FEM = load('TFMPackage/forceField (FEM)/forceField.mat','forceField');
forceField_FEM = forceField_FEM.forceField;
displField = load('TFMPackage/displacementField/displField.mat','displField');
displField = displField.displField;

%Ground Truth
load('groundTruth.mat')

%Reassemble Matrix
%FTTC
forceVecX_FTTC = forceField_FTTC(1).vec(:,1) + forceField_FTTC(1).vec(:,2);
forceVecY_FTTC = forceField_FTTC(2).vec(:,1) + forceField_FTTC(2).vec(:,2);
forceMatX_FTTC = zeros(length(force_x));
forceMatY_FTTC = zeros(length(force_y));
for k = 1:length(forceField_FTTC(1).pos(:,1))
forceMatX_FTTC(forceField_FTTC(1).pos(k,1),forceField_FTTC(1).pos(k,2))...
    = forceVecX_FTTC(k);
forceMatY_FTTC(forceField_FTTC(2).pos(k,1),forceField_FTTC(2).pos(k,2))...
    = forceVecY_FTTC(k);
end
%FastBEM
forceVecX_FastBEM = forceField_FastBEM(1).vec(:,1) + forceField_FastBEM(1).vec(:,2);
forceVecY_FastBEM = forceField_FastBEM(2).vec(:,1) + forceField_FastBEM(2).vec(:,2);
forceMatX_FastBEM = zeros(length(force_x));
forceMatY_FastBEM = zeros(length(force_y));
for k = 1:length(forceField_FastBEM(1).pos(:,1))
forceMatX_FastBEM(forceField_FastBEM(1).pos(k,1),forceField_FastBEM(1).pos(k,2))...
    = forceVecX_FastBEM(k);
forceMatY_FastBEM(forceField_FastBEM(2).pos(k,1),forceField_FastBEM(2).pos(k,2))...
    = forceVecY_FastBEM(k);
end
%FEM
forceVecX_FEM = forceField_FEM(1).vec(:,1) + forceField_FEM(1).vec(:,2);
forceVecY_FEM = forceField_FEM(2).vec(:,1) + forceField_FEM(2).vec(:,2);
forceMatX_FEM = zeros(length(force_x));
forceMatY_FEM = zeros(length(force_y));
for k = 1:length(forceField_FEM(1).pos(:,1))
forceMatX_FEM(forceField_FEM(1).pos(k,1),forceField_FEM(1).pos(k,2))...
    = forceVecX_FEM(k);
forceMatY_FEM(forceField_FEM(2).pos(k,1),forceField_FEM(2).pos(k,2))...
    = forceVecY_FEM(k);
end

%Plotting Results
heatMap_FTTC = generateHeatmapFromField(forceField_FTTC);
heatMap_FastBEM = generateHeatmapFromField(forceField_FastBEM);
heatMap_FEM = generateHeatmapFromField(forceField_FEM);

%Heatmaps
figure(1)
tiledlayout(1,3)
nexttile
imagesc(heatMap_FTTC); cbr1 = colorbar; 
set(cbr1,'YTick',0:50:351)
title('FTTC Force Reconstruction Magnitude')
nexttile
imagesc(heatMap_FastBEM); colorbar;
title('FastBEM Force Reconstruction Magnitude')
nexttile
imagesc(heatMap_FEM); colorbar;
title('FEM Force Reconstruction Magnitude')

%Profile Plot
figure(2)
title('Comparing Traction Profiles of Different Force Reconstruction Methods')
xlabel('Pixel Number')
ylabel('Traction Magnitude')
midPoint = [length(heatMap_FTTC)/2 length(heatMap_FastBEM)/2 length(heatMap_FEM)/2]; %1,2,3: FTTC,FastBEM,FEM
plot(heatMap_FTTC(midPoint(1),:),'LineWidth',2)
hold on
plot(heatMap_FastBEM(midPoint(2),:),'LineWidth',2)
hold on
plot(heatMap_FEM(midPoint(3),:),'LineWidth',2)
legend('FTTC','FastBEM','FEM')
grid on