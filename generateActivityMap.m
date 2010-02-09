%% CSUX560_Adcon_1a_062609
clear;
load('/home/sb234/Desktop/FA/CSUX560_Adcon_1a_062609/ch488/analysis/movieData.mat');
load([movieData.protrusion.directory filesep movieData.protrusion.samples.fileName]);

map = protrusionSamples.averageNormalComponent;
fprintf('%d nans over %d pixels\n', nnz(isnan(map)), numel(map));
m = min(map(:));
M = max(map(:));

map_n = (1 / (M - m)) * (2 * map - (M + m));

imagesc(map_n);
colorbar;
xlabel('Frame');
ylabel('Sector #');
title('Protrusion Velocity (pixel/frame)');

print(gcf,'-dtiff', '/home/sb234/Desktop/FA/CSUX560_Adcon_1a_062609/ch488/analysis/activityMap.tif');
print(gcf,'-depsc2', '/home/sb234/Desktop/FA/CSUX560_Adcon_1a_062609/ch488/analysis/activityMap.eps');

%% 062309_cre_CSUX_3
clear;
load('/home/sb234/Desktop/FA/062309_cre_CSUX_3/ch488/analysis/movieData.mat');
load([movieData.protrusion.directory filesep movieData.protrusion.samples.fileName]);

map = protrusionSamples.averageNormalComponent;
fprintf('%d nans over %d pixels\n', nnz(isnan(map)), numel(map));

m = min(map(:));
M = max(map(:));

map_n = (1 / (M - m)) * (2 * map - (M + m));

figure, imagesc(map_n);
colorbar;
xlabel('Frame');
ylabel('Sector #');
title('Protrusion Velocity (pixel/frame)');

print(gcf,'-dtiff', '/home/sb234/Desktop/FA/062309_cre_CSUX_3/ch488/analysis/activityMap.tif');
print(gcf,'-depsc2', '/home/sb234/Desktop/FA/062309_cre_CSUX_3/ch488/analysis/activityMap.eps');

clear;