%% 062309_cre_CSUX_3
load('/home/sb234/Desktop/FA/062309_cre_CSUX_3/ch488/analysis/movieData.mat');
load([movieData.protrusion.directory filesep movieData.protrusion.samples.fileName]);

map = protrusionSamples.averageNormalComponent;
map_s = sort(map(:));
m1 = map_s(ceil(.01 * numel(map_s)));
M1 = map_s(ceil((1 - .01) * numel(map_s)));
map(map < m1) = m1;
map(map > M1) = M1;

figure('Name', '062309_cre_CSUX_3'), imagesc(map);
colorbar;
xlabel('Frame');
ylabel('Sector #');
title('Protrusion Velocity (pixel/frame)');

print(gcf,'-dtiff', '/home/sb234/Desktop/FA/062309_cre_CSUX_3/ch488/analysis/activityMap.tif');
print(gcf,'-depsc2', '/home/sb234/Desktop/FA/062309_cre_CSUX_3/ch488/analysis/activityMap.eps');

% CSUX560_Adcon_1a_062609
load('/home/sb234/Desktop/FA/CSUX560_Adcon_1a_062609/ch488/analysis/movieData.mat');
load([movieData.protrusion.directory filesep movieData.protrusion.samples.fileName]);

map = protrusionSamples.averageNormalComponent;
map_s = sort(map(:));
m2 = map_s(ceil(.01 * numel(map_s)));
M2 = map_s(ceil((1 - .01) * numel(map_s)));
map(map < m2) = m2;
map(map > M2) = M2;

figure('Name', 'CSUX560_Adcon_1a_062609'), imagesc(map, [m1 M1]);
colorbar;
xlabel('Frame');
ylabel('Sector #');
title('Protrusion Velocity (pixel/frame)');

print(gcf,'-dtiff', '/home/sb234/Desktop/FA/CSUX560_Adcon_1a_062609/ch488/analysis/activityMap.tif');
print(gcf,'-depsc2', '/home/sb234/Desktop/FA/CSUX560_Adcon_1a_062609/ch488/analysis/activityMap.eps');
