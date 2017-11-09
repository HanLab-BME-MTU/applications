%===============================================================================
% 1) Set up paths
%===============================================================================

% TIRF data (RPE cells)
% pixel size: 65x65 nm

% Channels:
% 488: CLCa-EGFP
% 561: APPL
% 647: EEA1


dpath = '/Users/aguet/Documents/MATLAB/endocytosis/Carlos/4Francois (TIRF APPL EEA1)/140115 FL and AD/FL/';
froot = 'FL_ctrl_TIRF ';
midx = [1 2 3 6 7]; % good cells

chIdx = {'488', '561', '647'};
clear dataFL;
for i = 1:numel(midx)
    for c = 1:3
        dataFL(i).channels{c} = [dpath froot chIdx{c} '_Cell' num2str(midx(i)) '.tif'];
    end
    dataFL(i).results = [dpath froot '_Cell' num2str(midx(i)) '.mat'];
end

dpath = '/Users/aguet/Documents/MATLAB/endocytosis/Carlos/4Francois (TIRF APPL EEA1)/140115 FL and AD/AD/';
froot = 'AD_ctrl_TIRF ';
midx = [1 3:11]; % good cells
clear dataAD
for i = 1:numel(midx)
    for c = 1:3
        dataAD(i).channels{c} = [dpath froot chIdx{c} '_Cell' num2str(midx(i)) '.tif'];
    end
    dataAD(i).results = [dpath froot '_Cell' num2str(midx(i)) '.mat'];
end

%%
%===============================================================================
% 2) Run detection
%===============================================================================
% reset random number generator to ensure reproducibility
rng('default');
processFramesIF(dataFL, 'Overwrite', true);
processFramesIF(dataAD, 'Overwrite', true);

%%
%===============================================================================
% 3) Plot 'cell edge distance' histograms
%===============================================================================
fopts = {'Normalized', false, 'Axis', {[0 10 0 25],[0 10 0 150]},...
    'DisplayMode', 'print', 'names', {' APPL', ' EEA1'}};
outFL = calcDistFromCellEdgeIF(dataFL, 'Name', 'FL', fopts{:});
outAD = calcDistFromCellEdgeIF(dataAD, 'Name', 'AD', fopts{:});

% Clathrin only
foptsc = {'Channels', 1, 'Axis', {[0 10 0 30],[0 10 0 200]},...
    'Hues', 0.55, 'Names', []};
outFLc = calcDistFromCellEdgeIF(dataFL, 'Name', 'FL', fopts{:}, foptsc{:});
outADc = calcDistFromCellEdgeIF(dataAD, 'Name', 'AD', fopts{:}, foptsc{:});

%%
%===============================================================================
% 4) Plot FL & AD cell images
%===============================================================================
% idxFL = 4;  cell6 
idxFL = 1; % cell1
idxAD = 6; % cell7

ch1FL = double(imread(dataFL(idxFL).channels{1}));
ch2FL = double(imread(dataFL(idxFL).channels{2}));
ch3FL = double(imread(dataFL(idxFL).channels{3}));
resFL = load(dataFL(idxFL).results);

ch1AD = double(imread(dataAD(idxAD).channels{1}));
ch2AD = double(imread(dataAD(idxAD).channels{2}));
ch3AD = double(imread(dataAD(idxAD).channels{3}));
resAD = load(dataAD(idxAD).results);

dRange1 = [min(prctile(resFL.ps(1).c,25), prctile(resAD.ps(1).c,25)) ...
    max(prctile(ch1FL(:), 99.5), prctile(ch1AD(:), 99.5))];
dRange2 = [min(prctile(resFL.ps(2).c,25), prctile(resAD.ps(2).c,25)) ...
    max(prctile(ch2FL(:), 99.5), prctile(ch2AD(:), 99.5))];
dRange3 = [min(prctile(resFL.ps(3).c,25), prctile(resAD.ps(3).c,25)) ...
    max(prctile(ch3FL(:), 99.5), prctile(ch3AD(:), 99.5))];

tmp1 = scaleContrast(ch2FL, dRange2);
tmp2 = scaleContrast(ch3FL, dRange3);
rgb = cat(3, tmp1, tmp2, zeros(size(ch2FL)));

[ny,nx] = size(ch2FL);

B = bwboundaries(resFL.mask);
B = vertcat(B{:}); % [y x] coordinates
B(B(:,1)==1 | B(:,1)==ny,:) = [];
B(B(:,2)==1 | B(:,2)==nx,:) = [];

rgb = uint8(rgb);
figure('Units', 'Pixels', 'Position', [150 150 nx/2 ny/2], 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]);
imagesc(rgb); axis image off; colormap(gray(256));
hold on; plot(B(:,2), B(:,1), 'Color', 0.99*[1 1 1], 'LineWidth', 1);
plotScaleBar(5/0.065);

% all channels:
% rgb = cat(3, scaleContrast(ch2FL, dRange2), scaleContrast(ch3FL, dRange3),...
%     scaleContrast(ch1FL, dRange1));
% figure('Units', 'Pixels', 'Position', [150 150 nx/2 ny/2], 'PaperPositionMode', 'auto');
% axes('Position', [0 0 1 1]);
% imagesc(uint8(rgb)); axis image off; colormap(gray(256));
% hold on; plot(B(:,2), B(:,1), 'Color', 0.99*[1 1 1], 'LineWidth', 1);
% plotScaleBar(5/0.065);


% AD cell7
tmp1 = scaleContrast(ch2AD, dRange2);
tmp2 = scaleContrast(ch3AD, dRange3);
rgb = cat(3, tmp1, tmp2, zeros(size(ch2AD)));

[ny,nx] = size(ch2AD);

B = bwboundaries(resAD.mask);
B = vertcat(B{:});
B(B(:,1)==1 | B(:,1)==ny,:) = [];
B(B(:,2)==1 | B(:,2)==nx,:) = [];

rgb = uint8(rgb);
figure('Units', 'Pixels', 'Position', [150 150 nx/2 ny/2], 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]);
imagesc(rgb); axis image off;% colormap(gray(256)); colorbar;
hold on; plot(B(:,2), B(:,1), 'Color', 0.99*[1 1 1], 'LineWidth', 1);
plotScaleBar(5/0.065);

% all channels:
% rgb = cat(3, scaleContrast(ch2AD, dRange2), scaleContrast(ch3AD, dRange3),...
%     scaleContrast(ch1AD, dRange1));
% figure('Units', 'Pixels', 'Position', [150 150 nx/2 ny/2], 'PaperPositionMode', 'auto');
% axes('Position', [0 0 1 1]);
% imagesc(uint8(rgb)); axis image off; colormap(gray(256));
% hold on; plot(B(:,2), B(:,1), 'Color', 0.99*[1 1 1], 'LineWidth', 1);
% plotScaleBar(5/0.065);


%%
%===============================================================================
% 5) Colocalization analysis
%===============================================================================
% Channel order: clathrin, APPL, EEA1

R = 2; % search radius

nd = numel(dataFL);
pctCCPwAPPL = zeros(1,nd);
pctCCPwEEA1 = zeros(1,nd);
pctAPPLwEEA1 = zeros(1,nd);
pctEEA1wAPPL = zeros(1,nd);

pctCCPwAPPLrand = zeros(10,nd);
pctCCPwEEA1rand = zeros(10,nd);
pctAPPLwEEA1rand = zeros(10,nd);
pctEEA1wAPPLrand = zeros(10,nd);
for i = 1:nd
    res = load(dataFL(i).results);
    
    nClat = numel(res.ps(1).x);
    nAPPL = numel(res.ps(2).x);
    nEEA1 = numel(res.ps(3).x);
    
    % colocalization btw CCPs and APPL, % of CCPs
    idx = colocalizationLAP([res.ps(1).x' res.ps(1).y'], [res.ps(2).x' res.ps(2).y'], R);
    pctCCPwAPPL(i) = numel(idx)/nClat;
    
    % pct of CCPs containing EEA1
    idx = colocalizationLAP([res.ps(1).x' res.ps(1).y'], [res.ps(3).x' res.ps(3).y'], R);
    pctCCPwEEA1(i) = numel(idx)/nClat;
    
    % pct of APPL containing EEA1
    idx = colocalizationLAP([res.ps(2).x' res.ps(2).y'], [res.ps(3).x' res.ps(3).y'], R);
    pctAPPLwEEA1(i) = numel(idx)/nAPPL;
    % pct of EEA1 containing APPL
    pctEEA1wAPPL(i) = numel(idx)/nEEA1;
    
    for k = 0:9
        idx = colocalizationLAP([res.ps(1).x' res.ps(1).y'],...
            [res.psRand(2).x(k*nAPPL+(1:nAPPL))' res.psRand(2).y(k*nAPPL+(1:nAPPL))'], R);
        pctCCPwAPPLrand(k+1,i) = numel(idx)/nClat;
        
        idx = colocalizationLAP([res.ps(1).x' res.ps(1).y'],...
            [res.psRand(3).x(k*nEEA1+(1:nEEA1))' res.psRand(3).y(k*nEEA1+(1:nEEA1))'], R);
        pctCCPwEEA1rand(k+1,i) = numel(idx)/nClat;
        
        idx = colocalizationLAP([res.ps(2).x' res.ps(2).y'],...
            [res.psRand(3).x(k*nEEA1+(1:nEEA1))' res.psRand(3).y(k*nEEA1+(1:nEEA1))'], R);
        pctAPPLwEEA1rand(k+1,i) = numel(idx)/nAPPL;
        
        idx = colocalizationLAP([res.ps(3).x' res.ps(3).y'],...
            [res.psRand(2).x(k*nAPPL+(1:nAPPL))' res.psRand(2).y(k*nAPPL+(1:nAPPL))'], R);
        pctEEA1wAPPLrand(k+1,i) = numel(idx)/nEEA1;
    end
    pctCCPwAPPLrand = mean(pctCCPwAPPLrand,1);
    pctCCPwEEA1rand = mean(pctCCPwEEA1rand,1);
    pctAPPLwEEA1rand = mean(pctAPPLwEEA1rand,1);
    pctEEA1wAPPLrand = mean(pctEEA1wAPPLrand,1);
end
samplesFL = {pctCCPwAPPL, pctCCPwEEA1, pctAPPLwEEA1, pctEEA1wAPPL};
muFL = cellfun(@mean, samplesFL);
sFL = cellfun(@std, samplesFL);
samplesFLrand = {pctCCPwAPPLrand, pctCCPwEEA1rand, pctAPPLwEEA1rand, pctEEA1wAPPLrand};
muFLrand = cellfun(@mean, samplesFLrand);
sFLrand = cellfun(@std, samplesFLrand);


nd = numel(dataAD);
pctCCPwAPPL = zeros(1,nd);
pctCCPwEEA1 = zeros(1,nd);
pctAPPLwEEA1 = zeros(1,nd);
pctEEA1wAPPL = zeros(1,nd);

pctCCPwAPPLrand = NaN(10,nd);
pctCCPwEEA1rand = NaN(10,nd);
pctAPPLwEEA1rand = NaN(10,nd);
pctEEA1wAPPLrand = NaN(10,nd);
for i = 1:nd
    res = load(dataAD(i).results);
    
    nClat = numel(res.ps(1).x);
    nAPPL = numel(res.ps(2).x);
    nEEA1 = numel(res.ps(3).x);
    
    % pct of CCPs containing APPL
    idx = colocalizationLAP([res.ps(1).x' res.ps(1).y'], [res.ps(2).x' res.ps(2).y'], R);
    pctCCPwAPPL(i) = numel(idx)/nClat;
    
    % pct of CCPs containing EEA1
    idx = colocalizationLAP([res.ps(1).x' res.ps(1).y'], [res.ps(3).x' res.ps(3).y'], R);
    pctCCPwEEA1(i) = numel(idx)/nClat;
    
    % pct of APPL containing EEA1
    idx = colocalizationLAP([res.ps(2).x' res.ps(2).y'], [res.ps(3).x' res.ps(3).y'], R);
    pctAPPLwEEA1(i) = numel(idx)/nAPPL;
    % pct of EEA1 containing APPL
    pctEEA1wAPPL(i) = numel(idx)/nEEA1;
    
    for k = 0:9
        idx = colocalizationLAP([res.ps(1).x' res.ps(1).y'],...
            [res.psRand(2).x(k*nAPPL+(1:nAPPL))' res.psRand(2).y(k*nAPPL+(1:nAPPL))'], R);
        pctCCPwAPPLrand(k+1,i) = numel(idx)/nClat;
        
        idx = colocalizationLAP([res.ps(1).x' res.ps(1).y'],...
            [res.psRand(3).x(k*nEEA1+(1:nEEA1))' res.psRand(3).y(k*nEEA1+(1:nEEA1))'], R);
        pctCCPwEEA1rand(k+1,i) = numel(idx)/nClat;
        
        idx = colocalizationLAP([res.ps(2).x' res.ps(2).y'],...
            [res.psRand(3).x(k*nEEA1+(1:nEEA1))' res.psRand(3).y(k*nEEA1+(1:nEEA1))'], R);
        pctAPPLwEEA1rand(k+1,i) = numel(idx)/nAPPL;
        
        idx = colocalizationLAP([res.ps(3).x' res.ps(3).y'],...
            [res.psRand(2).x(k*nAPPL+(1:nAPPL))' res.psRand(2).y(k*nAPPL+(1:nAPPL))'], R);
        pctEEA1wAPPLrand(k+1,i) = numel(idx)/nEEA1;
    end
    pctCCPwAPPLrand = mean(pctCCPwAPPLrand,1);
    pctCCPwEEA1rand = mean(pctCCPwEEA1rand,1);
    pctAPPLwEEA1rand = mean(pctAPPLwEEA1rand,1);
    pctEEA1wAPPLrand = mean(pctEEA1wAPPLrand,1);
end

samplesAD = {pctCCPwAPPL, pctCCPwEEA1, pctAPPLwEEA1, pctEEA1wAPPL};
muAD = cellfun(@mean, samplesAD);
sAD = cellfun(@std, samplesAD);
samplesADrand = {pctCCPwAPPLrand, pctCCPwEEA1rand, pctAPPLwEEA1rand, pctEEA1wAPPLrand};
muADrand = cellfun(@mean, samplesADrand);
sADrand = cellfun(@std, samplesADrand);


ncat = 4;
hval = NaN(1,ncat);
pval = NaN(1,ncat);
for i = 1:ncat
    [hval(i), pval(i)] = permTest(samplesFL{i}, samplesAD{i}, 'CmpFunction', @mean);
end
hval
pval

%%
%===============================================================================
% 6) Plot colocalization statistics
%===============================================================================
% bar plot:
setupFigure('DisplayMode', 'print', 'AxesWidth', 8);
barplot2(100*[muFL; muFLrand; muAD; muADrand]', 100*[sFL; sFLrand; sAD; sADrand]',[],[],...
    'FaceColor', hsv2rgb([0.55 1 1; 0.55 0.3 1; 0.1 1 1; 0.1 0.3 1]),...
    'XTickLabel', {'CCS w/ APPL','CCS w/ EEA1', 'APPL w/ EEA1','EEA1 w/ APPL'});
ylabel( '% detected structures', 'FontSize', 12);
hl = legend('FL', 'FL, rand.', '\Delta\alphaAD', '\Delta\alphaAD, rand.', 'Location', 'NorthWest');
set(hl, 'Box', 'off');

% box plot:
% setupFigure('DisplayMode', 'print', 'AxesWidth', 8);
% boxplot2(cellfun(@(i) 100*i, [samplesFL; samplesFLrand; samplesAD; samplesADrand], 'unif', 0),...
%     'DetectOutliers', false, 'YLim', [0 15],...
%     'FaceColor', hsv2rgb([0.55 1 1; 0.55 0.3 1; 0.1 1 1; 0.1 0.3 1]),...
%     'XTickLabel', {'CCS w/ APPL','CCS w/ EEA1', 'APPL w/ EEA1','EEA1 w/ APPL'});
% ylabel( '% detected structures', 'FontSize', 12);
% hl = legend('FL', 'FL, rand.', '\Delta\alphaAD', '\Delta\alphaAD, rand.', 'Location', 'NorthWest');
% set(hl, 'Box', 'off');

setupFigure('DisplayMode', 'print');
barplot2([100*muFL; 100*muAD]', [100*sFL; 100*sAD]',[], [3 4; 5 6; 7 8],...
    'XTickLabel', {'CCS w/ APPL','CCS w/ EEA1', 'APPL w/ EEA1','EEA1 w/ APPL'},...
    'YLim', [0 15], 'YLabel', '% colocalization', 'FaceColor', hsv2rgb([0.55 1 0.9; 0.1 1 1]));
ylabel( '% detected structures', 'FontSize', 12);
hl = legend('FL','\Delta\alphaAD', 'Location', 'NorthWest');
set(hl, 'Box', 'off');
