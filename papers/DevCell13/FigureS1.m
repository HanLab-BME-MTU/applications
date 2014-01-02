%=========================================================================================
% Figure S1
%=========================================================================================
% This script generates the panels for Figure S1 from Aguet et al., Dev. Cell, 2013.


f1path = '/Users/aguet/Documents/MATLAB/endocytosis/NCBpaper/Figure S1 - LCa overexpression/';
fset = loadFigureSettings('print');

lg = {' Control', ' CLCa, high O/X', ' CLCa, low O/X', ' CLCb, high O/X', ' CLCb, low O/X'};
xa = [0 2 5 10];

% colors
b = 0.9;
s = 0.5;
cv = hsv2rgb([0 0 0;
              0.3 1 b;
              0.3 s 1;
              0.55 1 b;
              0.55 s 1]);

pos = get(0, 'DefaultFigurePosition');
pos(3:4) = [250 250];      

fset.tfont = [fset.fontName, 'FontSize', 8];
fset.sfont = [fset.fontName, 'FontSize', 10];
fset.lfont = [fset.fontName, 'FontSize', 12];
fset.axOpts = ['Layer', 'top', 'TickDir', 'out', 'LineWidth', 1.5, fset.sfont, 'TickLength', [0.02 0]];

          
%%
% order
idx = [1 3 2 5 4];
          
% HeLa cells
num = xlsread([f1path 'Book1.xlsx'], 'K37:K64');
num(isnan(num)) = [];
num = reshape(num, [4 5])';
num = num(idx,:); % reverse low/high order


figure('Position', pos, 'PaperPositionMode', 'auto');
% axes('Position', [0.25 0.25 0.7 0.7]);
axes('Units', 'pixels', 'Position', [60 60 175 175]);
hold on;
set(gca, fset.axOpts{:}, 'ColorOrder', cv, 'YLim', [0 200], 'XLim', [0 10.2], 'XTick', 0:2:10, 'YTick', 0:25:200);
plot(xa, num, '.-', 'LineWidth', 2, 'MarkerSize', 20);
hl = legend(lg, 'Location', 'SouthEast');
set(hl, fset.tfont{:}, 'Box', 'off');
xlabel('Uptake (min)', fset.lfont{:});
print('-depsc2', '-loose', [f1path 'HeLa_LCa,b_ox.eps']);

%%
% BSC1 cells
[num,txt,raw] = xlsread([f1path 'Book1.xlsx'], 'E3:E30');
num(isnan(num)) = [];
num = reshape(num, [4 5])';
num = num(idx,:); % reverse low/high order

figure('Position', pos, 'PaperPositionMode', 'auto');
axes('Units', 'pixels', 'Position', [60 60 175 175]);
hold on;
set(gca, 'ColorOrder', cv, fset.axOpts{:}, 'YLim', [0 200], 'XLim', [0 10.2], 'XTick', 0:2:10, 'YTick', 0:25:200, 'TickLength', [0.02 0]);
plot(xa, num, '.-', 'LineWidth', 2, 'MarkerSize', 20);
hl = legend(lg, 'Location', 'SouthEast');
set(hl, fset.tfont{:}, 'Box', 'off');
xlabel('Uptake (min)', fset.lfont{:});
% ylabel('Endocytic efficiency', fset.lfont{:});
print('-depsc2', '-loose', [f1path 'BSC1_LCa,b_ox.eps']);

%%
% C- vs. N-terminus label
[num,txt,raw] = xlsread([f1path 'TnfRendo_BSC1clones_multiround_summary_110707.xlsx'], 'L2:M23');

mu = num(:,1); mu(isnan(mu)) = []; mu = reshape(mu, [4 4]);
sigma = num(:,2); sigma(isnan(sigma)) = []; sigma = reshape(sigma, [4 4]);
xa = [0 2.5 5 10];

cv = hsv2rgb([0 0 0;
              0.3 1 b;
              0.99 1 b;
              0.55 1 b;
              0.55 s 1]);
%%

figure('Position', pos, 'PaperPositionMode', 'auto');
axes('Units', 'pixels', 'Position', [60 60 175 175]);
hold on;
idx = [1 2 4 3];
for i = 1:4;
    k = idx(i);
    he = errorbar(xa, mu(:,k), sigma(:,k), 'Color', cv(k,:), 'LineWidth', 1.5);
    setErrorbarStyle(he);
end

for i = 1:4;
    k = idx(i);
    hp(i) = plot(xa, mu(:,k), '.-', 'LineWidth', 2, 'MarkerSize', 20, 'Color', cv(k,:));
end
    
set(gca, fset.axOpts{:}, 'YLim', [0 200], 'XLim', [0 10.2], 'XTick', 0:2:10, 'YTick', 0:25:200, 'TickLength', [0.02 0]);

% lg = {'WT', 'rat LCa o/x', 'mk LCa o/x', 'rat LCa endo.'};
%lg = {'WT', 'rat EGFP-LCa o/x', 'mk LCa-RFP o/x', 'mk LCa-RFP endo.'};
lg = {' Control', ' rat EGFP-CLCa O/X', ' mk CLCa-RFP O/X', ' mk enCLCa-RFP'};

hl = legend(hp, lg(idx), 'Location', 'NorthWest');%, 'Position', [120 220 10 0.0004]);

set(hl, fset.tfont{:}, 'Box', 'off');
xlabel('Uptake (min)', fset.lfont{:});

print('-depsc2', '-loose', [f1path 'BSC1_ENvsOX.eps']);
%%
[num,txt,raw] = xlsread([f1path 'ClcA data.xlsx'], 'M149:N167');
mu = num(:,1); mu(isnan(mu)) = []; mu = reshape(mu, [4 4]);
sigma = num(:,2); sigma(isnan(sigma)) = []; sigma = reshape(sigma, [4 4]);
xa = 0:2:6;

mu = 100*mu;
sigma = 100*sigma;

mu = mu-repmat(mu(1,:), [4 1]);
mu = mu/mu(4,1)*100;
cv = hsv2rgb([0 0 0;
              0.55 0 0.5;
              0.99 1 b;
              0.3 1 b]);

figure('Position', pos, 'PaperPositionMode', 'auto');
axes('Units', 'pixels', 'Position', [60 60 175 175]);
hold on;


idx = [1:4];
for i = 1:4;
    k = idx(i);
    he = errorbar(xa, mu(:,k), sigma(:,k), 'Color', cv(k,:), 'LineWidth', 1.5);
    setErrorbarStyle(he, 0.2*6/10);
end

for i = 1:4;
    k = idx(i);
    hp(i) = plot(xa, mu(:,k), '.-', 'LineWidth', 2, 'MarkerSize', 20, 'Color', cv(k,:));
end
set(gca, fset.axOpts{:}, 'YLim', [0 200], 'XLim', [0 6+6/50], 'XTick', 0:2:10, 'YTick', 0:25:200, 'TickLength', [0.02 0]);


lg = {' Control', ' Untagged CLCa O/X', ' CLCa-EGFP O/X', ' EGFP-CLCa O/X'};
hl = legend(hp, lg, 'Location', 'NorthWest');
set(hl, fset.tfont{:}, 'Box', 'off');
xlabel('Uptake (min)', fset.lfont{:});
print('-depsc2', '-loose', [f1path 'HeLa_NvsC.eps']);


