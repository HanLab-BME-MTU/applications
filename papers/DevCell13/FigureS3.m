%====================================================
% Figure S3: Illustration of multi-step processes
%====================================================
fset = loadFigureSettings('print');
fpath = '/Users/aguet/Documents/MATLAB/endocytosis/CMEpaper/Figure S3 - Multi-step kinetics/';

gammaFct = @(t, k, n) k^n * exp(-k*t) .* t.^(n-1) / gamma(n);
t_fine = 0:0.1:150;

%----------------------------------------------------
% Panel A: exponential w/ var. rates
%----------------------------------------------------
figure(fset.fOpts{:})
axes(fset.axOpts{:});
hold on;
kv = [0.05 0.1 0.2];
ltext = cell(1,3);

h = linspace(0.55, 0.6, 3);
v = ones(1,3);
s = 0.6:0.2:1;
cmap = hsv2rgb([h' s' v']);

for i = 1:numel(kv)
    plot(t_fine, gammaFct(t_fine, kv(i), 1), '-', 'Color', cmap(3-i+1,:), 'LineWidth', 1);
    ltext{i} = [' k = ' num2str(kv(i)) ' s^{-1}'];
end
hl = legend(ltext{:});
set(hl, 'Box', 'off', fset.sfont{:}, 'Position', [5.5 4 1.5 1.25]);
axis([0 120 0 0.06]);
xlabel('Time (s)', fset.lfont{:});
ylabel('Frequency', fset.lfont{:});
% print('-depsc2', '-loose', [fpath 'multistep_exp.eps']);

%%
%----------------------------------------------------
% Panel B: Gamma distribution w/ var. steps
%----------------------------------------------------
k = 0.1;
h = linspace(0.25, 0.35, 5);
s = linspace(0.3, 0.9, 5);
v = linspace(0.3, 0.9, 5);
cmap = hsv2rgb([h' s' v']);

figure(fset.fOpts{:})
axes(fset.axOpts{:});
hold on;
for i = 1:5
    plot(t_fine, gammaFct(t_fine, k, i), '-', 'Color', cmap(i,:), 'LineWidth', 1);
end
axis([0 120 0 0.06]);
xlabel('Time (s)', fset.lfont{:});
ylabel('Frequency', fset.lfont{:});
% print('-depsc2', '-loose', [fpath 'multistep_gamma.eps']);

% Graphical abstract version
k = 0.1;
h = linspace(0.55, 0.6, 5);
s = linspace(0.2, 0.9, 5);
v = linspace(1, 1, 5);
cmap = hsv2rgb([h' s' v']);

figure(fset.fOpts{:})
axes(fset.axOpts{:}, 'Position', [1.5 1.5 6 2.5]);
hold on;
for i = 1:5
    plot(t_fine, gammaFct(t_fine, k, i), '-', 'Color', cmap(i,:), 'LineWidth', 1);
end
axis([0 120 0 0.06]);
xlabel('Lifetime (s)', fset.lfont{:});
ylabel('Frequency', fset.lfont{:});
% print('-depsc2', '-loose', ['multistep_gammaGA.eps']);


%%
%----------------------------------------------------
% Panel C: Gamma function w/ 3 steps and var. rate
%----------------------------------------------------
h = linspace(0.55, 0.6, 3);
v = ones(1,3);
s = 0.6:0.2:1;
cmap = hsv2rgb([h' s' v']);
kv = [0.05 0.1 0.2];

figure(fset.fOpts{:})
axes(fset.axOpts{:});
hold on;
for i = 1:numel(kv)
    plot(t_fine, gammaFct(t_fine, kv(i), 3), '-', 'Color', cmap(3-i+1,:), 'LineWidth', 1);
    ltext{i} = [' k = ' num2str(kv(i)) ' s^{-1}'];
end
hl = legend(ltext{:});
set(hl, 'Box', 'off', fset.sfont{:}, 'Position', [5.5 4 1.5 1.25]);
axis([0 120 0 0.06]);
xlabel('Time (s)', fset.lfont{:});
ylabel('Frequency', fset.lfont{:});

% print('-depsc2', '-loose', [fpath 'multistep_gamma_var.eps']);
