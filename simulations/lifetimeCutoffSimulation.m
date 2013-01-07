function lifetimeCutoffSimulation()

% this simulation is not complete, since it does not include cut-off tracks
% at beginning/end

lambda = 80;
k = 2;

% for PDF/CDF calc.
% tmax = lambda*nthroot(-log(1-0.9999), k)

tmin = lambda*nthroot(-log(0.01), k);
fprintf('Time interval to record 99th percentile: %f\n', tmin);

tmax = 250;
framerate = 4;
dt = framerate;
t = dt:dt:tmax;


nFrames = length(t);

% generate theoretical histogram by stepwise integration of the PDF
F = 20; % subsampling factor
dtf = dt/F;
ti = t(1):dtf:(tmax+dt-dtf);
PDF = k./lambda.*(ti/lambda).^(k-1).*exp(-(ti/lambda).^k);
histTH = arrayfun(@(k) sum(PDF(1+(k-1)*F:k*F)), 1:length(t));

% generate PDF for plotting
t_fine = 0:dtf:tmax;
PDF = k./lambda.*(t_fine/lambda).^(k-1).*exp(-(t_fine/lambda).^k);
%CDF = 1 - exp(-(t/lambda).^k);


N = 1e4; % new objects per frame
samples = cell(1,nFrames);
for n = 1:nFrames

    % generate lifetimes for events starting in current frame
    csamples = lambda*nthroot(-log(1-rand(1,N)), k);
    % cut off at lower frame bound
    csamples = csamples - mod(csamples, dtf);
    % throw out short (< 1 frame) tracks
    csamples(csamples==0) = []; 

    endTimes = t(n) + csamples;
    endTimes(endTimes > tmax) = [];
    samples{n} = endTimes - t(n);
end
samples = [samples{:}];


% binning for lifetime histogram
ni = histc(samples, t);


% truncate all histograms at last two frames, then normalize
histTH = histTH(1:end-2);
histTH = histTH/sum(histTH*dt);
ni = ni(1:end-2);
t = t(1:end-2);
t_bc = t(2)-dt/2:dt:t(end)+dt/2;



figure;
subplot(2,1,1);
h = bar(t_bc, ni/sum(ni*dt));
set(h, 'BarWidth', 1, 'FaceColor', [0 0.8 0], 'EdgeColor', [0 0.4 0]);
hold on;
plot(t_bc, histTH, 'r.');
plot(t_fine, PDF, 'r', 'LineWidth', 1);
xlim([0 tmax]);
%title('Measured lifetimes', 'FontName', 'Helvetica', 'FontSize', 16);
%xlabel('Time (s)', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('Relative frequency', 'FontName', 'Helvetica', 'FontSize', 16);
set(gca, 'FontName', 'Helvetica', 'FontSize', 14, 'LineWidth', 1.5);
legend('Measured lifetimes', 'PDF');

m = length(ni);
w = m./(m:-1:1);
wx = w;
% ni = ni.*w;
nix = ni.*w;
nix = nix / sum(nix*dt);
% SE = sum((nix-histTH).^2)

% weights as a function of frame #: w = nFrames ./ (nFrames-(1:nFrames)-1)
% lengths of 'framerate' and 'framerate-1' cannot be measured -> set to zero
w = nFrames ./ (nFrames-2:-1:1);
ni = ni.*w;

% normalize
ni = ni/sum(ni*dt);

% SE = sum((ni-histTH).^2)


subplot(2,1,2);
h = bar(t_bc, ni);
set(h, 'BarWidth', 1, 'FaceColor', [0 0.8 0], 'EdgeColor', [0 0.4 0]);
hold on;
plot(t_bc, ni, '.', 'Color', [0 0.4 0]);
plot(t_bc, histTH, 'ro');
plot(t_fine, PDF, 'r-', 'LineWidth', 1);
xlim([0 tmax]);
%title('Corrected lifetimes', 'FontName', 'Helvetica', 'FontSize', 16);
xlabel('Time (s)', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('Relative frequency', 'FontName', 'Helvetica', 'FontSize', 16);
set(gca, 'FontName', 'Helvetica', 'FontSize', 14, 'LineWidth', 1.5);
legend('Corrected lifetimes', 'PDF');
%print('-depsc2', '-r300', 'lifetimeCutoffCorrection.eps');
