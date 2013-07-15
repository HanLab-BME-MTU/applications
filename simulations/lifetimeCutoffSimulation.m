function lifetimeCutoffSimulation()

N = 1e5; % new objects per frame
% N = 5;

% lifetime distribution is modeled as a gamma function
lambda = 80;
k = 2;
lftPDF = @(lambda, k, t) k./lambda.*(t/lambda).^(k-1).*exp(-(t/lambda).^k);
lftCDF = @(lambda, k, t) 1 - exp(-(t/lambda).^k);

% for PDF/CDF calc.
% tmax = lambda*nthroot(-log(1-0.9999), k)

tmin = lambda*nthroot(-log(0.01), k);
fprintf('Time interval to record 99th percentile: %f\n', tmin);

tmax = 160;
framerate = 4;
dt = framerate;
t = dt:dt:tmax;

nFrames = length(t);

% generate theoretical histogram by stepwise integration of the PDF
histTH = diff(lftCDF(lambda, k, dt:dt:t(end)+dt));

% generate PDF for plotting
t_fine = 0:dt/20:tmax;
PDF = lftPDF(lambda, k, t_fine);


samples = cell(1,nFrames);
for n = 1:nFrames

    % generate lifetimes for events starting in current frame
    csamples = lambda*nthroot(-log(1-rand(1,N)), k);
    % cut off at lower frame bound
    %csamples = csamples - mod(csamples, dtf);
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
xlim([0 tmax+10]);
ylabel('Relative frequency', 'FontName', 'Helvetica', 'FontSize', 16);
set(gca, 'FontName', 'Helvetica', 'FontSize', 14, 'LineWidth', 1.5, 'Layer', 'top');
legend('Measured lifetimes', 'PDF');


% weights as a function of frame #: w = nFrames ./ (nFrames-(1:nFrames)-1)
% lengths of 'framerate' and 'framerate-1' cannot be measured -> set to zero
w = nFrames ./ (nFrames-2:-1:1);
ni = ni.*w;

% normalize
ni = ni/sum(ni*dt);


subplot(2,1,2);
h = bar(t_bc, ni);
set(h, 'BarWidth', 1, 'FaceColor', [0 0.8 0], 'EdgeColor', [0 0.4 0]);
hold on;
plot(t_bc, histTH, 'ro');
plot(t_fine, PDF, 'r-', 'LineWidth', 1);
xlim([0 tmax+10]);
xlabel('Time (s)', 'FontName', 'Helvetica', 'FontSize', 16);
ylabel('Relative frequency', 'FontName', 'Helvetica', 'FontSize', 16);
set(gca, 'FontName', 'Helvetica', 'FontSize', 14, 'LineWidth', 1.5, 'Layer', 'top');
legend('Corrected lifetimes', 'disc. PDF', 'cont. PDF');
%print('-depsc2', '-r300', 'lifetimeCutoffCorrection.eps');

figure; plot(t_bc, ni-histTH);

