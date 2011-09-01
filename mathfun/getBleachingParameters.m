% Francois Aguet, 08/31/2011

function res = getBleachingParameters(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('Display', 'off', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.parse(data, varargin{:});

nd = length(data);
res(1:nd) = struct('x', [], 'intensity', [], 'model', [], 'k', [], 'a', [], 'd', []);
for i = 1:nd
    if ~(exist([data(i).source 'Detection' filesep 'bleaching.mat'], 'file') == 2) || ip.Results.Overwrite
        res(i) = main(data(i), ip.Results.Display);
    else
        tmp = load([data(i).source 'Detection' filesep 'bleaching.mat']);
        res(i) = tmp.res;
    end
end


function res = main(data, display)

load([data.source 'Detection' filesep 'detection_v2.mat']);

mpath = [data.source 'Detection' filesep 'cellmask.tif'];
if ~(exist(mpath, 'file')==2)
    getCellMask(data);
end

mask = double(imread(mpath));
mask = mask~=0;

nf = data.movieLength;

int_mask = NaN(1,nf);
%int_c = NaN(1,nf);
int_A = NaN(1,nf);
fprintf('Bleaching estimation (%s):     ', getDirFromPath(data.source));
for f = 1:nf
    frame = double(imread(data.framePaths{1}{f}));
    int_mask(f) = mean(frame(mask));
    %int_c(f) = mean(frameInfo(f).c);
    int_A(f) = mean(frameInfo(f).A + frameInfo(f).c);
    fprintf('\b\b\b\b%3d%%', round(100*f/(nf)));
end
fprintf('\n');

x = 1:data.movieLength;

% on mean intensity
% [k a d] = fitDoubleExp(x, int_mask);
% bleachingModel = a(1)*exp(-k(1)*x) + a(2)*exp(-k(2)*x) + d;
% on CCP background
% [k a d] = fitDoubleExp(x, int_c);
% model_c = a(1)*exp(-k(1)*x) + a(2)*exp(-k(2)*x) + d;

[k a d] = fitDoubleExp(x, int_A);
model_A = a(1)*exp(-k(1)*x) + a(2)*exp(-k(2)*x) + d;

%A0 = bleachingModel(1);
A0 = model_A(1);


res.x = x;
res.intensity = int_A;
res.model = model_A;
res.k = k;
res.a = a;
res.d = d;

save([data.source 'Detection' filesep 'bleaching.mat'], 'res');


if strcmpi(display, 'on')
    [~,~] = mkdir([data.source 'Figures']);
    figure;
    hold on;
    %plot(x, int_mask, 'k');
    %plot(x, bleachingModel, 'r');
    %plot(x, int_A, 'k');
    %plot(x, model_A, 'g');
    
    plot(x, int_A/A0, 'k', 'LineWidth', 2);
    plot(x, model_A/A0, 'r', 'LineWidth', 2);
    set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'LineWidth', 1.5);
    
    xlabel('Time (frames)');
    ylabel('Mean CCP intensity (relative)');
    axis([1 data.movieLength 0.7 1.02]);
    print('-depsc2', [data.source 'Figures' filesep 'bleaching.eps']);
    
    dRange = frameInfo(1).dRange{1};
    dRange(2) = dRange(2)*0.8;
    px = data.pixelSize/data.M;
    
    h(1) = figure('Visible', 'off', 'Color', [1 1 1]);
    imagesc(double(imread(data.framePaths{1}{1})));
    caxis(dRange);
    plotScaleBar(5e-6/px);
    axis image off;
    colormap(gray(256)); colorbar;
    print('-depsc', [data.source 'Figures' filesep 'frame1_bleaching.eps']);
    
    h(2) = figure('Visible', 'off', 'Color', [1 1 1]);
    imagesc(double(imread(data.framePaths{1}{end})));
    caxis(dRange);
    plotScaleBar(5e-6/px);
    axis image off;
    colormap(gray(256)); colorbar;
    print('-depsc', [data.source 'Figures' filesep 'frame' num2str(data.movieLength) '_bleaching.eps']);
    
    close(h);
end




function [k a d] = fitDoubleExp(x, fx)

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-6, ...
    'Tolfun', 1e-6);

dy = fx(end);%mean(fx(end-20:end));
mu = sum(x.*(fx-dy))/sum(fx-dy);
init = [1/mu fx(1)-dy 1/mu fx(1)-dy dy];

lb = zeros(1,5);
ub = Inf(1,5);

[p,resnorm,~,~,~,~,J] = lsqnonlin(@cost, init, lb, ub, opts, x, fx);
k = p([1 3]);
a = p([2 4]);
d = p(5);

C = resnorm*full(inv(J'*J));
prmStd = sqrt(diag(C)/(numel(fx)-numel(p) - 1));


% parameters: k1 a1 k2 a2 dy
function v = cost(p, x, fx)
v = p(2)*exp(-p(1)*x) + p(4)*exp(-p(3)*x) + p(5) - fx;

