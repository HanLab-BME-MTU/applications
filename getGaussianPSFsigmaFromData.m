% getGaussianPSFsigmaFromData returns the standard deviation of the Gaussian
% approximation of the PSF estimated from point sources in the data.
%
% INPUTS     data     : data structure
%            {frames} : list of frames to include in the estimation. Default: 1.

% Francois Aguet, September 2010

function sigma = getGaussianPSFsigmaFromData(data, frames)

if nargin<2
    frames = 1;
end

ny = data.imagesize(1);
nx = data.imagesize(2);

frameList = dir([data.source '*.tif*']);

maskPath = [data.source 'Detection' filesep 'Masks' filesep];
maskList = dir([maskPath '*.tif']);

w = 6; % assuming sigma = 2, sufficient for estimate
sigmaCell = cell(1,length(frames));
for f = frames
    % load detection results
    load([data.source 'Detection' filesep 'detectionResults.mat']);
    frame = double(imread([data.source frameList(f).name]));
    mask = double(imread([maskPath maskList(f).name]));
    
    xi = round(frameInfo(1).xcom);
    yi = round(frameInfo(1).ycom);
    % detections within frame bounds
    idx = xi<=w | xi>nx-w | yi<=w | yi>ny-w;
    xi(idx) = [];
    yi(idx) = [];
    np = length(xi);
    
    sigmaVect = zeros(1,np);
    for k = 1:np
        window = frame(yi(k)-w:yi(k)+w, xi(k)-w:xi(k)+w);
        % binary mask
        maskWindow = ~mask(yi(k)-w:yi(k)+w, xi(k)-w:xi(k)+w);
        maskWindow(maskWindow~=0) = 1;
        % background estimate
        c = mean(mean(window(maskWindow)));
        [p] = fitGaussian2D(window, [0 0 max(window(:))-c 1.5 c], 'xyAs');
        sigmaVect(k) = p(4);
    end
    sigmaCell{f} = sigmaVect;
end
sigmaVect = [sigmaCell{:}];
figure;
[ni ti] = hist(sigmaVect);
bar(ti, ni, 'BarWidth', 1, 'FaceColor', [0.8 0 0], 'EdgeColor', [0.4 0 0]);
sigma = mean(sigmaVect);
fprintf('Sigma = %1.3f ± %.3f pixels\n', sigma, std(sigmaVect));