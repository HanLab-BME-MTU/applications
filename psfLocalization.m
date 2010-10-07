% psfLocalization applies Gaussian-approximated PSF fitting at the positions specified in
% the detection structure.
%
% INPUT     data        : experiment structure
%           sigma       : standard deviation of the Gaussian PSF approximation.
%                         Compute using 'getGaussianPSFsigma.m'
%           {overwrite} : 

% Francois Aguet, September 2010

function psfLocalization(data, sigma, overwrite)

if nargin<3
    overwrite = 0;
end
load([data.source 'Detection' filesep 'detectionResults.mat']);

if ~isfield(frameInfo, 'xloc') || overwrite
    
    ny = data.imagesize(1);
    nx = data.imagesize(2);
    
    frameList = dir([data.source '*.tif*']);
    maskPath = [data.source 'Detection' filesep 'Masks' filesep];
    maskList = dir([maskPath '*.tif']);
    
    w = ceil(3*sigma);
    
    % loop through frames
    fprintf('Progress:     ');
    for k = 1:data.movieLength
        
        frame = double(imread([data.source frameList(k).name]));
        mask = double(imread([maskPath maskList(k).name]));
        
        xi = round(frameInfo(k).xcom);
        yi = round(frameInfo(k).ycom);
        
        np = length(xi);
        xloc = zeros(1,np);
        yloc = zeros(1,np);
        A = zeros(1,np);
        cVect = zeros(1,np);
        cStd = zeros(1,np);
        
        for p = 1:np
            % detections within frame bounds
            if xi(p)<=w || xi(p)>nx-w || yi(p)<=w || yi(p)>ny-w
                xloc(p) = NaN;
                yloc(p) = NaN;
                A(p) = NaN;
                cVect(p) = NaN;
            else
                window = frame(yi(p)-w:yi(p)+w, xi(p)-w:xi(p)+w);
                
                % binary mask
                maskWindow = ~mask(yi(p)-w:yi(p)+w, xi(p)-w:xi(p)+w);
                maskWindow(maskWindow~=0) = 1;
                % background estimate
                c = mean(window(maskWindow));
                prm = fitGaussian2D(window, [0 0 max(window(:))-c sigma c], 'xyA');
                xloc(p) = xi(p) + prm(1);
                yloc(p) = yi(p) + prm(2);
                A(p) = prm(3);
                cVect(p) = c;
                cStd(p) = std(window(maskWindow));
            end
        end
        frameInfo(k).xloc = xloc;
        frameInfo(k).yloc = yloc;
        frameInfo(k).A = A;
        frameInfo(k).sigma = sigma;
        frameInfo(k).c = cVect;
        frameInfo(k).cStd = cStd;
        fprintf('\b\b\b\b%3d%%', round(100*k/(data.movieLength)));
    end
    fprintf('\n');
    
    save([data.source 'Detection' filesep 'detectionResults.mat'], 'frameInfo');
else
    fprintf('Localization has already been performed for dataset ''%s''\n', data.source);
end