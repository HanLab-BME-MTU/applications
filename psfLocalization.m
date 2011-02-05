% psfLocalization applies Gaussian-approximated PSF fitting at the positions specified in
% the detection structure.
%
% INPUT     data        : experiment structure
%           sigma       : standard deviation of the Gaussian PSF approximation.
%                         Compute using 'getGaussianPSFsigma.m'
%           {overwrite} : 

% Francois Aguet, September 2010

function psfLocalization(data, overwrite)

if nargin<2
    overwrite = 0;
end

h = load([data.source 'Detection' filesep 'detectionResults.mat']);
frameInfo = h.frameInfo;


if ~isfield(frameInfo, 'xloc') || overwrite

    sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{1});
    w = ceil(3*sigma);
    
    %=================================
    % Identify master/slave channels
    %=================================
    nCh = length(data.channels);
    masterChannel = regexp(data.source, data.channels);
    masterChannel = find([masterChannel{:}]);
    slaveChannels = setdiff(1:nCh, masterChannel);

    ny = data.imagesize(1);
    nx = data.imagesize(2);
    
    frameList = cell(1,nCh);
    for k = 1:nCh
        frameList{k} = dir([data.channels{k} '*.tif*']);
    end
    
    maskPath = [data.source 'Detection' filesep 'Masks' filesep];
    maskList = dir([maskPath '*.tif']);
        
    % loop through frames
    frames = zeros(ny,nx,nCh);
    fprintf('Progress:     ');
    for k = 1:data.movieLength
        
        xi = round(frameInfo(k).xcom);
        yi = round(frameInfo(k).ycom);
        np = length(xi);

        xloc = NaN(1,np);
        yloc = NaN(1,np);
        A = NaN(nCh,np);
        c = NaN(nCh,np);
        cStd = NaN(nCh,np);
        
        xStd = NaN(1,np);
        yStd = NaN(1,np);
        aStd = NaN(nCh,np);
        
        hval = NaN(1,np);
        pval = NaN(1,np);
        
        mask = double(imread([maskPath maskList(k).name]));

        for ch = 1:nCh
            frames(:,:,ch) = double(imread([data.channels{ch} frameList{ch}(k).name]));
        end
        
        
        for p = 1:np
            % detections within frame bounds
            if (xi(p)>w && xi(p)<=nx-w && yi(p)>w && yi(p)<=ny-w)
                ch = masterChannel;
                window = frames(yi(p)-w:yi(p)+w, xi(p)-w:xi(p)+w, ch);
                
                % binary mask
                maskWindow = ~mask(yi(p)-w:yi(p)+w, xi(p)-w:xi(p)+w);
                maskWindow(maskWindow~=0) = 1;
                % background estimate
                ci = mean(window(maskWindow));
                [prm, prmStd, ~, res] = fitGaussian2D(window, [0 0 max(window(:))-ci sigma ci], 'xyA');
                pval(p) = res.pval;
                xp = prm(1);
                yp = prm(2);
                if (xp > -w && xp < w && yp > -w && yp < w)
                    xloc(p) = xi(p) + xp;
                    yloc(p) = yi(p) + yp;
                    A(ch,p) = prm(3);
                    xStd(1,p) = prmStd(1);
                    yStd(1,p) = prmStd(2);
                    aStd(ch,p) = prmStd(3);
                    c(ch,p) = ci;
                    cStd(ch,p) = std(window(maskWindow));
                    for ch = slaveChannels
                        window = frames(yi(p)-w:yi(p)+w, xi(p)-w:xi(p)+w, ch);
                        ci = mean(window(maskWindow));
                        [prm prmStd] = fitGaussian2D(window, [xp yp max(window(:))-ci sigma ci], 'A');
                        A(ch,p) = prm(3);
                        aStd(ch,p) = prmStd(1);
                        c(ch,p) = ci;
                        cStd(ch,p) = std(window(maskWindow));
                    end
                end
                
            end
        end
        frameInfo(k).xloc = xloc;
        frameInfo(k).yloc = yloc;
        frameInfo(k).A = A;
        frameInfo(k).x_std = xStd;
        frameInfo(k).y_std = yStd;
        frameInfo(k).A_std = aStd;
        frameInfo(k).sigma = sigma;
        frameInfo(k).c = c;
        frameInfo(k).cStd = cStd;
        frameInfo(k).pval_KS = pval;
        valid = zeros(1,np);
        valid(pval>=0.05) = 1;
        frameInfo(k).valid = valid;
        fprintf('\b\b\b\b%3d%%', round(100*k/(data.movieLength)));
    end
    fprintf('\n');
    
    save([data.source 'Detection' filesep 'detectionResults.mat'], 'frameInfo');
else
    fprintf('Localization has already been performed for dataset ''%s''\n', data.source);
end