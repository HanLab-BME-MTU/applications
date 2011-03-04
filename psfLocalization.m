% psfLocalization applies Gaussian-approximated PSF fitting at the positions specified in
% the detection structure.
%
% INPUT     data        : experiment structure
%           sigma       : standard deviation of the Gaussian PSF approximation.
%                         Compute using 'getGaussianPSFsigma.m'
%           {overwrite} : 

% Francois Aguet, September 2010
% Last modified: February 23, 2011

function psfLocalization(data, overwrite)

if nargin<2
    overwrite = 0;
end

h = load([data.source 'Detection' filesep 'detectionResults.mat']);
frameInfo = h.frameInfo;


if ~isfield(frameInfo, 'xloc') || overwrite

    sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{1});
    w2 = ceil(2*sigma);
    w3 = ceil(3*sigma);
    w4 = ceil(4*sigma);
    dw = w4-w3;
    
    [x,y] = meshgrid(-w4:w4);
    r = sqrt(x.^2+y.^2);
    bandMask = zeros(size(r));
    bandMask(r<=w4 & r>=w3) = 1;
    
    [x,y] = meshgrid(-w3:w3);
    r = sqrt(x.^2+y.^2);
    diskMask = zeros(size(r));
    diskMask(r<=w3) = 1;
    
    %=================================
    % Identify master/slave channels
    %=================================
    nCh = length(data.channels);
    masterChannel = find(strcmp(data.source, data.channels));
    slaveChannels = setdiff(1:nCh, masterChannel);

    ny = data.imagesize(1);
    nx = data.imagesize(2);
    
    frameList = cell(1,nCh);
    for k = 1:nCh
        frameList{k} = dir([data.channels{k} '*.tif*']);
    end
    
    maskPath = [data.source 'Detection' filesep 'Masks' filesep];
    maskList = dir([maskPath '*.tif']);

    % indexes for cross-correlation coefficients
    n = 4;
    idx = pcombs(1:n);
    i = idx(:,1);
    j = idx(:,2);
    ij = i+n*(j-1);
    ii = i+n*(i-1);
    jj = j+n*(j-1);
    
    
    % loop through frames
    iRange = zeros(nCh,2);
    frames = zeros(ny,nx,nCh);
    fprintf('Localization progress:     ');
    for k = 1:data.movieLength
        
        xi = round(frameInfo(k).xcom);
        yi = round(frameInfo(k).ycom);
        np = length(xi);

        xloc = NaN(1,np);
        yloc = NaN(1,np);
        A = NaN(nCh,np);
        c = NaN(nCh,np);
        cStd = NaN(nCh,np);
        cStd_mask = NaN(nCh,np);
        
        xStd = NaN(1,np);
        yStd = NaN(1,np);
        aStd = NaN(nCh,np);
        
        pval = NaN(1,np);
        
        mask = double(imread([maskPath maskList(k).name]));
        % binarize
        mask(mask~=0) = 1;
        labels = bwlabel(mask);
        
        for ch = 1:nCh
            frames(:,:,ch) = double(imread([data.channels{ch} frameList{ch}(k).name]));
             % dynamic range
            iRange(ch,:) = [min(min(frames(:,:,ch))) max(max(frames(:,:,ch)))];
        end
        
        
        for p = 1:np
            % detections within frame bounds
            if (xi(p)>w4 && xi(p)<=nx-w4 && yi(p)>w4 && yi(p)<=ny-w4)
                ch = masterChannel;
                
                % label mask
                maskWindow = labels(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4);
                maskWindow(maskWindow==maskWindow(w4+1,w4+1)) = 0;
                
                % estimate background
                cmask = bandMask;
                cmask(maskWindow~=0) = 0;
                window = frames(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4, ch);
                ci = mean(window(cmask==1));
                %cs = std(window(cmask==1));
                % set any other components to NaN
                window(maskWindow~=0) = NaN;
                
                % standard deviation of the background within mask
                bgStd = nanstd(window(bandMask==1));
                
                % reduce to w = 3*sigma from w = 4*sigma
                window = window(dw+1:end-dw, dw+1:end-dw);
                
                % mask w/ 3*sigma disk
                window(diskMask==0) = NaN;
                
                % fit
                [prm, prmStd, C, res] = fitGaussian2D(window, [0 0 max(window(:))-ci sigma ci], 'xyAc');
                K = corrFromC(C,ij,ii,jj);
                
                pval(p) = res.pval;
                xp = prm(1);
                yp = prm(2);
                % eliminate points where localization failed or which are close to image border
                if (xp > -w2 && xp < w2 && yp > -w2 && yp < w2 && prm(3)<2*iRange(ch,2) && max(abs(K(:)))<0.8)
                    xloc(p) = xi(p) + xp; 
                    yloc(p) = yi(p) + yp;
                    A(ch,p) = prm(3);
                    xStd(1,p) = prmStd(1);
                    yStd(1,p) = prmStd(2);
                    aStd(ch,p) = prmStd(3);
                    c(ch,p) = prm(5);
                    cStd(ch,p) = prmStd(4);
                    cStd_mask(ch,p) = bgStd;
                    for ch = slaveChannels
                        window = frames(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4, ch);
                        ci = mean(window(cmask==1));
                        
                        window(maskWindow~=0) = NaN;
                        window = window(dw+1:end-dw, dw+1:end-dw);
                        window(diskMask==0) = NaN;

                        [prm prmStd] = fitGaussian2D(window, [xp yp max(window(:))-ci sigma ci], 'Ac');
                        A(ch,p) = prm(3);
                        aStd(ch,p) = prmStd(1);
                        c(ch,p) = prm(5);
                        cStd(ch,p) = prmStd(2);%std(window(cmask==1));
                    end
                end
                
            end
        end
        
        % localization parameters
        frameInfo(k).xloc = xloc;
        frameInfo(k).x_std = xStd;
        frameInfo(k).yloc = yloc;
        frameInfo(k).y_std = yStd;
        frameInfo(k).A = A;
        frameInfo(k).A_std = aStd;
        frameInfo(k).sigma = sigma;
        frameInfo(k).c = c;
        frameInfo(k).cStd = cStd;
        frameInfo(k).cStd_mask = cStd_mask;

        % test results
        frameInfo(k).pval_KS = pval;
        frameInfo(k).valid_KS = pval>=0.05;
        frameInfo(k).iRange = iRange;
        
        
        % use localization results in tracker
        %frameInfo(k).amp = zeros(np,2);
        %frameInfo(k).xCoord = zeros(np,2);
        %frameInfo(k).yCoord = zeros(np,2);
        
        %frameInfo(k).amp(:,1) = frameInfo(k).A(masterChannel,:);
        %frameInfo(k).xCoord(:,1) = frameInfo(k).xloc;
        %frameInfo(k).yCoord(:,1) = frameInfo(k).yloc;
        
        fprintf('\b\b\b\b%3d%%', round(100*k/(data.movieLength)));
    end
    fprintf('\n');
    
    save([data.source 'Detection' filesep 'detectionResults.mat'], 'frameInfo');
else
    fprintf('Localization has already been performed for dataset ''%s''\n', data.source);
end


function K = corrFromC(C,ij,ii,jj)
n = size(C,1);
K = zeros(n,n);

K(ij) = C(ij) ./ sqrt(C(ii).*C(jj));
% remaining components are not needed
% K = K + K';
% K(sub2ind([n n], 1:n, 1:n)) = 1;
