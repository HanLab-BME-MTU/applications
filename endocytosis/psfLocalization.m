% psfLocalization applies Gaussian-approximated PSF fitting at the positions specified in
% the detection structure.
%
% INPUT     data        : experiment structure
%           sigma       : standard deviation of the Gaussian PSF approximation.
%                         Compute using 'getGaussianPSFsigma.m'
%           {overwrite} : 

% Francois Aguet, September 2010
% Last modified: March 6, 2011

function psfLocalization(data, overwrite)

if nargin<2
    overwrite = 0;
end

h = load([data.source 'Detection' filesep 'detectionResults.mat']);
frameInfo = h.frameInfo;

alpha = 0.05;

kLevel = norminv(1-alpha/2.0, 0, 1); % ~2 std above background

if ~isfield(frameInfo, 'xloc') || overwrite
    
    ny = data.imagesize(1);
    nx = data.imagesize(2);
    
    %=================================
    % Identify master/slave channels
    %=================================
    nCh = length(data.channels);
    masterChannel = find(strcmp(data.source, data.channels));
    slaveChannels = setdiff(1:nCh, masterChannel);
    
    sigmaV = zeros(nCh, 1);
    for k = 1:nCh
        sigmaV(k) = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{k});
    end
    
    maskPath = [data.source 'Detection' filesep 'Masks' filesep];
    maskList = dir([maskPath '*.tif']);
    
    sigma = sigmaV(masterChannel);
    w2 = ceil(2*sigma);
    w3 = ceil(3*sigma);
    w4 = ceil(4*sigma);
    
    [x,y] = meshgrid(-w4:w4);
    r = sqrt(x.^2+y.^2);
    annularMask = zeros(size(r));
    annularMask(r<=w4 & r>=w3) = 1;
    
    % loop through frames
    iRange = zeros(nCh,2);
    frames = zeros(ny,nx,nCh);
    fprintf('Localization progress:     ');
    for k = 1:data.movieLength
        
        xi = round(frameInfo(k).xcom);
        yi = round(frameInfo(k).ycom);
        np = length(xi);

        x = NaN(1,np);
        y = NaN(1,np);
        A = NaN(nCh,np);
        c = NaN(nCh,np);
        
        x_pstd = NaN(1,np);
        y_pstd = NaN(1,np);
        A_pstd = NaN(nCh,np);
        c_pstd = NaN(nCh,np);
        
        sigma_r = NaN(nCh,np);
        SE_sigma_r = NaN(nCh,np);
        pval_Ar = NaN(nCh,np);
        pval_KS = NaN(1,np);
  
        mask = double(imread([maskPath maskList(k).name]));
        % binarize
        mask(mask~=0) = 1;
        labels = bwlabel(mask);
        
        for ch = 1:nCh
            frames(:,:,ch) = double(imread(data.framePaths{ch}{k}));
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
                cmask = annularMask;
                cmask(maskWindow~=0) = 0;
                window = frames(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4, ch);
                
                % initial background value
                ci = mean(window(cmask==1));
                
                % exlude remaining masked components from fit
                window(maskWindow~=0) = NaN;
                
                npx = sum(isfinite(window(:)));
                
                % fit
                [prm, prmStd, ~, res] = fitGaussian2D(window, [0 0 max(window(:))-ci sigma ci], 'xyAc');
                pval_KS(p) = res.pval;
                
                dx = prm(1);
                dy = prm(2);
                % eliminate points where localization failed or which are close to image border
                if (dx > -w2 && dx < w2 && dy > -w2 && dy < w2 && prm(3)<2*diff(iRange(ch,:)))
                    
                    x(p) = xi(p) + dx; 
                    y(p) = yi(p) + dy;
                    A(ch,p) = prm(3);
                    c(ch,p) = prm(5);
                    x_pstd(1,p) = prmStd(1);
                    y_pstd(1,p) = prmStd(2);
                    A_pstd(ch,p) = prmStd(3);
                    c_pstd(ch,p) = prmStd(4);
                    
                    sigma_r(ch,p) = res.std;
                    SE_sigma_r(ch,p) = res.std/sqrt(2*(npx-1));
                    SE_r = SE_sigma_r(ch,p) * kLevel;
                    pval_KS(p) = res.pval;
                    
                    df2 = (npx-1) * (A_pstd(ch,p).^2 + SE_r.^2).^2 ./...
                        (A_pstd(ch,p).^4 + SE_r.^4);
                    scomb = sqrt((A_pstd(ch,p).^2 + SE_r.^2)/npx);
                    T = (A(ch,p) - res.std*kLevel) ./ scomb;
                    pval_Ar(ch,p) = tcdf(T, df2);
                    
                    for ch = slaveChannels
                        window = frames(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4, ch);
                        ci = mean(window(cmask==1));
                        window(maskWindow~=0) = NaN;
                        
                        [prm, prmStd, ~, res] = fitGaussian2D(window, [dx dy max(window(:))-ci sigmaV(ch) ci], 'Ac');
                        A(ch,p) = prm(3);
                        A_pstd(ch,p) = prmStd(1);
                        c(ch,p) = prm(5);
                        c_pstd(ch,p) = prmStd(2);
                        
                        sigma_r(ch,p) = res.std;
                        SE_sigma_r(ch,p) = res.std/sqrt(2*(npx-1));
                        SE_r = SE_sigma_r(ch,p) * kLevel;
                        
                        df2 = (npx-1) * (A_pstd(ch,p).^2 + SE_r.^2).^2 ./...
                            (A_pstd(ch,p).^4 + SE_r.^4);
                        scomb = sqrt((A_pstd(ch,p).^2 + SE_r.^2)/npx);
                        T = (A(ch,p) - res.std*kLevel) ./ scomb;
                        pval_Ar(ch,p) = tcdf(T, df2);                    
                    end
                end
            end
        end
        
        % localization parameters
        frameInfo(k).x = x;
        frameInfo(k).y = y;
        frameInfo(k).A = A;
        frameInfo(k).c = c;        
        frameInfo(k).x_pstd = x_pstd;
        frameInfo(k).y_pstd = y_pstd;
        frameInfo(k).A_pstd = A_pstd;
        frameInfo(k).c_pstd = c_pstd;
        
        frameInfo(k).sigma_r = sigma_r;
        frameInfo(k).SE_sigma_r = SE_sigma_r;
        
        frameInfo(k).pval_KS = pval_KS;
        frameInfo(k).pval_Ar = pval_Ar;
        frameInfo(k).sigma = sigmaV;

        % test results
        frameInfo(k).isPSF = pval_KS >= 0.05;
        frameInfo(k).iRange = iRange;
        
        % to use localization results in tracker
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
