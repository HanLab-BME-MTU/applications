%newDetection(data) detects CCPs using a combination of model-based (PSF) fitting and statistical tests
%
% Inputs:      data : data/movie structure
%     {'overwrite'} : true | {false}

% Francois Aguet, April 2011 (last modified 05/24/2011)

function runDetection(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('overwrite', false, @islogical);
ip.parse(data, varargin{:});
overwrite = ip.Results.overwrite;

for i = 1:length(data)
    if ~(exist([data(i).source 'Detection' filesep 'detection_v2.mat'], 'file') == 2) || overwrite
        fprintf('Running detection for %s ...', getShortPath(data(i)));
        main(data(i));
        fprintf(' done.\n');
    else
        fprintf('Detection has already been run for %s\n', getShortPath(data(i)));
    end
end



function main(data)

% master channel
mCh = strcmp(data.channels, data.source);
nCh = length(data.channels);

sigma = arrayfun(@(k) getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{k}), 1:nCh);

frameInfo(1:data.movieLength) = struct('x', [], 'y', [], 'A', [], 's', [], 'c', [],...
    'x_pstd', [], 'y_pstd', [], 'A_pstd', [], 's_pstd', [], 'c_pstd', [],...
    'x_init', [], 'y_init', [], 'A_mask', [], ...
    'sigma_r', [], 'SE_sigma_r', [], 'RSS', [], 'pval_KS', [], 'pval_Ar', [], 'isPSF', [],...
    'xCoord', [], 'yCoord', [], 'amp', []);

nx = data.imagesize(2);
ny = data.imagesize(1);

fmt = ['%.' num2str(ceil(log10(data.movieLength))) 'd'];
[~,~] = mkdir([data.source 'Detection']);
[~,~] = mkdir([data.source 'Detection' filesep 'Masks']);

%fprintf('Progress:     ');
parfor k = 1:data.movieLength
    img = double(imread(data.framePaths{mCh}{k}));
    
    [pstruct, mask] = pointSourceDetection(img, sigma(mCh));
    np = numel(pstruct.x);
    
    % expand structure for slave channels
    %fnames = {'A', 'A_pstd', 'c', 'c_pstd', 'sigma_r', 'SE_sigma_r', 'RSS', 'pval_Ar'};
    fnames = {'x', 'x_pstd', 'y', 'y_pstd', 'A', 'A_pstd', 'c', 'c_pstd', 'sigma_r', 'SE_sigma_r', 'RSS', 'pval_Ar'};
    for f = 1:length(fnames)
        tmp = NaN(nCh, np);
        tmp(mCh,:) = pstruct.(fnames{f});
        pstruct.(fnames{f}) = tmp;
    end
    
    % retain only mask regions containing localizations    
    CC = bwconncomp(mask);
    labels = labelmatrix(CC);
    loclabels = labels(sub2ind([ny nx], pstruct.y_init, pstruct.x_init));
    idx = setdiff(1:CC.NumObjects, loclabels);
    CC.PixelIdxList(idx) = [];
    CC.NumObjects = length(CC.PixelIdxList);

    % clean mask
    labels = labelmatrix(CC);
    mask = labels~=0;

    % update labels
    loclabels = labels(sub2ind([ny nx], pstruct.y_init, pstruct.x_init));
    
    % get component intensity for each detection
    compInt = cellfun(@(i) sum(img(i))/numel(i), CC.PixelIdxList);
    pstruct.A_mask = compInt(loclabels);   
    
    for ci = setdiff(1:nCh, mCh)
        img = double(imread(data.framePaths{ci}{k}));
        pstructSlave = fitGaussians2D(img, pstruct.x(mCh,:), pstruct.y(mCh,:), [], sigma(ci)*ones(1,np), [], 'Ac');
        
        % localize, and compare intensities & (x,y)-coordinates. Use localization result if it yields better contrast
        pstructSlaveLoc = fitGaussians2D(img, pstruct.x(mCh,:), pstruct.y(mCh,:), pstructSlave.A, sigma(ci)*ones(1,np), pstructSlave.c, 'xyAc');
        idx = sqrt((pstruct.x(mCh,:)-pstructSlaveLoc.x).^2 + (pstruct.y(mCh,:)-pstructSlaveLoc.y).^2) < 3*sigma(mCh) & pstructSlaveLoc.A > pstructSlave.A;
        
        % fill slave channel information
        for f = 1:length(fnames)
            pstruct.(fnames{f})(ci,~idx) = pstructSlave.(fnames{f})(~idx);
            pstruct.(fnames{f})(ci,idx) = pstructSlaveLoc.(fnames{f})(idx);
        end
    end
    
    % add fields for tracker
    pstruct.xCoord = [pstruct.x(mCh,:)' zeros(np,1)];
    pstruct.yCoord = [pstruct.y(mCh,:)' zeros(np,1)];
    pstruct.amp = [pstruct.A(mCh,:)' zeros(np,1)];
    
    frameInfo(k) = orderfields(pstruct, fieldnames(frameInfo(k)));
    
    maskPath = [data.source 'Detection' filesep 'Masks' filesep 'dmask_' num2str(k, fmt) '.tif'];
    imwrite(uint8(255*mask), maskPath, 'tif', 'compression' , 'lzw');
    %fprintf('\b\b\b\b%3d%%', round(100*k/data.movieLength));
end
%fprintf('\n');

save([data.source 'Detection' filesep 'detection_v2.mat'], 'frameInfo', 'data');
