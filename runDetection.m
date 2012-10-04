%runDetection(data) detects CCPs using a combination of model-based (PSF) fitting and statistical tests
%
% Inputs:      data : data/movie structure
%         {'Sigma'} : standard deviation of the Gaussian used for fitting
%     {'Overwrite'} : true | {false}

% Francois Aguet, April 2011 (last modified 05/24/2011)

function runDetection(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('Sigma', [], @(x) numel(x)==length(data(1).channels));
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('Master', [], @isnumeric);
ip.parse(data, varargin{:});
overwrite = ip.Results.Overwrite;
mCh = ip.Results.Master;
if isempty(mCh)
    mCh = unique(arrayfun(@(i) find(strcmpi(i.channels, i.source)), data));
    if numel(mCh)>1
        error('Master channel index is inconsistent accross data sets.');
    end
end

for i = 1:length(data)
    if ~(exist([data(i).channels{mCh} 'Detection' filesep 'detection_v2.mat'], 'file') == 2) || overwrite
        fprintf('Running detection for %s ...', getShortPath(data(i)));
        main(data(i), ip.Results.Sigma, mCh);
        fprintf(' done.\n');
    else
        fprintf('Detection has already been run for %s\n', getShortPath(data(i)));
    end
end



function main(data, sigma, mCh)

% master channel
nCh = length(data.channels);

if isempty(sigma)
    sigma = arrayfun(@(k) getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{k}), 1:nCh);
end

frameInfo(1:data.movieLength) = struct('x', [], 'y', [], 'A', [], 's', [], 'c', [],...
    'x_pstd', [], 'y_pstd', [], 'A_pstd', [], 'c_pstd', [],...
    'x_init', [], 'y_init', [], 'maskA', [], 'maskN', [],...
    'sigma_r', [], 'SE_sigma_r', [], 'RSS', [], 'pval_Ar', [],  'mask_Ar', [], 'hval_Ar', [],  'hval_AD', [], 'isPSF', [],...
    'xCoord', [], 'yCoord', [], 'amp', [], 'dRange', []);

nx = data.imagesize(2);
ny = data.imagesize(1);

fmt = ['%.' num2str(ceil(log10(data.movieLength+1))) 'd'];
[~,~] = mkdir([data.channels{mCh} 'Detection']);
[~,~] = mkdir([data.channels{mCh} 'Detection' filesep 'Masks']);

% double fields, multi-channel
dfields = {'x', 'y', 'A', 'c', 'x_pstd', 'y_pstd', 'A_pstd', 'c_pstd', 'sigma_r', 'SE_sigma_r', 'RSS', 'pval_Ar'};
% logical fields, multi-channel
lfields = {'hval_Ar', 'hval_AD', 'isPSF'};
% slave channel fields
sfields = [dfields {'hval_Ar', 'hval_AD'}]; % 'isPSF' is added later

rmfields = [dfields lfields {'x_init', 'y_init', 'maskA', 'maskN', 'mask_Ar'}];

parfor k = 1:data.movieLength
    img = double(imread(data.framePaths{mCh}{k})); %#ok<PFBNS>
    
    [pstruct, mask] = pointSourceDetection(img, sigma(mCh)); %#ok<PFBNS>
    
    if ~isempty(pstruct)
        pstruct.s = sigma;
        pstruct = rmfield(pstruct, 's_pstd');
        
        pstruct.dRange{mCh} = [min(img(:)) max(img(:))];
        np = numel(pstruct.x);
        
        % expand structure for slave channels
        for f = 1:length(dfields)
            tmp = NaN(nCh, np);
            tmp(mCh,:) = pstruct.(dfields{f});
            pstruct.(dfields{f}) = tmp;
        end
        for f = 1:length(lfields)
            tmp = false(nCh, np);
            tmp(mCh,:) = pstruct.(lfields{f});
            pstruct.(lfields{f}) = tmp;
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
        
        % get component size and intensity for each detection
        compSize = cellfun(@(i) numel(i), CC.PixelIdxList);
        pstruct.maskN = compSize(loclabels);
        compInt = cellfun(@(i) sum(img(i))/numel(i), CC.PixelIdxList);
        pstruct.maskA = compInt(loclabels);
        
        for ci = setdiff(1:nCh, mCh)
            img = double(imread(data.framePaths{ci}{k}));
            pstruct.dRange{ci} = [min(img(:)) max(img(:))];
            pstructSlave = fitGaussians2D(img, pstruct.x(mCh,:), pstruct.y(mCh,:), [], sigma(ci)*ones(1,np), [], 'Ac');
            
            % localize, and compare intensities & (x,y)-coordinates. Use localization result if it yields better contrast
            pstructSlaveLoc = fitGaussians2D(img, pstruct.x(mCh,:), pstruct.y(mCh,:), pstructSlave.A, sigma(ci)*ones(1,np), pstructSlave.c, 'xyAc');
            idx = sqrt((pstruct.x(mCh,:)-pstructSlaveLoc.x).^2 + (pstruct.y(mCh,:)-pstructSlaveLoc.y).^2) < 3*sigma(mCh) & pstructSlaveLoc.A > pstructSlave.A;
            
            % fill slave channel information
            for f = 1:length(sfields)
                pstruct.(sfields{f})(ci,~idx) = pstructSlave.(sfields{f})(~idx);
                pstruct.(sfields{f})(ci,idx) = pstructSlaveLoc.(sfields{f})(idx);
            end
            
            nanIdx = isnan(pstructSlave.x); % points within slave channel border, remove from detection results
            for f = 1:length(rmfields)
                pstruct.(rmfields{f})(:,nanIdx) = [];
            end
            np = size(pstruct.x,2);
            
            pstruct.isPSF(ci,:) = ~pstruct.hval_AD(ci,:);
        end
        
        % add fields for tracker
        pstruct.xCoord = [pstruct.x(mCh,:)' pstruct.x_pstd(mCh,:)'];
        pstruct.yCoord = [pstruct.y(mCh,:)' pstruct.y_pstd(mCh,:)'];
        pstruct.amp =    [pstruct.A(mCh,:)' pstruct.A_pstd(mCh,:)'];
        frameInfo(k) = orderfields(pstruct, fieldnames(frameInfo(k))); %#ok<PFOUS>
    else
        frameInfo(k).dRange{mCh} = [min(img(:)) max(img(:))];
        for ci = setdiff(1:nCh, mCh)
            img = double(imread(data.framePaths{ci}{k}));
            frameInfo(k).dRange{ci} = [min(img(:)) max(img(:))];
        end
    end
    
    maskPath = [data.channels{mCh} 'Detection' filesep 'Masks' filesep 'dmask_' num2str(k, fmt) '.tif'];
    imwrite(uint8(255*mask), maskPath, 'tif', 'compression' , 'lzw');
end

save([data.channels{mCh} 'Detection' filesep 'detection_v2.mat'], 'frameInfo');
