%[psnr] = getPSNRDistribution(data, varargin) returns the PSNR distribution for the detections in 'data'.
% The PSNR is calculated as A^2*ni/RSS where 'ni' is the number of pixels in the support
% used for fitting.

% Francois Aguet, 12/18/12

function [psnr] = getPSNRDistribution(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('frameIdx', [], @isnumeric);
ip.addParamValue('Channel', 1);
ip.addParamValue('OutputMode', [], @(x) strcmpi(x, 'dB'));
ip.addParamValue('Cutoff_f', 5);
ip.addParamValue('Mode', 'all', @(x) any(strcmpi(x, {'all', 'max'})));
ip.parse(data, varargin{:});
ch = ip.Results.Channel;

switch ip.Results.Mode
    case 'all'
        frameIdx = ip.Results.frameIdx;
        if isempty(frameIdx)
            frameIdx = 1:data.movieLength;
        end
        
        load([data.source 'Detection' filesep 'detection_v2.mat']);
        w = ceil(4*frameInfo(1).s(1));
        ni = (2*w+1)^2; % support used for PSF fit
        
        nf = numel(frameIdx);
        psnr = cell(1,nf);
        
        parfor k = 1:nf
            f = frameIdx(k);
            psnr{k} = frameInfo(f).A(ch,:).^2*ni ./ frameInfo(f).RSS(ch,:);
        end
        psnr = [psnr{:}];
        if strcmpi(ip.Results.OutputMode, 'dB')
            psnr = 10*log10(psnr);
        end
    case 'max'
        sigma = getGaussianPSFsigma(data(1).NA, data(1).M, data(1).pixelSize, data(1).markers{1});
        w = ceil(4*sigma);
        ni = (2*w+1)^2; % support used for PSF fit

        lftData = getLifetimeData(data, 'Cutoff_f', ip.Results.Cutoff_f);
        nd = numel(data);
        
        psnr = cell(1,nd);
        parfor i = 1:numel(data)
            [maxA, idx] = nanmax(lftData(i).A(:,:,ch),[],2);
            [nt,nf,nc] = size(lftData(i).A);
            rssAtMaxA = lftData(i).RSS(sub2ind([nt,nf,nc], (1:nt)', idx, ones(nt,1)));
            psnr{i} = maxA.^2*ni ./ rssAtMaxA;
        end
end
