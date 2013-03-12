% Francois Aguet, 12/18/12

function psnr = getPSNRDistribution(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('frameIdx', [], @isnumeric);
ip.addParamValue('Channel', 1);
ip.addParamValue('OutputMode', [], @(x) strcmpi(x, 'dB'));
ip.addParamValue('Cutoff_f', 5);
ip.addParamValue('Mode', 'max');
ip.parse(data, varargin{:});
ch = ip.Results.Channel;

switch ip.Results.Mode
    case 'max'
        lftData = getLifetimeData(data, 'Cutoff_f', ip.Results.Cutoff_f);
        nd = numel(data);
        
        sigma = getGaussianPSFsigma(data(1).NA, data(1).M, data(1).pixelSize, data(1).markers{1});
        w = ceil(4*sigma);
        ni = (2*w+1)^2; % support used for PSF fit
        
        psnr = cell(1,nd);
        parfor i = 1:numel(data)
            [maxA, idx] = nanmax(lftData(i).A(:,:,ch),[],2);
            [nt,nf,nc] = size(lftData(i).A);
            rssAtMaxA = lftData(i).RSS(sub2ind([nt,nf,nc], (1:nt)', idx, ones(nt,1)));
            psnr{i} = maxA.^2*ni ./ rssAtMaxA;
        end
        %psnr = vertcat(psnr{:});
    otherwise
       
        frameIdx = ip.Results.frameIdx;
        if isempty(frameIdx)
            frameIdx = 1:data.movieLength;
        end
        
        sigma = getGaussianPSFsigma(data.NA, data.M, data.pixelSize, data.markers{1});
        w = ceil(4*sigma);
        ni = (2*w+1)^2; % support used for PSF fit
        
        load([data.source 'Detection' filesep 'detection_v2.mat']);
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
end
