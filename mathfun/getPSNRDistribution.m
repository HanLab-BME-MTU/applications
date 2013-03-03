% Francois Aguet, 12/18/12

function psnr = getPSNRDistribution(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('frameIdx', [], @isnumeric);
ip.addParamValue('Channel', 1);
ip.addParamValue('OutputMode', [], @(x) strcmpi(x, 'dB'));
ip.parse(data, varargin{:});
ch = ip.Results.Channel;

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
