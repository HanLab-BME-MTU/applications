% Francois Aguet, 12/18/12

function psnr = getPSNRDistribution(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('frameIdx', [], @isnumeric);
ip.parse(data, varargin{:});

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


for k = 1:nf
    f = frameIdx(k);
    
    psnr{k} = 10*log10(frameInfo(f).A(1,:).^2*ni ./ frameInfo(f).RSS);
    
end
psnr = [psnr{:}];
