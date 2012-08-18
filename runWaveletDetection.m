% Runs the wavelet-based detection algorithm used in the Loerke et al. 2009 paper
% The 'PostProcLevel' parameter is the same as the old 'iclean' parameter, set to 1 by default.

% Francois Aguet, 08/18/12

function runWaveletDetection(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('FileName', 'detection_v1.mat');
ip.addParamValue('PostProcLevel', 1);
ip.parse(data, varargin{:});
overwrite = ip.Results.Overwrite;

for i = 1:numel(data)
    if ~(exist([data(i).source 'Detection' filesep ip.Results.FileName], 'file') == 2) || overwrite
        fprintf('Running wavelet detection for %s ...', getShortPath(data(i)));
        main(data(i), ip.Results.FileName, ip.Results.PostProcLevel);
        fprintf(' done.\n');
    else
        fprintf('Detection has already been run for %s\n', getShortPath(data(i)));
    end
end


function main(data, fileName, postProcLevel)

mpath = [data.source 'Detection' filesep 'WaveletMasks' filesep];
[~,~] = mkdir(mpath);

fmt = ['%.' num2str(ceil(log10(data.movieLength))) 'd'];

nf = data.movieLength;

frameInfo(1:nf) = struct('ymax', [], 'xmax', [], 'inn', [], 'yav', [], 'xav', [], 'intot', [],...
    'csize', [], 'lxm', [], 'labl', [], 'num', [], 'nmax', [], 'xCoord', [], 'yCoord', [], 'amp', []);
parfor f = 1:nf

    frame = double(imread(data.framePaths{1}{f})); %#ok<PFBNS>

    [iFrameInfo mask] = main283AUTO_standalone(frame, postProcLevel);
    Z = zeros(iFrameInfo.num,1);
    iFrameInfo.xCoord = [iFrameInfo.xav Z];
    iFrameInfo.yCoord = [iFrameInfo.yav Z];
    intot = zeros(iFrameInfo.num,1);
    intot(iFrameInfo.labl) = iFrameInfo.intot;
    iFrameInfo.amp = [intot Z];
    frameInfo(f) = iFrameInfo; %#ok<PFOUS>
    maskPath = [data.source 'Detection' filesep 'WaveletMasks' filesep 'wmask_' num2str(f, fmt) '.tif'];
    imwrite(mask, maskPath, 'tif', 'compression' , 'lzw');    
end
save([data.source 'Detection' filesep fileName], 'frameInfo');
