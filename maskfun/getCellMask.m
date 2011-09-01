% Francois Aguet, 08/31/2011

function getCellMask(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('Overwrite', false, @islogical);
ip.parse(data, varargin{:});

for i = 1:length(data)
    if ~(exist([data(i).source 'Detection' filesep 'cellmask.tif'], 'file') == 2) || ip.Results.Overwrite
        computeMask(data(i))
    else
        fprintf('Cell mask has already been computed for %s\n', getShortPath(data(i)));
    end
end



function computeMask(data)
nf = data.movieLength;
proj = zeros(data.imagesize);
fprintf('Computing cell mask for %s:     ', getShortPath(data));
for f = 1:nf
    proj = proj + double(imread(data.framePaths{1}{f}));
    fprintf('\b\b\b\b%3d%%', round(100*f/(nf)));
end
proj = proj/nf;
fprintf('\n');
imwrite(uint16(proj), [data.source 'Detection' filesep 'cellAIP.tif'], 'tif', 'compression' , 'lzw');

% smoothing
proj = filterGauss2D(proj, 2);
try
    t = thresholdFluorescenceImage(proj);
    mask = proj > t;
catch e
    fprintf('No background found (%s).\n', e.identifier);
    mask = true(data.imagesize);
end

CC = bwconncomp(mask,8);
csize = cellfun(@numel, CC.PixelIdxList);
CC.PixelIdxList = CC.PixelIdxList(csize == max(csize));
CC.NumObjects = 1;
mask = labelmatrix(CC);

imwrite(uint8(255*mask), [data.source 'Detection' filesep 'cellmask.tif'], 'tif', 'compression' , 'lzw');

