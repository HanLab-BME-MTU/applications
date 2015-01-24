function [label]=markedWatershed(vol,scales,thresh,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('vol', @isnumeric);
ip.addRequired('scales', @isnumeric);
ip.addRequired('thresh', @isnumeric);
ip.addParamValue('WindowSize', []); % Default: 2*sigma, see fitGaussians3D.
ip.parse(vol, scales,thresh, varargin{:});

if ~isa(vol, 'double')
    vol = double(vol);
end

if numel(scales)==1
    scales = [scales scales];
end

ws = ip.Results.WindowSize;
if isempty(ws)
    ws = ceil(2*scales);
elseif numel(ws)==1
    ws = [ws ws];
end

%-------------------------------------------------------------------------------------------
% Convolutions
%-------------------------------------------------------------------------------------------
% right-hand side of symmetric kernels
gx = exp(-(0:ws(1)).^2/(2*scales(1)^2));
gz = exp(-(0:ws(2)).^2/(2*scales(2)^2));
fg = conv3fast(vol, gx, gx, gz);
fu =  conv3fast(vol,    ones(1,ws(1)+1), ones(1,ws(1)+1), ones(1,ws(2)+1));
fu2 = conv3fast(vol.^2, ones(1,ws(1)+1), ones(1,ws(1)+1), ones(1,ws(2)+1));

% Laplacian of Gaussian-filtered input
gx2 = (0:ws(1)).^2 .*gx;
gz2 = (0:ws(2)).^2 .*gz;
fgx2 = conv3fast(vol, gx2, gx, gz);
fgy2 = conv3fast(vol, gx, gx2, gz);
fgz2 = conv3fast(vol, gx, gx, gz2);
imgLoG = (2/scales(1)^2+1/scales(2)^2)*fg - ((fgx2+fgy2)/scales(1)^4 + fgz2/scales(2)^4);
clear fgx2 fgy2 fgz2;

imseriesshow(imgLoG)

label=watershed(-imgLoG);
label(smooth3(vol,'gaussian',[3 3 3],scales(1))<ip.Results.thresh)=0;
pstruct=[];
