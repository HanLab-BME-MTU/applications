function [ count, bwmasked ] = bwneighbors( bw, connectivity )
%bwneighbors Counts the number of connected pixels for each pixel
%
% INPUT
% bw a black and white image
% connectivity 4 or 8 depending on the kind of connectivity
%
% OUTPUT
% count, same size as bw, with integer values
%
% See bwNneighbors

if(nargin < 2)
    connectivity = 8;
end

if(connectivity == 8)
%     nfilter = ones(3);
%     nfilter(2,2) = 0;
    nfilter = [ 1 1 1
                1 0 1 
                1 1 1 ];
elseif(connectivity == 4)
%     nfilter = zeros(3);
%     nfilter(:) = ~mod(1:9,2);
    nfilter = [ 0 1 0 
                1 0 1
                0 1 0 ];
end

bw = uint8(bw);

count = imfilter(bw, nfilter);
bwmasked = count .* bw;

end

