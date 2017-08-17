function [ maxima, coords, response] = linearizeMaxima( response, maxima, K )
%linearizeMaxima Linearize maxima by getting rid of NaNs and replicating
%the response as needed
%
% INPUT
% maxima - orientation local maxima in R x C x M form with NaNs
% response - OrientationSpaceResponse object or response matrix
%          - R x C x N (N usually is 2*K+1)
% K - scalar or R x C matrix
%
% OUTPUT
% maxima - orientation local maxima 1 x MRC, no NaNs
% coords - struct with the following fields
% .r - row number, 1 x MRC
% .c - column number, 1 x MRC
% .m - local maxima number, 1 x MRC
% .K - 1 x MRC, indicating K parameter for the response
% response - response matrix N x MRC (replicating response for each maxima)

if(isa(response,'OrientationSpaceResponse'))
    R = response;
    response = real(R.a);
    if(nargin < 2)
        maxima = R.getRidgeOrientationLocalMaxima;
        maxima = maxima*2;
    end
    if(nargin < 3)
        K = repmat(R.filter.K,size(maxima));
    end
end
if(nargin > 2)
    if(isscalar(K))
        K = repmat(K,size(maxima));
    end
end
        
    

maxima = shiftdim(maxima,2);
maxima_size = size(maxima);

coords.r = repmat(1:maxima_size(2),[maxima_size(1) 1 maxima_size(3)]);
coords.c = repmat(shiftdim(1:maxima_size(3),-1),[maxima_size(1) maxima_size(2) 1]);
coords.m = repmat(shiftdim(1:maxima_size(1),1),[1 maxima_size(2) maxima_size(3)]);

nanMap = ~isnan(maxima);

maxima = maxima(nanMap).';

coords.r = coords.r(nanMap).';
coords.c = coords.c(nanMap).';
coords.m = coords.m(nanMap).';

if(nargout > 2)
    response = permute(response,[3 4 1 2]);
    response = repmat(response,[1 maxima_size(1) 1 1]);
    response = response(:,nanMap);
end
if(exist('K','var'))
    coords.K = shiftdim(K,-1);
    coords.K = repmat(coords.K,[maxima_size(1) 1 1]);
    coords.K = coords.K(nanMap).';
end

end

