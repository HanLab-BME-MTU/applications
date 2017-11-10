function [ extended , extended_raw] = extendUsingOrientationField( bw, pts, g, maxPixels, connectivity , mask)
%extendUsingOrientationField Extend endpts using the orientation field from
%a steerable filter
%
% INPUT
%
% bw - logical image containing existing skeleton... usually nms
% pts - endpoints either as a binary image the same size as bw or a two
% column vector indiciating the column and row of endpoints
%   Alternatively a vector of linear indices for bw
%   Alternatively a logical vector
% g - structure containing steerable filter output
%   g.response - filter response, double 2D matrix
%   g.nms - non-maximal suppression, double 2D matrix
%   g.theta - best orientation, double 2D matrix
%   g.a - angular response, double 3D matrix, :dimensions[size(bw) nAngles]
%   g.l - logical matrix size of g.a indicating where local maxima are
% mask (optional) - logical matrix, size of bw, indicating permissible
% paths
%
% OUTPUT
%
% extended - extended image, skeletonized with spurs cut off
% extended_raw - extended image, without further processing
%
% See also extendVectoruntilConnected
%
% @author Mark Kittisopikul
% @date April 6th, 2015

if(nargin < 6)
    mask = true(size(bw));
end

if(size(pts,2) == 2)
    c = pts(:,1);
    r = pts(:,2);
elseif(isvector(pts))
    if(islogical(pts))
        pts = find(pts);
    end    
    [r,c] = ind2sub(size(bw),pts);
else
    [r,c] = find(pts);
end

if(~isfield(g,'l'))
    g.l = findRobustPeaksOnCircle(g.a);
end

la = g.l.*(g.a-repmat(min(g.a,[],3),1,1,size(g.a,3)));
[M,MI] = sort(la,3);
N = size(M,3);

% round to the nearest to 45 degree angle
% thetas = round((MI - 1)*4/N)*pi/4;
% or not

thetas = (MI-1)*pi/N;
% set angles to NaN if there is no or negative response
thetas(M <= 0) = NaN;

% order so that we can index by endpts last
thetas = shiftdim(thetas,2);

% use local orientation to determine initial direction
local = lamins.functions.getEndPointVector(bw,'local');
% normalize local orientation
% local = local ./ repmat(hypot(local(:,1),local(:,2)),1,2);
initTheta = atan2(local(:,2),local(:,1));


assert(~isscalar(connectivity) || connectivity == 4 || connectivity == 8,'Connectivity must either be 4, 8 or a filter');

% setup donut filter to find number of neighbors each pixel has in 8 space
neighborsFilter = ones(3);
neighborsFilter(2,2) = 0;

% extended starts with bw as a double
extended = double(bw);

% boolean array indicating whether each pt should be extended still
keepExtending = true(size(r));
extend_idx = sub2ind(size(bw),r,c);

% r_i = r;
% c_i = c;

nAngles = 3;
considerNormalDirection = true;
thetas = thetas(end-nAngles+1:end,:);
g.response(~mask) = NaN;

for i=1:maxPixels
%     keepExtending(vvr(:,i) < 1) = false;
%     keepExtending(uuc(:,i) < 1) = false;
%     keepExtending(vvr(:,i) > 1024) = false;
%     keepExtending(uuc(:,i) > 1024) = false;
%     extend_idx = sub2ind(size(bw),vvr(keepExtending,i),uuc(keepExtending,i));

    % try the best two angles and their 180 degree rotations
    anglesOfInterest = thetas(:,extend_idx(keepExtending));
    % consider the angle in either the positive or negative direction
    candidateAngles = [ anglesOfInterest+pi/2; 
                        anglesOfInterest-pi/2];
    if(considerNormalDirection)
        % perhaps the normal direction or it's negative is also a candidate
        candidateAngles = [ candidateAngles;
                            anglesOfInterest;
                            anglesOfInterest-pi];
    end
    candidateAngles = candidateAngles';
    deltaAngle = candidateAngles - repmat(initTheta(keepExtending),1,size(candidateAngles,2));
    % find the corresponding response of all the candidate angle moves
    % consider interpolating rather than rounding
    idx = sub2ind(size(g.response), ...
        round(cos(candidateAngles)+repmat(c(keepExtending),1,size(candidateAngles,2))), ...
        round(sin(candidateAngles)+repmat(r(keepExtending),1,size(candidateAngles,2))) ...
        );
    candidateResponse = NaN(size(candidateAngles));
    
    candidateResponse(~isnan(idx)) = g.response(idx(~isnan(idx)));
    % get smallest angle between angles
    deltaAngle = atan2(sin(deltaAngle),cos(deltaAngle));
    [~,minDeltaAngleIdx] = nanmin(abs(deltaAngle)./candidateResponse,[],2);
%     [~,minDeltaAngleIdx] = nanmin(abs(deltaAngle),[],2);
    % select the best angle with the lowest score
    angle = candidateAngles(sub2ind(size(candidateAngles),(1:size(candidateAngles,1))',minDeltaAngleIdx));
    initTheta(keepExtending) = angle;
    r(keepExtending) = r(keepExtending) + sin(angle);
    c(keepExtending) = c(keepExtending) + cos(angle);
    
    r_i = round(r);
    c_i = round(c);
    
    % stop extending at the image boundary
    keepExtending(r_i < 1) = false;
    keepExtending(c_i < 1) = false;
    keepExtending(r_i > size(bw,1)) = false;
    keepExtending(c_i > size(bw,2)) = false;
    % stop extending if NaN
    keepExtending(isnan(r_i)) = false;
    
    extend_idx = sub2ind(size(bw),r_i,c_i);

    extended(extend_idx(keepExtending)) = 1;
    neighbors = imfilter(extended,neighborsFilter);
    extend_neighbors = zeros(size(keepExtending));
    extend_neighbors(keepExtending) = neighbors(extend_idx(keepExtending));
    keepExtending = keepExtending & extend_neighbors < 3;
    
    
end

if(nargout > 1)
    extended_raw = extended;
end

extended = bwmorph(extended,'skel',Inf);
extended = bwmorph(extended,'spur',Inf);

%figure;
%imshowpair(extended,bw);

%figure;
%regions = bwlabel(~extended,4);
%regions(regions == 1) = 0;
%showPropMatrix(regions);

end

