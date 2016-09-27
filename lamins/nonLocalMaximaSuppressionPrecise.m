function [ nlms ] = nonLocalMaximaSuppressionPrecise( rotationResponse, theta , suppressionValue, interpMethod)
%nonLocalMaximaSuppression Suppress pixels which are not local maxima in
%both orientation and filter response
%
% INPUT
% rotationResponse: matrix Y by X by R, R = # of rotation angles
%             (4th output of steerableDetector or steerableVanGinkelFilter)
% theta: (optional) rotation angles. Default: (0:R-1)*pi/R
% suppressionValue: (optional) value to give suppressed pixels (def: 0)
%
% OUTPUT
% nlms: rotationResponse with non-local maxima set to suppressedValue

% Mark Kittisopikul, 2015
% UT Southwestern

% TODO: full-backwards compatability with nonMaximumSuppression?

nO = size(theta,3);

if(nargin < 2 || isempty(theta))
    % default value is the rotation planes correspond to angles that
    % equally divide pi
    theta = 0:nO-1;
    theta = theta*pi/nO;
    %TODO: allow for theta of different sizes
end
if(nargin < 3 || isempty(suppressionValue))
    suppressionValue = 0;
end
if(nargin < 4 || isempty(interpMethod))
    interpMethod = 'cubic';
end
% if(nargin < 4)
%     distance = 1;
% end

% nlms = theta;
rotationResponseSize = size(rotationResponse);
ny = rotationResponseSize(1);
nx = rotationResponseSize(2);
nAngles = rotationResponseSize(3);
period = nAngles*3;

% See Boyd, Chebyshev and Fourier Spectral Methods, Second Edition
% (Revised), Dover, 2001. Toronto. ISBN 0-486-41183-4, page 198
rotationResponse = interpft(rotationResponse,period,3);
rotationResponse = padarrayXT(rotationResponse,[1 1 0] ,'symmetric');
rotationResponse = padarrayXT(rotationResponse,[0 0 1] ,'circular');

% Map angles in radians (0:2*pi) to angle index (1:nAngles)
% angleIdx = theta/(2*pi)*(nAngles*3-1)+1;

% Map angles in from [-pi/2:pi/2) to (0:pi]
angleIdx = theta;
angleIdx(angleIdx < 0) = angleIdx(angleIdx < 0) + pi; % now (0:pi)
% Map angles in radians (0:pi) to angle index (2:nAngles*3+1)
angleIdx = angleIdx/pi*(period)+2;
% Negative values should be advanced by one period
% angleIdx(angleIdx < 0) = angleIdx(angleIdx < 0)+period;

% Offset by 1 due to padding
[x,y] = meshgrid(2:nx+1,2:ny+1);

x_offset = cos(theta);
y_offset = sin(theta);

Xplus = bsxfun(@plus,x,x_offset);
Yplus = bsxfun(@plus,y,y_offset);

Xminus = bsxfun(@minus,x,x_offset);
Yminus = bsxfun(@minus,y,y_offset);

x = cat(4,Xminus,repmat(x,[1 1 nO]),Xplus);
y = cat(4,Yminus,repmat(y,[1 1 nO]),Yplus);
angleIdx = repmat(angleIdx,[1 1 1 3]);

clear Xplus Yplus Xminus Yminus x_offset y_offset theta;

A = interp3(rotationResponse,x,y,angleIdx,interpMethod,0);

nlms = A(:,:,:,2);
nlms(nlms < A(:,:,:,1) | nlms < A(:,:,:,3)) = suppressionValue;

% TODO: optimize later
% for o = 1:nO
% %     nlms(:,:,o) = nonMaximumSuppression(rotationResponse(:,:,o),theta(o));
% 
% %     res = padarrayXT(rrotationResponse(:,:,o), [1 1 0], 'symmetric');
% 
%     x_offset = cos(theta(:,:,o));
%     y_offset = sin(theta(:,:,o));
% 
%     % +1 interp
%     A1 = interp3(rotationResponse, x+x_offset, y+y_offset,angleIdx(:,:,o),interpMethod,0);
% 
%     % -1 interp
%     A2 = interp3(rotationResponse, x-x_offset, y-y_offset,angleIdx(:,:,o),interpMethod,0);
% 
%     % TODO: We only need to interpolate where not NaN
%     temp = interp3(rotationResponse, x, y, angleIdx(:,:,o),interpMethod,0);
%     temp(temp < A1 | temp < A2) = suppressionValue;
%     nlms(:,:,o) = temp;
% 
% %     res(res<A1 | res<A2) = 0;
% end

%% since the input theta already should represent the local maxima in
%% in orientation, we no longer need to suppress in the orientation dimension
% suppression response if less than the previous or next orientation
% response
% nlms(nlms < nlms(:,:,[2:end 1]) | nlms < nlms(:,:,[end 1:end-1])) = suppressionValue;
% nlms(rotationResponse < rotationResponse(:,:,[2:end 1]) | rotationResponse < rotationResponse(:,:,[end 1:end-1])) = suppressionValue;
% what if equal on either side or both sides?

end % end of function
