function [ scaleFilter ] = steerableScaleFilter( scale, sz )
%steerableScaleFilter A Steerable Filter to detect scale 

if(nargin < 1)
    scale = 1/2/pi/2;
end
if(nargin < 2)
    sz = 201;
end

if(~isscalar(scale))
    scaleFilter(numel(scale)) = OrientationSpaceFilter;
    for jj=1:numel(scale)
        scaleFilter(jj) = orientationSpace.steerableScaleFilter(scale(jj),sz);
    end
    scaleFilter = reshape(scaleFilter,size(scale));
    return;
end

angularOrder = 8;
nAngles = angularOrder*2+1;
angles = 0:angularOrder;
angles = angles/nAngles*pi;
%angles = (0:1/nAngles/2:(1/2-1/nAngles/2))*pi;
% angles = 0:1/nAngles/2*pi:acos(scale);
% angles = linspace(0,acos(scale),nAngles);
% angles = angles(1:end-1);
scales = abs(scale./cos(angles));
disp(1/2/pi/max(scales));

% nAngles= length(angles);

F = OrientationSpaceFilter(scales,scales,angularOrder,'none');
F = F(logical(eye(length(F))));
F.setupFilter(sz);
F(angularOrder+2:nAngles) = F(end:-1:2);

scaleFilter = zeros(size(F(1).F));

for ii=1:nAngles
    scaleFilter = scaleFilter + circshift(F(ii).F,ii-1,3);
end

F(1).F = scaleFilter;
scaleFilter = F(1);

end

