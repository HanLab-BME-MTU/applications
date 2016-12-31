function [ scaleFilter ] = steerableScaleFilter( scale, orientation, sz )
%steerableScaleFilter A Steerable Filter to detect scale 

if(nargin < 1)
    scale = 1/2/pi/2;
end
if(nargin < 2)
    orientation = 0;
end
if(nargin < 3)
    sz = 201;
end

angularOrder = 8;
nAngles = 8;
angles = 0:nAngles-1;
angles = angles/nAngles/2*pi;
%angles = (0:1/nAngles/2:(1/2-1/nAngles/2))*pi;
% angles = 0:1/nAngles/2*pi:acos(scale);
% angles = linspace(0,acos(scale),nAngles);
% angles = angles(1:end-1);
scales = scale./cos(angles);

% nAngles= length(angles);

F = OrientationSpaceFilter(scales,scales,angularOrder,'none');
F = F(logical(eye(nAngles)));
F.setupFilter(sz);

scaleFilter = F(1).getFilterAtAngle(angles(1)+orientation);
for ii=2:nAngles
    scaleFilter = scaleFilter + F(ii).getFilterAtAngle( angles(ii)+orientation);
    scaleFilter = scaleFilter + F(ii).getFilterAtAngle(-angles(ii)+orientation);
end;

end

