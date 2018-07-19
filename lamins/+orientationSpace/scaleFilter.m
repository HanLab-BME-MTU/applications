function F = scaleFilter(scale,angle,fuzz)

coords = orientationSpace.getFrequencySpaceCoordinates(201);
if(nargin < 2)
    angle = 0;
end
scale = shiftdim(scale(:),-2);
angle = shiftdim(angle(:),-3);

scale = 1/2/pi./scale;
% coords.x = coords.f.*cos(coords.theta+angle);
% coords.y = coords.f.*sin(coords.theta+angle);
coords.x = bsxfun(@times,coords.f,cos(bsxfun(@plus,coords.theta,angle)));
coords.y = bsxfun(@times,coords.f,sin(bsxfun(@plus,coords.theta,angle)));
% F = exp(-coords.x.^2./2./(scale).^2).*abs(coords.x).*exp(-coords.y.^2/2/(1/2/pi/2).^2)/scale;

if(nargin < 3 || fuzz)
    F = abs(coords.x).*exp(-coords.y.^2/2/(1/2/pi/2));
else
    F = abs(coords.x);
end

F = bsxfun(@rdivide,F,scale);
F = exp(-bsxfun(@rdivide,coords.x.^2/2,scale.^2)).*F;

end