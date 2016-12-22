function F = scaleFilter(scale,angle)

coords = orientationSpace.getFrequencySpaceCoordinates(201);
if(nargin < 2)
    angle = 0;
end
scale = 1/2/pi/scale;
coords.x = coords.f.*cos(coords.theta+angle);
coords.y = coords.f.*sin(coords.theta+angle);
F = exp(-coords.x.^2./2./(scale).^2).*abs(coords.x).*exp(-coords.y.^2/2/(1/2/pi/2).^2)/scale;

end