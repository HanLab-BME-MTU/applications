function [ s, rho, I ] = synthesizeIntersection( lineAngle, noise)
%synthesizeIntersection Gets steerable filter information about an
%intersection

if(nargin < 2)
    noise = 1e-3;
end
    import intersections.*;

    I = drawTwoLines(lineAngle);
    I = imfilter(I,fspecial('gaussian',10,2));
    I = imnoise(I,'gaussian',0,noise);
    [s.res, s.theta, s.nms, s.a] = steerableDetector(double(I),4,2);
    rho = s.a( 50, 50 , : );
end

