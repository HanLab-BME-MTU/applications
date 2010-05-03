function [params, Im] = focalAdhesionDetector(I, sigmaPSF, minSize)
% [params, Im] = focalAdhesionDetector(I, sigmaPSF, minSize)
% THIS FILE IS NOT USED ANYMORE

% Make sure I is type double
I = double(I);

% Get initial segment parameters
[params, I] = getInitialSegmentParams(I,sigmaPSF,minSize);

% % Define bounds
% lb = zeros(size(params));
% % min value of center coordinates
% lb(:,1:2) = 1;
% % min value of segment amplitude
% lb(:,3) = 0;
% % min length
% lb(:,4) = 0;
% % min orientation value
% lb(:,5) = -inf;
% 
% ub = zeros(size(params));
% % max x-coordinate
% ub(:,1) = size(I,2);
% % max y-coordinate
% ub(:,2) = size(I,1);
% % max value of signal intensity
% ub(:,3) = +inf;
% % max length
% ub(:,4) = +inf;
% % max orientation
% ub(:,5) = +inf;
% 
% % Set options for lsqnonlin
% options = optimset('Jacobian', 'on', 'MaxFunEvals', 1e3, 'MaxIter', 1e3, ...
%     'Display', 'off', 'TolX', 1e-4, 'Tolfun', 1e-4);
% 
% fun = @(x) subResSegment2DFit(x, I, sigmaPSF);
% [params, ~,residual] = lsqnonlin(fun, params, lb, ub, options);
% Im = I - reshape(residual, size(I));
Im = subResSegment2DImageModel(params,sigmaPSF,size(I));