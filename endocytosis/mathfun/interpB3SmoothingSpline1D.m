%[fi, fi_dx, fi_d2x] = interpB3SmoothingSpline1D(xi, f, varargin)
%
% Inputs:
%           xi : interpolation coordinates
%            f : function to interpolate
%
% Optional inputs:
%            n : degree of the spline, 1-3. Default: 3
%     boundary : boundary conditions: 'symmetric' (default) or 'periodic'
%
% Outputs:
%           fi : interpolated values at xi
%        fi_dx : derivative of f at xi 
%       fi_d2x : 2nd derivative of f at xi

% Francois Aguet, 2009 (Last modified: 10/05/2011)

function [fi, fi_dx, fi_d2x] = interpB3SmoothingSpline1D(xi, f, varargin)

ip = inputParser;
ip.addRequired('xi', @isvector);
ip.addRequired('f', @isvector);
ip.addOptional('lambda', 0, @(x) isscalar(x) && x>=0);
ip.addOptional('boundary', 'symmetric', @(x) any(strcmpi(x, {'symmetric', 'periodic'})));
ip.parse(xi, f, varargin{:});

c = computeBSplineCoefficients(f, ip.Results.lambda, 'Boundary', ip.Results.boundary);

fi = arrayfun(@(i) interpBSplineValue(i, c, 3, ip.Results.boundary), xi);
if nargout>1
    fi_dx = arrayfun(@(i) interpBSplineValue(i+0.5, c, 2, ip.Results.boundary) - ...
        interpBSplineValue(i-0.5, c, 2, ip.Results.boundary), xi);
end
if nargout>2
    fi_d2x = arrayfun(@(i) interpBSplineValue(i+1, c, 1, ip.Results.boundary) - ...
        2*interpBSplineValue(i, c, 1, ip.Results.boundary) + ...
        interpBSplineValue(i-1, c, 1, ip.Results.boundary), xi);
end
