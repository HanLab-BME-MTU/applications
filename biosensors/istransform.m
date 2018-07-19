function isIt = istransform(xFormIn)
%ISTRANSFORM checks if the input variable is an image toolbox image transform
%
% isIt = istransform(xFormIn)
% 
% 
% Input:
%   
%   xFormIn - The variable to be tested.
% 
% Output:
% 
%  iIt - True if the input xFormIn was an Image Processing Toolbox
%  transform, and false otherwise.
% 
% Hunter Elliott
%

if isfield(xFormIn,'ndims_in') && isfield(xFormIn,'ndims_out') && isfield(xFormIn,'forward_fcn') && isfield(xFormIn,'inverse_fcn')...
        && isfield(xFormIn,'tdata')
    
    isIt = true;
    
else
    isIt = false;
end