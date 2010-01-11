function isIt = istransform(xFormIn)

%Checks if the input variable is an image toolbox image transform

if isfield(xFormIn,'ndims_in') && isfield(xFormIn,'ndims_out') && isfield(xFormIn,'forward_fcn') && isfield(xFormIn,'inverse_fcn')...
        && isfield(xFormIn,'tdata')
    
    isIt = true;
    
else
    isIt = false;
end