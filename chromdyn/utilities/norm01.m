function [in,add,mult] = norm01(in,quantile)
%NORM01 norms an array to 0..1
%
% SYNOPSIS: out = norm01(in)
%
% INPUT in: any array. uint8 and uint16 are supported.
%       quantile: [min% max%] will set the lowest min% pixels to 0 and the
%           highest max% pixels to 1. Default: [0 100]
%
% OUTPUT in: same array, but with minimum 0 and maximum 1
%        add, mult: If you need to backtransform the image:
%           originalImage = (outputImage*mult)+add
%
% REMARKS
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 22-Sep-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(in) || ~any(isfinite(in(:)))  || all(in(:)==0)
    return
end

arrayClass = class(in);
% check for uint8, uint16. Everything else will become double. This could
% cause memory issues, so one should try to convert logical to uint8 before
% calling this function!
switch arrayClass
    case 'uint8'
        maxMax = 255;
    case 'uint16'
        maxMax = 65535;
    otherwise
        maxMax = 1;
        in = double(in); 
end



    if nargin < 2 || isempty(quantile)
        minMovie = nanmin(in(:));
        maxMovie = (nanmax(in(:))-minMovie)/maxMax;

        in = in - minMovie;
        in = in/(maxMovie);
        if nargout > 1
            mult = maxMovie;
            add = minMovie;
        end
    else
        if length(quantile) == 1
            quantile = [quantile 100-quantile];
        end
        minMovie = prctile(in(:),min(quantile));
        maxMovie = (prctile(in(:),max(quantile))-minMovie);

        in = max(in,minMovie);
        in = in - minMovie;
        in = min(in,maxMovie);
        in = in/(double(maxMovie)/double(maxMax));
        if nargout > 1
            mult = double(maxMovie)*double(maxMax);
            add = minMovie;
        end
    end

