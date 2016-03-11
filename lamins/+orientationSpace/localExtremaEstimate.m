function [ maxima, maxima_r, minima, minima_r ] = localExtremaEstimate( orientationMatrix, scale, doSort, doCompact )
%orientationSpace.localExtremaEstimate Estimate location of local extrema
    import orientationSpace.*;
    s = size(orientationMatrix);
%     M = real(reshape(orientationMatrix,s(1)*s(2),s(3)));
    n = (s(3))/2;
    
    if(nargin < 2)
        scale = 1;
    end
    if(nargin < 3)
        doSort = false;
    end
    if(nargin < 4)
        doCompact = true;
    end
    
    % new
    % estimate zeros in derivative
    deriv = orientationSpace.derivative(real(orientationMatrix),pi/s(3)/scale);
    derivDiff = -diff(deriv(:,:,[1:end 1]),1,3);
    offset = deriv./derivDiff;
    maxima = bsxfun(@plus,offset,shiftdim(0:s(3)*scale-1,-1));
    maxima(offset < 0 | offset > 1) = NaN;
    scaleFactor = pi/s(3)/scale;
    if(nargout > 2)
        minima = maxima;
        minima(derivDiff > 0) = NaN;
        if(doCompact)
            minima = compactNaN(minima);
%         else
%             NaNmap = isnan(minima);
%             minima(minima == 0) = NaN;
%             minima(NaNmap) = 0;
%             minima = sparse(minima(:));
        end
        
        if(doSort)
            [minima, minima_r] = sortExtrema(minima,Inf);
        elseif(nargout > 3)
            minima_r = orientationSpace.interpolate(orientationMatrix,minima);
        end
    end
    maxima(derivDiff < 0 ) = NaN;
    if(doCompact)
        maxima = compactNaN(maxima);
%     else
%         NaNmap = isnan(maxima);
%         maxima(maxima == 0) = NaN;
%         maxima(NaNmap) = 0;
%         maxima = sparse(maxima(:));
    end
    
    if(doSort)
        [maxima , maxima_r] = sortExtrema(maxima,-Inf);
    elseif(nargout > 1)
        maxima_r = orientationSpace.interpolate(orientationMatrix,maxima);
    end
    
    function out =  compactNaN(in)
        notNaNMap = ~isnan(in);
        notNaNCount = sum(notNaNMap,3);
        maxNotNaNCount = max(notNaNCount(:));

        % pad maxima in the 3rd dimension with extra NaNs
        in = cat(3,in,NaN(s(1),s(2),maxNotNaNCount));
        in = shiftdim(in,2);

        % select elements in the padding to pad final selection to
        % maxNotNanCount
        nanPad = bsxfun(@lt,notNaNCount,shiftdim(1:maxNotNaNCount,-1));
        notNaNMap = cat(3,notNaNMap,nanPad);
        notNaNMap = shiftdim(notNaNMap,2);

        in = in(notNaNMap);
        in = reshape(in,maxNotNaNCount,s(1),s(2));
        out = shiftdim(in,1);
        out = out*scaleFactor;
        out = wraparoundN(out,-pi/2,pi/2);
    end
    function [out,sr] = sortExtrema(in,nanValue)
        r = orientationSpace.interpolate(orientationMatrix,in);
        r(~isfinite(in)) = nanValue;
        [sr,~,out] = cosort(r,3,'descend',in);
    end

end

