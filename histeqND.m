function [imAdjusted] = histeqND( im, numHistogramBins )

    if ~exist( numHistogramBins )
        numHistogramBins = max(im(:)) - min(im(:));        
    end
    
    hist


end