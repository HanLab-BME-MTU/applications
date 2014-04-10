function [ intensityRange ] = ComputeImageDynamicRange( im, cover_percent )

    [p,x] = hist( double(im), 255 );   
    p = p / sum(p);
    
    min_xlow = [];
    min_xhigh = [];
    min_xdiff = [];
    min_xcover = [];
    
    for i = 1:numel(x)
        for j = i+1:numel(x)
    
            if sum( p(i:j) ) < 0.01 * cover_percent
                continue;
            end
            
            if isempty(min_xdiff) || (x(j) - x(i)) < min_xdiff
                min_xlow = x(i);
                min_xhigh = x(j);
                min_xdiff = x(j) - x(i);              
            end
        end
    end
    
    w = 0.5 * (x(2) - x(1));
    intensityRange = [min_xlow-w, min_xhigh+w];
    
end