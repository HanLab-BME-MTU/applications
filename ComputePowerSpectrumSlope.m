function [powerSlope] = ComputePowerSpectrumSlope( im, spacing )

    if ~exist( 'spacing', 'var' )
        spacing = ones(1, ndims(im));
    end
    
    imsize = size(im);
    for i = 1:ndims(im)
       if mod(imsize(i), 2) ~= 0
           imsize(i) = imsize(i) - 1;
       end
    end
    
    freqComp = cell(1, ndims(im) );
    imFreq = cell(1, ndims(im) );
    maxRadVec = zeros(1, ndims(im));
    
    for i = 1:ndims(im)
        freqComp{i} = [0:(imsize(i)/2-1), -imsize(i)/2:-1] / spacing(i);
        maxRadVec(i) = max( freqComp{i} );
    end
    [imFreq{:}] = ndgrid( freqComp{:} );

    imRadFreq = zeros(size(imFreq{1}));
    for i = 1:ndims(im)
        imRadFreq = imRadFreq + imFreq{i}.^2;
    end
    imRadFreq = sqrt(imRadFreq);
    imFFTMag = abs( fftn(im - mean(im(:)), imsize) ) / prod(imsize);  
    
    flagValidFreq = true( size(imRadFreq) );
    flagValidFreq(1,1) = false; % remove dc component
    
    p = polyfit( log(imRadFreq(flagValidFreq)), log(imFFTMag(flagValidFreq)), 1 );
    pRobust = robustfit( log(imRadFreq(flagValidFreq)), log(imFFTMag(flagValidFreq)) );
    
%     figure, plot( log(imRadFreq(flagValidFreq)), log(imFFTMag(flagValidFreq)), '.' )
%     hold on;
%     plot( log(imRadFreq(flagValidFreq)), p(1) * log(imRadFreq(flagValidFreq)) + p(2), 'r' );
%     plot( log(imRadFreq(flagValidFreq)), pRobust(2) * log(imRadFreq(flagValidFreq)) + pRobust(1), 'g' );
%     hold off;
    
    powerSlope = p(1);
    
end