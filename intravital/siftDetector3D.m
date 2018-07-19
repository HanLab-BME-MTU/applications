function points = siftInterestPointDetector(im, varargin)

    p = inputParser;    
    p.addRequired('im', @(x) (isnumeric(x) && ismember( ndims(x), [2,3] )));
    p.addParamValue('spacing', ones( 1, ndims(im) ), @(x) (isnumeric(x) && numel(x) == ndims(im)));
    p.addParamValue('numOctaves', 4, @isscalar);
    p.addParamValue('numScalesPerOctave', 3, @isscalar);
    p.addParamValue('dogResponseCutoff', 0.0, @isscalar);
    p.addParamValue('blobnessCutoff', 5, @isscalar);
    p.parse(im, varargin{:});

    
end