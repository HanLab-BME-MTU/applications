function [featureStruct, varargout] = ComputeCellStateClassificationFeatures_G1SG2_M( imageData, imValidROIMask, imCellMask, spacing, varargin )

    [featureStruct, featureCompParameters] = ComputeCellStateClassificationFeatures( imageData, imValidROIMask, imCellMask, spacing, varargin{:} );
    
    if nargout > 1
        varargout{1} = featureCompParameters;
    end
    
end
