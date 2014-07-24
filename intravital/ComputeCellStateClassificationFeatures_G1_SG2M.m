function [featureStruct, varargout] = ComputeCellStateClassificationFeatures_G1_SG2M( imageData, imValidROIMask, imCellMask, spacing, varargin )

    [featureStruct, featureCompParameters] = ComputeCellStateClassificationFeatures( imageData, imValidROIMask, imCellMask, spacing, varargin{:} );
    
    if nargout > 1
        varargout{1} = featureCompParameters;
    end
    
end