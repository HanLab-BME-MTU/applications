classdef LoGResponseThresholdingModule < ForegroundSegmentationModule
% LoGResponseThresholdingModule Module for segmentation
% of foreground by blobness (Laplacian of Gaussian Response)
% thresholding.
%
%   The foreground is segmented by applying a cut off to the response
%   of the Laplacian of Gaussian Filter which is a measure of blobness. 
%   This was found  to work well for segmenting nuclei foreground in 
%   zebrafish and drosophila data where the nuclei are uniform in size
%   and more or less round in shape.
%
    
    properties
        flagDebugMode
        flagParallelize
    end
    
    methods
        
        function obj = LoGResponseThresholdingModule(varargin)
            
            p = inputParser;    
            p.addParamValue('debugMode', false, @(x) (isscalar(x) && islogical(x)));    
            p.addParamValue('flagParallelize', false, @(x) (isscalar(x) && islogical(x)));
            p.parse(varargin{:});
            
            obj.flagParallelize = p.Results.flagParallelize;
            obj.flagDebugMode = p.Results.debugMode;
            
            % default parameters
            obj.parameters.blobDiameter = 20;
            obj.parameters.choice_of_threshold = 'MinimumError';            
            
        end
        
        function [strName] = getName(obj)
            strName = 'Blobness (Laplacian of Gaussian Response) Thresholding';
        end
        
        function [imForegroundMask] = segmentImage(obj, imInput, spacing)            
            
            if ~exist('spacing', 'var')
                spacing = ones(1, ndims(imInput));
            end
            
            imForegroundMask = thresholdBlobness(imInput, obj.parameters.blobDiameter, ...
                                                'spacing', spacing, ... 
                                                'choice_of_threshold', obj.parameters.choice_of_threshold, ...
                                                'flagDebugMode', obj.flagDebugMode);

        end
        
        function uiSetParameters(obj)
        
            thAlgoList = [obj.parameters.choice_of_threshold, ...
                          setdiff({'Otsu', 'Rosin', 'MinimumError'}, obj.parameters.choice_of_threshold)];
            
            % parameter arguments         
            pArgs = {{'Approx Blob Diameter (um)', 'blobDiameter'}, obj.parameters.blobDiameter, ...
                     {'Thresholding Method', 'choice_of_threshold'}, thAlgoList, ...
                    };
            
            strDesc = ['The foreground is segmented by applying ' ...
                       'a cut off to the response of the Laplacian ' ...
                       'of Gaussian Filter which is a measure of ' ...
                       'blobness. This was found  to work well for ' ...
                       'segmenting nuclei foreground in zebrafish ' ...
                       'and drosophila data where the nuclei are ' ...
                       'uniform in size and more or less round in ' ...
                       'shape.'];
            
            [usel, bsel] = settingsdlg('Description', sprintf('Parameters for %s\n\n%s', obj.getName(), strDesc), ...
                                        pArgs{:});
              
            if ~strcmpi(bsel, 'ok')
                return;
            end
            
            obj.parameters = usel;
            
        end
        
    end    
    
end

