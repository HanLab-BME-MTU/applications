classdef LocalSBRThresholdingModule < ForegroundSegmentationModule
% LocalSBRThresholdingModule - Module for segmentation
% of foreground by local signal to background ratio thresholding
%
%   The foreground is segmented by applying a cut off to the local 
%   signal to background ratio estimate. The local bacground signal is
%   estimated using grayscale morphological opening. This was found 
%   to work well for segmenting nuclei foreground in most of the 
%   intravital microscopy data and 3D cell culture data.
%
    
    properties(SetAccess = private)
        flagDebugMode
        flagParallelize
    end

    methods

        function obj = LocalSBRThresholdingModule(varargin)

            p = inputParser;    
            p.addParamValue('debugMode', false, @(x) (isscalar(x) && islogical(x)));    
            p.addParamValue('flagParallelize', false, @(x) (isscalar(x) && islogical(x)));
            p.parse(varargin{:});

            obj.flagParallelize = p.Results.flagParallelize;
            obj.flagDebugMode = p.Results.debugMode;

            % default parameters
            obj.parameters.maxObjectRadius = 20;
            obj.parameters.sbrCutoff = 2.5;            
            obj.parameters.kernelDimensions = 2;
            obj.parameters.minObjectRadius = 0;
            obj.parameters.downsamplingFactor = 1.0;

        end

        function [strName] = getName(obj)
            strName = 'Local Signal To Background Ratio Thresholding';
        end

        function [imForegroundMask] = segmentImage(obj, imInput, spacing)            

            if ~exist('spacing', 'var')
                spacing = ones(1, ndims(imInput));
            end

            imForegroundMask = thresholdSBR(imInput, obj.parameters.maxObjectRadius, obj.parameters.sbrCutoff, ...
                                            'spacing', spacing, ... 
                                            'minObjectRadius', obj.parameters.minObjectRadius, ...
                                            'kernelDimensions', obj.parameters.kernelDimensions, ...
                                            'downsamplingFactor', obj.parameters.downsamplingFactor, ...
                                            'flagDebugMode', obj.flagDebugMode);

        end

        function uiSetParameters(obj)

            kDimList = [obj.parameters.kernelDimensions, setdiff([2,3], obj.parameters.kernelDimensions)];

            % parameter arguments         
            pArgs = {{'Max Object Radius (um)', 'maxObjectRadius'}, obj.parameters.maxObjectRadius, ...
                     {'Signal To Background Ratio Cutoff', 'sbrCutoff'}, obj.parameters.sbrCutoff, ...
                     {'Min Object Radius', 'minObjectRadius'}, obj.parameters.minObjectRadius, ...
                     {'Opening Kernel Dimension', 'kernelDimensions'}, num2cell(kDimList), ...
                     {'Downsampling Factor (0-1)', 'downsamplingFactor'}, obj.parameters.downsamplingFactor
                    };

            strDesc = ['The foreground is segmented by applying a cut off to the local ' ...
                       'signal to background ratio estimate. The local bacground signal is ' ...
                       'estimated using grayscale morphological opening. This was found ' ...
                       'o work well for segmenting nuclei foreground in most of the ' ...
                       'intravital microscopy data and 3D cell culture data.'];
            
            [usel, bsel] = settingsdlg('Description', sprintf('Parameters for %s\n\n%s', obj.getName(), strDesc), ...
                                        pArgs{:});
              
            if ~strcmpi(bsel, 'ok')
                return;
            end
            
            obj.parameters = usel;

        end

    end        
    
end