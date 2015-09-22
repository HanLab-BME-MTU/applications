classdef LocalMinimumErrorThresholdingModule2D < ForegroundSegmentationModule
% LocalMinimumErrorThresholdingModule2D Module for segmentation
% of foreground using the local minimum error thresholding algorithm 
%
%   The minimum error thresholding algorithm is applied in a locally
%   adaptive fashion. If a 3D image/volume is provided then the algorithm
%   is applied on each slice separately.
    
    properties
        flagDebugMode
        flagParallelize
    end
    
    methods
        
        function obj = LocalMinimumErrorThresholdingModule2D(varargin)
            
            p = inputParser;    
            p.addParamValue('debugMode', false, @(x) (isscalar(x) && islogical(x)));    
            p.addParamValue('flagParallelize', false, @(x) (isscalar(x) && islogical(x)));
            p.parse(varargin{:});
            
            obj.flagParallelize = p.Results.flagParallelize;
            obj.flagDebugMode = p.Results.debugMode;
            
            % default parameters
            obj.parameters.model = 'poisson';
            obj.parameters.localThresholdWindowRadiusPhysp = 30;
            obj.parameters.localWindowPaceFraction = 1/3;
            obj.parameters.minLocalGlobalThresholdRatio = 0.6;
            obj.parameters.minSliceToStackThresholdRatio = 0.4;
            obj.parameters.numHistogramBins = 256;
            
        end
        
        function [strName] = getName(obj)
            strName = 'Slice By Slice Local Minimum Error Thresholding';
        end
        
        function [imForegroundMask] = segmentImage(obj, imInput, spacing)            
            
            if ~exist('spacing', 'var')
                spacing = ones(1, ndims(imInput));
            end
            
            localWindowRadius = round(obj.parameters.localThresholdWindowRadiusPhysp ./ spacing(1));
            localWindowPace = round(localWindowRadius * obj.parameters.localWindowPaceFraction);
            
            imForegroundMask = segmentCellForegroundUsingLocalMinError( imInput, localWindowRadius, ...
                                                                        'model', obj.parameters.model, ...  
                                                                        'localWindowPace', localWindowPace, ...
                                                                        'minLocalGlobalThresholdRatio', obj.parameters.minLocalGlobalThresholdRatio, ...
                                                                        'minSliceToStackThresholdRatio', obj.parameters.minSliceToStackThresholdRatio, ...
                                                                        'numHistogramBins', obj.parameters.numHistogramBins, ...
                                                                        'debugMode', obj.flagDebugMode, ...
                                                                        'flagParallelize', obj.flagParallelize );
            
        end
        
        function uiSetParameters(obj)
        
            modelList = [obj.parameters.model, ...
                         setdiff({'poisson', 'gaussian'}, obj.parameters.model)];
             
            % parameter arguments         
            pArgs = {{'Foreground/background intensity model', 'model'}, modelList, ...
                     {'Window Radius (um)', 'localThresholdWindowRadiusPhysp'}, obj.parameters.localThresholdWindowRadiusPhysp, ...
                     {'Window Step Fraction', 'localWindowPaceFraction'}, obj.parameters.localWindowPaceFraction, ...
                     {'Minimum Local Global Threshold Ratio (0-1)', 'minLocalGlobalThresholdRatio'}, obj.parameters.minLocalGlobalThresholdRatio, ...
                     {'Minimum Slice To Stack Threshold Ratio (0-1)', 'minSliceToStackThresholdRatio'}, obj.parameters.minSliceToStackThresholdRatio, ...
                     {'No of Histogram Bins', 'numHistogramBins'}, obj.parameters.numHistogramBins, ...
                    };
            
            strDesc = ['The minimum error thresholding algorithm is applied in a locally ' ...
                       'adaptive fashion. On 3D images/volumes the algorithm is applied ' ...
					   'on each slice separately.'];
			
            [usel, bsel] = settingsdlg('Description', sprintf('Parameters for %s\n\n%s', obj.getName(), strDesc), ...
                                        pArgs{:});
              
            if ~strcmpi(bsel, 'ok')
                return;
            end
            
            obj.parameters = usel;
            
        end
        
    end    
    
end

