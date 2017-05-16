classdef StrainEnergyCalculationProcess < DataProcessingProcess
    % Concrete process for calculating a force field
    %
    % Sebastien Besson, Aug 2011
    properties (SetAccess = protected)  
        tMapLimits_
        dELimits_
        distBeadMapLimits_
    end
    
    methods
        function obj = StrainEnergyCalculationProcess(owner,varargin)
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = StrainEnergyCalculationProcess.getName;
                super_args{3} = @calculateMovieStrainEnergy;
                if isempty(funParams)
                    funParams=StrainEnergyCalculationProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            
            obj = obj@DataProcessingProcess(super_args{:});
            
        end
    end
%         function status = checkChannelOutput(obj,varargin)
%             
%             status = logical(exist(obj.outFilePaths_{1},'file'));
%             
%         end
        
    methods (Static)
        function name =getName()
            name = 'Strain Energy Calculation';
        end
        function h = GUI()
            h= @strainEnergyCalculationProcessGUI;
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'strainEnergy'];
            funParams.exportCSV=true;
            funParams.performForceBlobAnalysis=true;
            funParams.useFOV=true;
            funParams.useCellMask=true;
        end
%         function units = getUnits(varargin)
%             units = 'Traction (Pa)';
%         end
    end
end