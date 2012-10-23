classdef CorrEBtoKinDynamicsProcess < Process
    % A class for correlating kEB signal with kinetochore dynamics
    %
    % Khuloud Jaqaman, October 2012
    
    methods (Access = public)
        
        function obj = CorrEBtoKinDynamicsProcess(owner, varargin)
            
            if nargin == 0
                super_args = {};
                funParams = [];
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) (isa(x,'MovieData')||isa(x,'MovieList')));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                super_args{1} = owner;
                super_args{2} = CorrEBtoKinDynamicsProcess.getName;
            end
            
            obj = obj@Process(super_args{:});
            
            obj.funName_ = @corrMoviesEBtoKinDynamics;
            if isempty(funParams)  % Default funParams
                funParams = CorrEBtoKinDynamicsProcess.getDefaultParams(owner,outputDir);
            end
            obj.funParams_ = funParams;
            
        end
                
    end
    
    methods (Static)
        
        function name = getName()
            name = 'kEB - kin dynamics correlation';
        end
        
        function funParams = getDefaultParams(owner,varargin)
            
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) (isa(x,'MovieData')||isa(x,'MovieList')));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'EBtoKinDynamicsCorr'];
            funParams.minDisp = 0.5;
            
        end
        
    end
    
end
