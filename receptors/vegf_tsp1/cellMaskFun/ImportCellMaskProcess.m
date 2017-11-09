classdef ImportCellMaskProcess < DataProcessingProcess
    %Class definition for cell mask import
        
    methods(Access = public)   
        
        function obj = ImportCellMaskProcess(owner,varargin)
            
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
                super_args{2} = ImportCellMaskProcess.getName;
                super_args{3} = @importCellMaskOO;                               
                if isempty(funParams)                                       
                    funParams = ImportCellMaskProcess.getDefaultParams(owner,outputDir);                               
                end
                super_args{4} = funParams;                    
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
        end        
        
    end
    
    methods(Static)
        
        function name =getName()
            name = 'Cell Mask Import';
        end

        function funParams = getDefaultParams(owner,varargin)
            
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:});
            outputDir=ip.Results.outputDir;
            
            % Define default process parameters
            funParams.OutputDirectory = [outputDir  filesep 'ImportedCellMask'];
            funParams.fileName = cell(size(owner.channels_));
            funParams.filePath = cell(size(owner.channels_));
            funParams.ChannelIndex = 1;
            funParams.askUser = true;
        end
        
    end
    
end
