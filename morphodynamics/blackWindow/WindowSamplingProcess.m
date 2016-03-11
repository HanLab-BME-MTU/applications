classdef WindowSamplingProcess < ImageSamplingProcess
    %Process
    %
    % Hunter Elliott
    % 7/2010
    %
    
    methods (Access = public)
        
        function obj = WindowSamplingProcess(owner,varargin)
            
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
                super_args{2} = WindowSamplingProcess.getName;
                super_args{3} = @sampleMovieWindows;
                if isempty(funParams)
                    funParams=WindowSamplingProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
                
            end
            
            obj = obj@ImageSamplingProcess(super_args{:});
        end
        

    end
    methods (Static)
        function name =getName()
            name = 'Window Sampling';
        end
        function name= GUI()
            name =@windowSamplingProcessGUI;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1:numel(owner.channels_);%Default is to sample all channels
            funParams.ProcessIndex = [];%Default is to use raw images
            funParams.SegProcessIndex = [];%Default is to use masks which were used in windowing.
            funParams.MaskChannelIndex = [];%Default is to use channel which was used for windowing.
            funParams.OutputName = '';%Default is to use raw images
            funParams.OutputDirectory = [outputDir  filesep 'window_sampling'];
            funParams.BatchMode = false;
        end
        function samplableInput = getSamplableInput()
            % List process output that can be sampled
            processNames = horzcat('Raw images','DoubleProcessingProcess',...
                repmat({'KineticAnalysisProcess'},1,3),'FlowAnalysisProcess');
            samplableOutput = {'','','netMap','polyMap','depolyMap','speedMap'};
            samplableInput=cell2struct(vertcat(processNames,samplableOutput),...
                {'processName','samplableOutput'});  
        end
        
        
    end
end