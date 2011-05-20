classdef WindowSamplingProcess < ImageAnalysisProcess
%Process     
%     
% Hunter Elliott
% 7/2010    
%

    methods (Access = public)
        
        function obj = WindowSamplingProcess(owner,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                
                
                super_args{1} = owner;

                super_args{2} = WindowSamplingProcess.getName;
                super_args{3} = @sampleMovieWindows;                               
                
                if nargin < 2 || isempty(funParams)                                       
                    
                    %----Defaults----%      
                    nChan = numel(owner.channels_);
                    funParams.ChannelIndex = 1:nChan;%Default is to sample all channels                    
                    funParams.ProcessIndex = [];%Default is to use raw images
                    funParams.OutputDirectory = ...
                        [owner.outputDirectory_  filesep 'window_sampling'];
                    funParams.BatchMode = false;                                                                                
                                    
                end
                
                super_args{4} = funParams;    
                                
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end   
        
        function samp = loadChannelOutput(obj,iChan)
            
            if nargin < 2 || isempty(iChan)
                error('You must specify a channel number to load window samples for!');
            elseif round(abs(iChan)) ~= iChan || iChan > numel(obj.owner_.channels_)
                error('The channel number must be a positive integer less than or equal to the number of channels in the movie!')
            end
                        
            tmp = load(obj.outFilePaths_{iChan});
            fNames = fieldnames(tmp);
            if numel(fNames) ~= 1
                error('Invalid window sample file !');
            end
            samp = tmp.(fNames{1});
            
            
        end
                
    end
    methods (Static)
        function name =getName()
            name = 'Window Sampling';
        end
    end 
    
end
       
    