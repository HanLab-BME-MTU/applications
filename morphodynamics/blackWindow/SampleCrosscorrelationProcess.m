classdef SampleCrosscorrelationProcess < DataProcessingProcess
%Process class definition for calculateWindowSampleAutocorrelation.m
%     
% Hunter Elliott
% 8/2011
%
    methods (Access = public)
        
        function obj = SampleCrosscorrelationProcess(owner,outputDir,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                
                
                super_args{1} = owner;

                super_args{2} = SampleCrosscorrelationProcess.getName;
                super_args{3} = @calculateWindowSampleCrosscorrelation;                               
                
                if nargin < 3 || isempty(funParams,outputDir)
                    if nargin < 2, outputDir = owner.outputDirectory_; end
                    funParams = SampleCrosscorrelationProcess.getDefaultParams(owner,outputDir);
                    
                end
                
                super_args{4} = funParams;    
                                
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
        end   
        
        function samp = loadChannelOutput(obj,iChan)                        
            
        end               
           
        function OK = checkChannelOutput(obj,iChan)
            
        end
            
    end
    
    methods (Static)
        function name =getName()
            name = 'Window Sample Crosscorrelation';
        end   
        function funParams=getDefaultParams(owner,outputDir)
            %----Defaults----%
            nChan = numel(owner.channels_);
            funParams.ChannelIndex = 1:nChan;%Default is all channels
            funParams.MaxLag = floor(owner.nFrames_/4);%Rule-of-thumb for autocorrelation                    
            funParams.UseBands =[];%Default is to use all bands.
            funParams.OutputDirectory = [outputDir  filesep 'sample_crosscorrelation'];                    
            funParams.BatchMode = false;                              
            
        end
    end 
    
end    