classdef SampleAutocorrelationProcess < DataProcessingProcess
%Process class definition for calculateWindowSampleAutocorrelation.m
%     
% Hunter Elliott
% 8/2011
%
    methods (Access = public)
        
        function obj = SampleAutocorrelationProcess(owner,outputDir,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                
                
                super_args{1} = owner;

                super_args{2} = SampleAutocorrelationProcess.getName;
                super_args{3} = @calculateWindowSampleAutocorrelation;                               
                
                if nargin < 3 || isempty(funParams)                                       
                    
                    %----Defaults----%
                    nChan = numel(owner.channels_);
                    funParams.ChannelIndex = 1:nChan;%Default is all channels
                    funParams.MaxLag = floor(owner.nFrames_/4);%Rule-of-thumb for autocorrelation
                    funParams.DetrendMethod = 'linear';%Default is linear trend removal.
                    funParams.OutputDirectory = [outputDir  filesep 'sample_autocorrelation'];
                    funParams.BatchMode = false;                              

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
            name = 'Window Sample Autocorrelation';
        end                
    end 
    
end           
    