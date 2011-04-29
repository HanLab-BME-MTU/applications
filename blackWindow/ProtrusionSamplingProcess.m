classdef ProtrusionSamplingProcess < ImageAnalysisProcess
%
% Process Class for sampling the protrusion vectors which correspond with
% each window for a given movie.
%     
% Hunter Elliott
% 1/2011    
%

    methods (Access = public)
        
        function obj = ProtrusionSamplingProcess(owner,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                
                
                super_args{1} = owner;
                super_args{2} = 'ProtrusionSampling';
                super_args{3} = @sampleMovieProtrusion;                               
                
                if nargin < 2 || isempty(funParams)                                       
                    
                    %----Defaults----%      
                    
                    funParams.OutputDirectory = ...
                        [owner.outputDirectory_  filesep 'protrusion_samples'];
                    funParams.BatchMode = false;                                                      
                                    
                end
                
                super_args{4} = funParams;    
                                
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
                
        function setOutFilePath(obj,filePath)
            %Overloads the method from ImageAnalysisProcess because there
            %is only one set of vectors for all channels, which is stored
            %as a single file
                                    
           if ~exist(filePath,'file')
               error('lccb:set:fatal',...
                   'The file specified specified as output for the function is invalid!') 
           else
               obj.outFilePaths_ = filePath;                               
           end
                        
        end 
        
        function status = checkChannelOutput(obj)
            %Overrides the generic function - there is only one set of prot
            %samples for all channels.            
            status = false;
            if exist(obj.outFilePaths_,'file')                                
                status = true;                                                   
            end
        end
        
        function prot = loadChannelOutput(obj)
           
            %Make sure the prot samples are ok
            if ~checkChannelOutput(obj)
                error('Cannot load the protrusion samples - they could not be found!')
            end
            
            prot = load(obj.outFilePaths_);
            fn = fieldnames(prot);
            if numel(fn) > 1
                error('Invalid protrusion sample file!')
            end            
            prot = prot.(fn{1});
                                    
        end
        
    end    
end
    