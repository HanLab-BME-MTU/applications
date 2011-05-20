classdef ProtrusionProcess < ImageAnalysisProcess
%
% Process Class for calculating protrusion vectors using the
% getMovieProtrusion.m wrapper function.
%     
% Hunter Elliott
% 8/2010    
%

    methods (Access = public)
        
        function obj = ProtrusionProcess(owner,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                
                
                super_args{1} = owner;
                super_args{2} = ProtrusionProcess.getName;
                super_args{3} = @getMovieProtrusion;                               
                
                if nargin < 2 || isempty(funParams)                                       
                    
                    %----Defaults----%      
                    
                    nChan = numel(owner.channels_);
                    funParams.ChannelIndex = 1:nChan;%Default is to combine masks from all channels
                    funParams.SegProcessIndex = [];%No default.
                    funParams.DownSample = 50;
                    funParams.SplineTolerance = 30;%This is the default in protrusionAnalysis, so I use it also.
                    funParams.OutputDirectory = ...
                        [owner.outputDirectory_  filesep 'protrusion'];
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
                   'The file specified as output for the function is invalid!') 
           else
               obj.outFilePaths_ = filePath;                               
           end
                        
        end 
        
        function status = checkChannelOutput(obj)
            %Overrides the generic function - there is only one set of prot
            %vectors for all channels.            
            status = false;
            if exist(obj.outFilePaths_,'file')                                
                tmp = load(obj.outFilePaths_);
                if isfield(tmp,'protrusion') && isfield(tmp,'normals') ...
                        && isfield(tmp,'smoothedEdge')
                    status = true;
                                   
                end
            end
        end
        
        function prot = loadChannelOutput(obj)
           
            %Make sure the prot vectors are ok
            if ~checkChannelOutput(obj)
                error('Cannot load the protrusion vectors - they could not be found!')
            end
            
            prot = load(obj.outFilePaths_);
                        
            
        end
        
    end
    methods (Static)
        function name =getName()
            name = 'Protrusion';
        end
    end
end
    