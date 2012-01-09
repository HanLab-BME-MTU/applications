classdef SegmentationProcess3D < SegmentationProcess
    
    %Class for segmenting 3D movies using segment3Dmovie
    
    methods (Access = public)
        
        function obj = SegmentationProcess3D(owner,funParams)
             
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = '3DSegmentation';
                super_args{3} = @segment3DMovie;                           
                
                if nargin < 2 || isempty(funParams)                                       
                   
                    funParams = SegmentationProcess3D.getDefaultParams(owner);
                   
                end
                %Make sure the input parameters are legit??
                super_args{4} = funParams;                    
            end
            
            obj = obj@SegmentationProcess(super_args{:});
            
        end               
            
    end
    methods(Static)
        function name =getName()
            name = 'SegmentationProcess3D';
        end
        function h = GUI()   
            h = [];
        end
        function procClasses = getConcreteClasses()
            procClasses = [];
        end        
        function funParams = getDefaultParams(owner)
            
             %----Defaults----%
            funParams.OutputDirectory = ...
                [owner.outputDirectory_  filesep 'masks'];                          
            funParams.ChannelIndex = 1 : numel(owner.channels_);    
            funParams.ThresholdValue = []; %Default is auto-thresh
            funParams.PreFilterSig = 0;
            funParams.Method = 'HuntThresh';
            funParams.FixJumps = true;
            funParams.PostProcess = true;                    
            funParams.MinVolume = 25;
            funParams.NumObjects = 1; 
            funParams.ClosureRadius = 2;
            funParams.FillHoles = 2;
            funParams.BatchMode = false;                                              
        end
        
        
    end
        
end