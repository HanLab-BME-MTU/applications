classdef MaskGeometry3DProcess < ImageAnalysisProcess
%MASKGEOMETRY3DPROCESS process for analyzing 3D mask geometry with analyze3DMovieMaskGeometry.m
%
%
% Hunter Elliott
% 3/2011
%

    properties (SetAccess = protected, GetAccess = public)
        
        
    end
    
    methods (Access = public)
        
        function obj = MaskGeometry3DProcess(owner,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                                                
                
                super_args{1} = owner;
                super_args{2} = 'MaskGeometry';
                super_args{3} = @analyze3DMovieMaskGeometry;                               
                
                if nargin < 2 || isempty(funParams)                                       
                    
                    funParams = MaskGeometry3DProcess.getDefaultParams(owner);                 
                                    
                end
                
                super_args{4} = funParams;    
                                
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
        function OK = checkChannelOutput(obj,iChan)
            
           %Checks if the selected channel has valid output files           
           if nargin < 2 || isempty(iChan) || numel(iChan) > 1 || ... 
                   ~all(obj.checkChanNum(iChan));
               error('You must specify a valid channel number!')
           end
           %Makes there are the right numbr of files in each output
           %directory
           OK =  arrayfun(@(x)(exist(obj.outFilePaths_{x},'dir') && ...
                             numel(dir([obj.outFilePaths_{x} filesep '*.mat']))==obj.owner_.nFrames_),iChan);
        end
        
        
        function mg = loadChannelOutput(obj,iChan,iFrame)
            
            if nargin < 3 || isempty(iFrame) || isempty(iChan)
                error('You must input a frame number to load output for!')
            elseif ~obj.checkChanNum(iChan) || numel(iChan) > 1
                error('You must specify a single, valid channel number!')
            elseif iFrame < 1 || iFrame > obj.owner_.nFrames_
                error('Invalid frame number!')
            elseif ~obj.checkChannelOutput(iChan)
                error('Specified channel does not have valid mask geometry files!')
            end
            
            fileNames = dir([obj.outFilePaths_{iChan} filesep '*.mat']);
            
            mg = load([obj.outFilePaths_{iChan} filesep fileNames(iFrame).name]);
            
            fNames = fieldnames(mg);
            if numel(fNames) ~=1
                error('Invalid mask geometry file!')
            end
            mg = mg.(fNames{1});
        end
        
        function sanityCheck(obj)
            
        end
           
    end
    methods(Static)
        function getName
        end        
        function h = GUI()   
            h = [];
        end
        function procClasses = getConcreteClasses()
            procClasses = [];
        end        
        function funParams = getDefaultParams(owner)
              %----Defaults----%      

                    funParams.ChannelIndex = 1;
                    funParams.SmoothIter =[];%Use the analyze3DMaskGeometry.m defaults.                    
                    funParams.PhysicalUnits = false;
                    funParams.OutputDirectory = ...
                        [owner.outputDirectory_  filesep 'mask_geometry'];
                    funParams.BatchMode = false;                                                   
        end
    end            
    
end