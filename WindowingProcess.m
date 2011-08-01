classdef WindowingProcess < ImageAnalysisProcess
%WINDOWINGPROCESS is a process for creating sampling windows with getMovieWindows.m     
%     
%     
% Hunter Elliott
% 7/2010    
%

    properties (SetAccess = protected, GetAccess = public)
        
        %Window number statistics
        nBandMax_
        nSliceMax_                
        
    end

    methods (Access = public)
        
        function obj = WindowingProcess(owner,outputDir,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                
                
                super_args{1} = owner;
                super_args{2} = WindowingProcess.getName;
                super_args{3} = @getMovieWindows;                               
                
                if nargin < 3 || isempty(funParams)                                       
                    
                    %----Defaults----%      
                    
                    nChan = numel(owner.channels_);
                    funParams.ChannelIndex = 1:nChan;%Default is to combine masks from all channels
                    funParams.SegProcessIndex = []; %No Default.
                    funParams.MethodName = 'ConstantNumber';   
                    funParams.PDEPar = []; %No default.
                    funParams.NonLinearSolver = 'off';
                    funParams.ParaSize = 10;
                    funParams.MeshQuality = []; %Use function default
                    funParams.PerpSize = 10;
                    funParams.ReInit = Inf;                    
                    funParams.StartPoint = []; %No default
                    funParams.MinSize = 10; %Minimum number of pixels a mask object must have to be windowed.
                    funParams.OutputDirectory = [outputDir  filesep 'windows'];
                    funParams.BatchMode = false;                                                      
                                    
                end
                
                super_args{4} = funParams;    
                                
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
        
        function wins = loadChannelOutput(obj,iFrame)
            
            if nargin < 2 || isempty(iFrame)
                error('You must specify a frame number to load windows for!');
            elseif round(abs(iFrame)) ~= iFrame || iFrame > obj.owner_.nFrames_
                error('The frame number must be a positive integer less than or equal to the number of frames in the movie!')
            end
            
            winNames = dir([obj.outFilePaths_ filesep '*.mat']);
            
            if numel(winNames) < iFrame
                error('The window file for the specified frame does not exist!')
            end
            
            tmp = load([obj.outFilePaths_ filesep winNames(iFrame).name]);
            fNames = fieldnames(tmp);
            if numel(fNames) ~= 1
                error('Invalid window file !');
            end
            wins = tmp.(fNames{1});
            
        end
        
        function setOutFilePath(obj,filePath)
            %Overloads the method from ImageAnalysisProcess because there
            %is only one set of windows for all channels.            
                                    
           if ~exist(filePath,'dir')
               error('lccb:set:fatal',...
                   'The directory specified for output is invalid!') 
           else
               obj.outFilePaths_ = filePath;     
               
           end
            
            
        end 
        function status = checkChannelOutput(obj)
           %Overloads the generic function - only one set for all channels
           %Make sure the windows exist and are valid
           if ~exist(obj.outFilePaths_,'dir') || ...
                   numel(dir([obj.outFilePaths_ filesep '*.mat'])) ~= obj.owner_.nFrames_;
               status = false;
           else
               status = true;
           end
            
        end
        
        function setWinStats(obj,nSliceMax,nBandMax)
            
            if nargin < 3 || isempty(nBandMax) || ...
                    isempty(nSliceMax)  || ...
                    any([nBandMax nSliceMax]<1)
                error('Invalid window numbers!')
            end
                        
            obj.nSliceMax_ = nSliceMax;
            obj.nBandMax_ = nBandMax;                        
            
        end
    end
    methods (Static)
        function name =getName()
            name = 'Windowing';
        end
        function h = GUI()
            h= @windowingProcessGUI;
        end
    end
end
    