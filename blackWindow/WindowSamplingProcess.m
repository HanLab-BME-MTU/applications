classdef WindowSamplingProcess < ImageAnalysisProcess
%Process     
%     
% Hunter Elliott
% 7/2010    
%

    methods (Access = public)
        
        function obj = WindowSamplingProcess(owner,outputDir,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                
                
                super_args{1} = owner;

                super_args{2} = WindowSamplingProcess.getName;
                super_args{3} = @sampleMovieWindows;                               
                
                if nargin < 3 || isempty(funParams)                                       
                    
                    %----Defaults----%      
                    nChan = numel(owner.channels_);
                    funParams.ChannelIndex = 1:nChan;%Default is to sample all channels                    
                    funParams.ProcessIndex = [];%Default is to use raw images
                    funParams.OutputDirectory = [outputDir  filesep 'window_sampling'];
                    funParams.BatchMode = false;                                                                                
                                    
                end
                
                super_args{4} = funParams;    
                                
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end   
        
        function samp = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'ImageAnalysisProcess'));
            ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
            ip.addOptional('iFrame',@(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParamValue('output',[],@ischar);            
            ip.parse(obj,iChan,varargin{:})      
                        
            tmp = load(obj.outFilePaths_{iChan});
            fNames = fieldnames(tmp);
            if numel(fNames) ~= 1
                error('Invalid window sample file !');
            end
            samp = tmp.(fNames{1});
            
            
        end
        
        function OK = checkChannelOutput(obj,iChan)
            
           %Checks if the selected channels have valid output files
           nChanTot = numel(obj.owner_.channels_);
           if nargin < 2 || isempty(iChan)
               iChan = 1:nChanTot;
           end
           %Makes sure there's at least one .mat file in the speified
           %directory
           OK =  arrayfun(@(x)(x <= nChanTot && ...
                             x > 0 && isequal(round(x),x) && ...
                             exist(obj.outFilePaths_{x},'file')),iChan);
        end                
        
        function output = getDrawableOutput(obj)
            output(1).name='Images';
            output(1).var='';
            output(1).formatData=[];
            output(1).type='image';
            output(1).defaultDisplayMethod=@(x) ImageDisplay('Colormap','jet',...
                'Colorbar','on','Units','','CLim',obj.getIntensityLimits(x));
        end
            
    end
    
    methods (Access=protected)
        function limits = getIntensityLimits(obj,iChan)
            ratioImages=arrayfun(@(x)loadChannelOutput(obj,iChan,x),1:obj.owner_.nFrames_,...
                'UniformOutput',false);
            allRatioImages = vertcat(ratioImages{:});
            limits=[min(allRatioImages(:)) max(allRatioImages(:))];
        end   
    end
    methods (Static)
        function name =getName()
            name = 'Window Sampling';
        end
        function name= GUI()
            name =@windowSamplingProcessGUI;
        end
        
    end 
    
end
       
    