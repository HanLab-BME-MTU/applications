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
        
        
        function h=draw(obj,iFrame,varargin)
            % Function to draw process output (template method)
            
            if ~ismember('getDrawableOutput',methods(obj)), h=[]; return; end
            outputList = obj.getDrawableOutput();
            ip = inputParser;
            ip.addRequired('obj',@(x) isa(x,'Process'));
            ip.addRequired('iFrame',@isnumeric);
            ip.addParamValue('output',outputList(1).var,@(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
            ip.KeepUnmatched = true;
            ip.parse(obj,iFrame,varargin{:})
			
            data=obj.loadChannelOutput(iFrame,'output',ip.Results.output);
            iOutput= find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
            if ~isempty(outputList(iOutput).formatData),
                data=outputList(iOutput).formatData(data);
            end
            try
                assert(~isempty(obj.displayMethod_{iOutput}));
            catch ME
                obj.displayMethod_{iOutput}=...
                    outputList(iOutput).defaultDisplayMethod();
            end
            
            % Delegate to the corresponding method
            tag = [obj.getName '_output' num2str(iOutput)];
            drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                2*numel(fieldnames(ip.Unmatched)),1);
            h=obj.displayMethod_{iOutput}.draw(data,tag,drawArgs{:});
        end
        
        
        function output = getDrawableOutput(obj)
            output(1).name='Images';
            output(1).var='';
            output(1).formatData=[];
            output(1).type='movieImage';
            output(1).defaultDisplayMethod=@(x) ImageDisplay('Colormap','jet',...
                'Colorbar','on','Units','','CLim',obj.getIntensityLimits());
        end
            
    end
    
    methods (Access=protected)
        function limits = getIntensityLimits(obj)
            allImages=arrayfun(@(x)loadChannelOutput(obj,x),1:obj.owner_.nFrames_,...
                'UniformOutput',false);
            allImages = vertcat(allImages{:});
            limits=[min(allImages(:)) max(allImages(:))];
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
       
    