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
            outputList = {'','avg'};
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'ImageAnalysisProcess'));
            ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
            ip.addParamValue('output','',@(x) all(ismember(x,outputList)));          
            ip.parse(obj,iChan,varargin{:})      
            output=ip.Results.output;      
                  
            tmp = load(obj.outFilePaths_{iChan});
            fNames = fieldnames(tmp);
            if numel(fNames) ~= 1
                error('Invalid window sample file !');
            end
            samp = tmp.(fNames{1});
            if ~isempty(output), samp=samp.(output); end
            
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
        
        function h=draw(obj,iChan,varargin)
            % Function to draw process output (template method)
            
            if ~ismember('getDrawableOutput',methods(obj)), h=[]; return; end
            outputList = obj.getDrawableOutput();
            ip = inputParser;
            ip.addRequired('obj',@(x) isa(x,'Process'));
            ip.addRequired('iChan',@isnumeric);
            ip.addParamValue('output',outputList(1).var,@(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
            ip.KeepUnmatched = true;
            ip.parse(obj,iChan,varargin{:})
            
            data=obj.loadChannelOutput(iChan,'output',ip.Results.output);
            iOutput= find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
            if ~isempty(outputList(iOutput).formatData),
                data=outputList(iOutput).formatData(data);
            end
            try
                assert(~isempty(obj.displayMethod_{iOutput,iChan}));
            catch ME
                obj.displayMethod_{iOutput,iChan}=...
                    outputList(iOutput).defaultDisplayMethod(iChan);
            end
            
            % Delegate to the corresponding method
            tag = [obj.getName '_channel' num2str(iChan) '_output' num2str(iOutput)];
            drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                2*numel(fieldnames(ip.Unmatched)),1);
            h=obj.displayMethod_{iOutput,iChan}.draw(data,tag,drawArgs{:});
        end
        
        
        function output = getDrawableOutput(obj)
            if isempty(obj.funParams_.ProcessIndex)
                output(1).name='Raw image map';
            else
                output(1).name='Activity map';
            end
            output(1).var='avg';
            output(1).formatData=[];
            output(1).type='graph';
            output(1).defaultDisplayMethod=@(x) ScalarMapDisplay('Colormap','jet',...
                'CLim',obj.getIntensityLimits(x),'Labels',{'Frame number','Window depth','Window number'});
        end
            
    end
    
    methods (Access=protected)
        function limits = getIntensityLimits(obj,iChan)
            data=obj.loadChannelOutput(iChan,'output','avg');
            limits=[min(data(:)) max(data(:))];
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
       
    