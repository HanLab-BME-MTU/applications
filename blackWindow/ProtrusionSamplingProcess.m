classdef ProtrusionSamplingProcess < ImageAnalysisProcess
%
% Process Class for sampling the protrusion vectors which correspond with
% each window for a given movie.
%     
% Hunter Elliott
% 1/2011    
%

    methods (Access = public)
        
        function obj = ProtrusionSamplingProcess(owner,outputDir,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                
                
                super_args{1} = owner;
                super_args{2} = ProtrusionSamplingProcess.getName;
                super_args{3} = @sampleMovieProtrusion;                               
                
                if nargin < 3 || isempty(funParams)                                       
                    
                    %----Defaults----%      
                    
                    funParams.OutputDirectory = [outputDir  filesep 'protrusion_samples'];
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
               obj.outFilePaths_ = {filePath};                               
           end
                        
        end 
        
        function status = checkChannelOutput(obj)
            %Overrides the generic function - there is only one set of prot
            %samples for all channels.            
            status = logical(exist(obj.outFilePaths_{1},'file'));
        end
        
        function prot = loadChannelOutput(obj,varargin)
           
            %Make sure the prot samples are ok
            if ~checkChannelOutput(obj)
                error('Cannot load the protrusion samples - they could not be found!')
            end
            
            outputList = {'','avgNormal'};
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'ProtrusionSamplingProcess'));        
            ip.addParamValue('output','',@(x) all(ismember(x,outputList)));
            ip.parse(obj,varargin{:})
            output = ip.Results.output;
            
            prot = load(obj.outFilePaths_{1});
            fn = fieldnames(prot);
            if numel(fn) > 1, error('Invalid protrusion sample file!'); end            
            prot = prot.(fn{1});
            if ~isempty(output), prot=prot.(output); end
                                    
        end
               
        function h=draw(obj,varargin)
            % Function to draw process output (template method)
            
            if ~ismember('getDrawableOutput',methods(obj)), h=[]; return; end
            outputList = obj.getDrawableOutput();
            ip = inputParser;
            ip.addRequired('obj',@(x) isa(x,'Process'));
            ip.addParamValue('output',outputList(1).var,@(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
            ip.KeepUnmatched = true;
            ip.parse(obj,varargin{:})
			
            data=obj.loadChannelOutput('output',ip.Results.output);
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
        
        
    end
    methods (Static)
        function name =getName()
            name = 'Protrusion Sampling';
        end
        function name =GUI()
            name =@protrusionSamplingProcessGUI;
        end
        function output = getDrawableOutput()
            output(1).name='Protrusion map';
            output(1).var='avgNormal';
            output(1).formatData=[];
            output(1).type='movieGraph';
            output(1).defaultDisplayMethod=@(x)ScalarMapDisplay('Colormap','jet',...
                'ScaleLabel','pixels/frame','Labels',{'Frame number','Window number'});
        end
        
        
    end 
end
    