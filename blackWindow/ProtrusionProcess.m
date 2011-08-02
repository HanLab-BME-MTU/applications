classdef ProtrusionProcess < ImageAnalysisProcess
    %
    % Process Class for calculating protrusion vectors using the
    % getMovieProtrusion.m wrapper function.
    %
    % Hunter Elliott
    % 8/2010
    %
    
    methods (Access = public)
        
        function obj = ProtrusionProcess(owner,outputDir,funParams)
            
            if nargin == 0
                super_args = {};
            else
                
                super_args{1} = owner;
                super_args{2} = ProtrusionProcess.getName;
                super_args{3} = @getMovieProtrusion;
                
                if nargin < 3 || isempty(funParams)
                    
                    %----Defaults----%
                    
                    nChan = numel(owner.channels_);
                    funParams.ChannelIndex = 1:nChan;%Default is to combine masks from all channels
                    funParams.SegProcessIndex = [];%No default.
                    funParams.DownSample = 50;
                    funParams.SplineTolerance = 30;%This is the default in protrusionAnalysis, so I use it also.
                    funParams.OutputDirectory = [outputDir filesep 'protrusion'];
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
        
        function varargout = loadChannelOutput(obj,varargin)
            
            %Make sure the prot vectors are ok
            if ~checkChannelOutput(obj)
                error('Cannot load the protrusion vectors - they could not be found!')
            end
            
%             prot = load(obj.outFilePaths_);
            
            outputList = {'protrusion','normals','smoothedEdge'};
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'ProtrusionProcess'));        
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,...
                @(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParamValue('output',{},@(x) all(ismember(x,outputList)));
            ip.parse(obj,varargin{:})
            iFrame = ip.Results.iFrame;
            
            % Data loading
            output = ip.Results.output;
            if ischar(output), output = {output}; end
            s = load(obj.outFilePaths_,output{:});
            
            if numel(ip.Results.iFrame)>1,
                if isempty(output)
                   varargout{1}=s;
                else
                    for i=1:numel(output),
                        varargout{i}=s.(output{i});
                    end
                end
            else
                varargout{1} = [s.smoothedEdge{iFrame} s.smoothedEdge{iFrame}+s.protrusion{iFrame}];
            end
                        
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
        
    end
    methods (Static)
        function name =getName()
            name = 'Protrusion';
        end
        function h = GUI()
            h= @protrusionProcessGUI;
        end
        
        function output = getDrawableOutput()
            output(1).name='Protrusion vectors';
            output(1).var={'smoothedEdge','protrusion'};
            output(1).formatData=@(x) x(:,[2 1 4 3]);
            output(1).type='movieOverlay';
            output(1).defaultDisplayMethod=@(x) VectorFieldDisplay('Color','r');
        end
    end
end