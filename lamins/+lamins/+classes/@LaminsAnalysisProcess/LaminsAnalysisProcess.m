classdef LaminsAnalysisProcess < ImageAnalysisProcess
    %LaminsAnalysisProcess Analyzes lamin curvilinear structures
    
    properties
    end
    
    methods
        function obj = LaminsAnalysisProcess(owner)
%             obj = obj@ImageAnalysisProcess(owner,'LaminsAnalysisProcess',@lamins.functions.analyzeLaminsForProcess);
              obj = obj@ImageAnalysisProcess(owner, ...
                  'LaminsAnalysisProcess', ...
                  @lamins.functions.analyzeLaminsForProcess, ...
                  struct());
            
        end
        function run(obj,varargin)
            basePath = [obj.owner_.outputDirectory_ filesep 'LaminsAnalysisProcess'];
            
            changed = false;
            for c=1:length(obj.owner_.channels_)
                outPath = [basePath filesep 'channel_' num2str(c)];
                if(isempty(obj.outFilePaths_{c}))
                    obj.setOutFilePaths(outPath,c);
                    changed = true;
                end
            end
            if(changed)
                obj.owner_.save();
            end
            
            params = obj.getParameters();
            if(~isfield(params,'outFilePaths') || isempty(params.outFilePaths))
                params.outFilePaths = obj.outFilePaths_;
            end
            params.process = obj;
            
            % Call the superclass runner
            obj.run@ImageAnalysisProcess(params,varargin{:});
        end
    end
    methods (Static)
        function name = getName(varargin)
            name = 'Lamins Analysis Process';
        end
        function params = getDefaultParams(varargin)
            ip = lamins.functions.analyzeLaminsForProcessParameters();
            ip.parse(varargin{:});
            params = ip.Results;
        end
    end
    
end

