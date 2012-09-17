classdef WindowSamplingProcess < ImageAnalysisProcess
    %Process
    %
    % Hunter Elliott
    % 7/2010
    %
    
    methods (Access = public)
        
        function obj = WindowSamplingProcess(owner,varargin)
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = WindowSamplingProcess.getName;
                super_args{3} = @sampleMovieWindows;
                if isempty(funParams)
                    funParams=WindowSamplingProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
                
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
        
        function samp = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'','avg'};
            nOutput = numel(obj.funParams_.ProcessIndex);
            ip =inputParser;
            ip.addRequired('iChan',@(x) obj.checkChanNum(x));
            ip.addOptional('iOutput',1,@(x) ismember(x,1:nOutput));
            ip.addParamValue('output','',@(x) all(ismember(x,outputList)));
            ip.parse(iChan,varargin{:})
            iOutput =ip.Results.iOutput;
            output=ip.Results.output;
            
            s = load(obj.outFilePaths_{iOutput,iChan});
            fNames = fieldnames(s);
            assert(numel(fNames) == 1,'Invalid window sample file !');
            samp = s.(fNames{1});
            
            if ~isempty(output), samp=samp.(output); end
        end
        
        function status = checkChannelOutput(obj,varargin)
            
            %Checks if the selected channels have valid output files
            nChanTot = numel(obj.owner_.channels_);
            nOutput = numel(obj.funParams_.ProcessIndex);
            ip = inputParser;
            ip.addOptional('iChan',1:nChanTot,@(x) obj.checkChanNum(x));
            ip.addOptional('iOutput',1:nOutput,@(x) ismember(x,1:nOutput));
            ip.parse(varargin{:});
            iChan=ip.Results.iChan;
            iOutput=ip.Results.iOutput;
            
            %Makes sure there's at least one .mat file in the speified
            %directory
            status =  cellfun(@(x)logical(exist(x,'file')),obj.outFilePaths_(iOutput,iChan));
        end
        
        function h=draw(obj,iChan,varargin)
            % Function to draw process output (template method)
            
            if ~ismember('getDrawableOutput',methods(obj)), h=[]; return; end
            outputList = obj.getDrawableOutput();
            ip = inputParser;
            ip.addRequired('iChan',@isscalar);
            ip.addOptional('iOutput',1,@isscalar);
            ip.KeepUnmatched = true;
            ip.parse(iChan,varargin{:})
            iOutput = ip.Results.iOutput;
            
            data=obj.loadChannelOutput(iChan,iOutput,'output','avg');
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
            tag = ['process' num2str(obj.getIndex) '_channel' num2str(iChan) '_output' num2str(iOutput)];
            drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                2*numel(fieldnames(ip.Unmatched)),1);
            h=obj.displayMethod_{iOutput,iChan}.draw(data,tag,drawArgs{:});
        end
        
        function output = getSampledOutput(obj)
            % Construct linear output of time series (to be used by
            % integrator)
            processIndex = obj.funParams_.ProcessIndex;
            channelIndex = obj.funParams_.ChannelIndex;
            outputName = obj.funParams_.OutputName;
            if ~iscell(processIndex), processIndex={processIndex}; end
            if ~iscell(channelIndex), processIndex={processIndex}; end
            if ~iscell(outputName), outputName={outputName}; end
            
            nChanPerOutput = cellfun(@numel,channelIndex);
            nOutputTot = sum(nChanPerOutput);
            output(nOutputTot,1)=struct();
            for i=1:numel(processIndex)
                procId = processIndex{i};
                if isempty(procId)
                    processName = '';
                    name='Raw images';
                    
                else
                    processName = class(obj.owner_.processes_{procId});
                    parentOutput = obj.owner_.processes_{procId}.getDrawableOutput;
                    iOutput = strcmp(outputName{i},{parentOutput.var});
                    name=parentOutput(iOutput).name;
                end
                
                for j =1:numel(channelIndex{i})
                    iInput = sum(nChanPerOutput(1:i-1))+j;
                    output(iInput).processName = processName;
                    output(iInput).processIndex = procId;
                    output(iInput).channelIndex = channelIndex{i}(j);
                    output(iInput).var = 'avg';
                    output(iInput).name = [name ' - Channel ' num2str(channelIndex{i}(j))];
                end
            end
        end
        
        function drawableOutput = getDrawableOutput(obj)
            % Build the list of drawable output from 
            processIndex = obj.funParams_.ProcessIndex;
            outputName = obj.funParams_.OutputName;
            if ~iscell(processIndex), processIndex={processIndex}; end
            if ~iscell(outputName), outputName={outputName}; end
            
            % Initialize drawable output array
            nOutput = numel(processIndex);
            drawableOutput(nOutput,1)=struct();
            
            timeInterval = obj.owner_.timeInterval_;
            if ~isempty(timeInterval)
                scaling = [1 timeInterval 1];
                xlabel = 'Time (s)';
            else
                scaling = [1 1 1];
                xlabel = 'Frame number';
            end
            
            for i=1:nOutput
                procId = processIndex{i};
                if isempty(procId)
                    % No procId - sampled raw images 
                    drawableOutput(i).name = 'Raw images';
                else
                    % Read output name from parent process & ouput variable
                    parentOutput = obj.owner_.processes_{procId}.getDrawableOutput;
                    iOutput = strcmp(outputName{i},{parentOutput.var});
                    drawableOutput(i).name=parentOutput(iOutput).name;
                end
                
                %
                drawableOutput(i).var='avg';
                drawableOutput(i).formatData=@(x) permute(x,[1 3 2]);
                drawableOutput(i).type='sampledGraph';
                
                % Get process-defined colormap if applicable
                cmap = @(varargin)jet(2^8);
                if ~isempty(procId) && ismember('getColormap',methods(obj.owner_.processes_{procId}))
                    cmap=@(x,i)obj.owner_.processes_{procId}.getColormap(iOutput,obj.getIntensityLimits(x,i));
                end
                
                % Get process-defined units if applicable
                units =@(varargin) '';
                if ~isempty(procId) && ismember('getUnits',methods(obj.owner_.processes_{procId}))
                    units=@(i)obj.owner_.processes_{procId}.getUnits(i);
                end
                
                drawableOutput(i).defaultDisplayMethod = @(x) ScalarMapDisplay(...
                    'Colormap', cmap(x,i), 'Units', units(i),...
                    'CLim', obj.getIntensityLimits(x,i), 'Scaling', scaling,...
                    'Labels', {xlabel, 'Window depth', 'Window number'});
            end
        end
        
        

        
        function limits = getIntensityLimits(obj,iChan,iOutput)
            data=obj.loadChannelOutput(iChan,iOutput,'output','avg');
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
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1:numel(owner.channels_);%Default is to sample all channels
            funParams.ProcessIndex = [];%Default is to use raw images
            funParams.OutputName = '';%Default is to use raw images
            funParams.OutputDirectory = [outputDir  filesep 'window_sampling'];
            funParams.BatchMode = false;
        end
        function samplableInput = getSamplableInput()
            % List process output that can be sampled
            processNames = horzcat('Raw images','DoubleProcessingProcess',...
                repmat({'KineticAnalysisProcess'},1,3),'FlowAnalysisProcess');
            samplableOutput = {'','','netMap','polyMap','depolyMap','speedMap'};
            samplableInput=cell2struct(vertcat(processNames,samplableOutput),...
                {'processName','samplableOutput'});  
        end
        
        
    end
end