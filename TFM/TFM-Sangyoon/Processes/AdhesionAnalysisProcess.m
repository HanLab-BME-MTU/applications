classdef AdhesionAnalysisProcess < DataProcessingProcess %& DataProcessingProcess
    
    methods (Access = public)
    
        function obj = AdhesionAnalysisProcess(owner, varargin)
    
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir', owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = AdhesionAnalysisProcess.getName;
                super_args{3} = @analyzeAdhesionMaturation;
                
                if isempty(funParams)
                    funParams = AdhesionAnalysisProcess.getDefaultParams(owner,outputDir);
                end
                
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
                        
        end

        function sanityCheck(obj)
            
            sanityCheck@DataProcessingProcess(obj);
            
            % Cell Segmentation Check
            if obj.funParams_.ApplyCellSegMask                
                iProc = obj.funParams_.SegCellMaskProc;
                % Check Mask is available
                if ~isempty(iProc)
                    assert(iProc < length(obj.owner_.processes_), 'Invalid Process # for Cell Mask Process');
                    assert(isa(obj.owner_.getProcess(iProc), 'MaskRefinementProcess'), ['Process: ' num2str(iProc) ' not a MaskRefinementProcess!']);
                    maskProc = obj.owner_.getProcess(iProc);
                    assert(maskProc.checkChannelOutput(obj.funParams_.ChannelIndex), 'Cell Segmentation Mask Output Not found');
                    % iProc = obj.owner_.getProcessIndex('MaskProcess', 'askUser', false, 'nDesired', Inf);
                else
                    iProc = obj.owner_.getProcessIndex('MaskRefinementProcess');   
                    assert(~isempty(iProc), 'MaskRefinementProcess Process not found cannot Apply Cell Mask!');
                    disp('Setting Cell Segmentation Mask Process index');
                    obj.funParams_.SegCellMaskProc = iProc;
                end
            else
                warning('You do not have segmentation process run for a mask. Using entire image field ...')
            end
            
            %% Sanity check for Detecting FAs 
            iProc = obj.funParams_.detectedNAProc;
            if ~isempty(iProc)
                assert(iProc < length(obj.owner_.processes_), ['Invalid Process #' num2str(iProc) 'for FA Detection Process']);
                assert(isa(obj.owner_.getProcess(iProc), 'DetectionProcess'));
            else
                iProc = obj.owner_.getProcessIndex('DetectionProcess');
                assert(~isempty(iProc), 'FA DetectionProcess Process not found!');
                disp('Setting FA DetectionProcess index');
                obj.funParams_.detectedNAProc = iProc;
            end
            
            %% Sanity check for Tracking FAs 
            iProc = obj.funParams_.trackFAProc;
            if ~isempty(iProc)
                assert(iProc < length(obj.owner_.processes_), ['Invalid Process #' num2str(iProc) 'for FA Tracking Process']);
                assert(isa(obj.owner_.getProcess(iProc), 'TrackingProcess'));
            else
                iProc = obj.owner_.getProcessIndex('TrackingProcess');
                assert(~isempty(iProc), 'FA TrackingProcess Process not found!');
                disp('Setting FA TrackingProcess index');
                obj.funParams_.trackFAProc = iProc;
            end
            
            %% Add Sanity check for FAsegmentationProcess Mask
            iProc = obj.funParams_.FAsegProc;
            if ~isempty(iProc)
                assert(iProc < length(obj.owner_.processes_), ['Invalid Process #' num2str(iProc) 'for FA Seg Process']);
                assert(isa(obj.owner_.getProcess(iProc), 'FocalAdhesionSegmentationProcess'));
            else
                iProc = obj.owner_.getProcessIndex('FocalAdhesionSegmentationProcess');
                assert(~isempty(iProc), 'FocalAdhesionSegmentationProcess Process not found!');
                disp('Setting FocalAdhesionSegmentationProcess index');
                obj.funParams_.FAsegProc = iProc;
            end
        end
        
        function varargout = loadChannelOutput(obj, iChan, varargin)
            % Input check
            outputList = {'trackFC','trackNA','trackFA','detBA',...
                          'detectedFA','detFA','detFC','detNA',...
                          'adhboundary_FA', 'adhboundary_FC'};

            ip =inputParser;
            ip.addRequired('obj');
            ip.addRequired('iChan', @(x) obj.checkChanNum(x));
            ip.addOptional('iFrame', [] ,@(x) obj.checkFrameNum(x));
            ip.addParameter('useCache',true,@islogical);
            ip.addParameter('output', outputList{3}, @(x) all(ismember(x,outputList)));
            ip.parse(obj,iChan,varargin{:})
            output = ip.Results.output;
            varargout = cell(numel(output), 1);
            iFrame = ip.Results.iFrame;
            if ischar(output),output={output}; end
            
            % Data loading
            s = cached.load(obj.outFilePaths_{1,iChan}, '-useCache', ip.Results.useCache, 'tableTracksNA');
%             st = cached.load(obj.outFilePaths_{1,iChan}, '-useCache', ip.Results.useCache, 'tracksNA');
            % Note, could do a stack.
            
            %% Check struct vs table loading           
            if isstruct(s)
                s = s.tableTracksNA;
            else
                disp('loaded as table');
            end

            nTracks = length(s.xCoord(:,iFrame));
            number = (1:length(s.xCoord(:,iFrame)))';
            state = categorical(s.state(:,iFrame));
            iiformat = ['%.' '3' 'd'];
            
            for iout = 1:numel(output)
                switch output{iout}
                    case 'detectedFA'  
                        varargout{1} = t;
                    case 'detBA' 
                        validState = state == 'BA';
                    case {'detNA', 'trackNA'}
                        validState = state == 'NA';
                    case {'detFC', 'trackFC', 'adhboundary_FC'}
                        validState = state == 'FC';
                    case {'detFA', 'trackFA', 'adhboundary_FA'}
                        validState = state == 'FA';
                    case 'staticTracks'
                    otherwise
                        error('Incorrect Output Var type');
                end   
                if ~isempty(strfind(output{iout}, 'det'))
                
                    t = table(s.xCoord(:,iFrame), s.yCoord(:,iFrame));
                    varargout{iout} = t{validState,:};                                 
                
                elseif ~isempty(strfind(output{iout},'track'))
                    
                    vars = {'xCoord', 'yCoord', 'number'};
                    validTracks = validState & s.startingFrameExtra <= iFrame & s.endingFrameExtra >= iFrame;                    
                    st = table(s.xCoord(:,1:iFrame), s.yCoord(:,1:iFrame), number, ...
                               'VariableNames', {'xCoord', 'yCoord', 'number'});                    
                    
                    varargout{iout}(nTracks, 1) = struct('xCoord', [], 'yCoord', [], 'number', []);
                    varargout{iout}(validTracks, :) = table2struct(st(validTracks, vars));
                
                elseif ~isempty(strfind(output{iout},'adhboundary'))                    
                
%                     adhBoundary = cellfun(@(x) x{iFrame}, s{validState, 'adhBoundary'}, 'UniformOutput', false);                         
                    p=obj.funParams_;
                    labelTifPath = [p.OutputDirectory filesep 'labelTifs'];
                    labelAdhesion = imread(strcat(labelTifPath,'/label',num2str(iFrame,iiformat),'.tif'));
                    labelAdhesion = bwlabel(labelAdhesion>0,4);
                    maxLabel=max(labelAdhesion(:));
                    adhBound = cell(1,maxLabel);
                    for ii=1:maxLabel
                        curAdhBound = bwboundaries(labelAdhesion==ii,4,'noholes');
                        adhBound{ii} = curAdhBound{1}; % strongly assumes each has only one boundary
                    end
                    validAdhState = cellfun(@(x) x(iFrame),s.refineFAID(validState));
                    
                    varargout{iout} = adhBound(validAdhState); %table2struct(table(adhBoundary, number(validState),'VariableNames',{'adhBoundary', 'number'}));                                 
                else
                    varargout{iout} = [];
                end
            end
        end      
    end


    methods (Static)
        function name = getName()
            name = 'Focal Adhesion Analysis';
        end

        function h = GUI()
            h = @focalAdhesionAnalysisProcessGUI;
        end
        
        function output = getDrawableOutput()
            ii = 10;
%             i = ii-1; output(i).name='Before Adhesion Detection'; 
%             output(i).var='detBA';
%             output(i).formatData=[];
%             output(i).type='overlay';
%             output(i).defaultDisplayMethod=@(x) LineDisplay('Marker','.',...
%                 'LineStyle', 'none', 'LineWidth', .6, 'Color', 'g',...
%                 'MarkerSize', 5);            
%             i = ii-2; output(i).name='Nascent Adhesion Detection'; 
%             output(i).var='detNA';
%             output(i).formatData=[];
%             output(i).type='overlay';
%             output(i).defaultDisplayMethod=@(x) LineDisplay('Marker','o',...
%                 'LineStyle','none', 'LineWidth', .8, 'Color', 'r'); 
%             i = ii-3; output(i).name='Focal Contact Detection'; 
%             output(i).var='detFC';
%             output(i).formatData=[];
%             output(i).type='overlay';
%             output(i).defaultDisplayMethod=@(x) LineDisplay('Marker','o',...
%                 'LineStyle','none', 'LineWidth', .8, 'Color', [255/255 153/255 51/255]); 
%             i = ii-4; output(i).name='Focal Adhesion Detection'; 
%             output(i).var='detFA';
%             output(i).formatData=[];
%             output(i).type='overlay';
%             output(i).defaultDisplayMethod=@(x) LineDisplay('Marker','o',...
%                 'LineStyle','none', 'LineWidth', .8, 'Color', 'b'); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Tracks Display
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            i = ii-5; output(i).name='Nascent Adhesion Tracks'; 
            output(i).var='trackNA';
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) FATracksDisplay('Linewidth', 1, 'Color', 'r'); 
            i = ii-6; output(i).name='Focal Contact Tracks'; 
            output(i).var='trackFC';
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) FATracksDisplay('Linewidth', 1, 'Color', [255/255 153/255 51/255]); 
            i = ii-7; output(i).name='Focal Adhesion Tracks'; 
            output(i).var='trackFA';
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) FATracksDisplay('Linewidth', 1, 'Color', 'b'); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Adhesion Boundaries
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
            i = ii-9; output(i).name='Focal Adhesion Boundary';
            output(i).var='adhboundary_FA';
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) AdhBoundaryDisplay('Color', 'b');
            
            i = ii-8; output(i).name='Focal Contact Boundary';
            output(i).var='adhboundary_FC';
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) AdhBoundaryDisplay('Color', [255/255 153/255 51/255]);

        end       

        function funParams = getDefaultParams(owner, varargin)

            % MD.getPackage(MD.getPackageIndex('FocalAdhesionPackage')).outputDirectory_
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.addOptional('iChan', 1:numel(owner.channels_),...
               @(x) all(owner.checkChanNum(x)));
            ip.parse(owner, varargin{:});
            
            % Set default parameters
            funParams.OutputDirectory = [ip.Results.outputDir filesep 'AdhesionAnalysis'];
            funParams.ChannelIndex = ip.Results.iChan;

            funParams.ApplyCellSegMask = true;
            funParams.SegCellMaskProc = []; % Specify Process with cell mask output
            funParams.detectedNAProc = []; % Specify FA detection Process index
            funParams.trackFAProc = []; % Specify FA tracking Process index
            funParams.FAsegProc = []; % Specify FA segmentation Process index           

            funParams.onlyEdge = false; 
            funParams.matchWithFA = true; 
            funParams.minLifetime = 5;  % For tracks
            funParams.reTrack = true;
            funParams.getEdgeRelatedFeatures = true;
            funParams.bandwidthNA = 7;
            
            %% TODO - likely will remove this.
            funParams.backupOldResults = true;           
        end

    end
end
