classdef InitialRiseTimeLagCalculationProcess < DataProcessingProcess
    methods (Access = public)
        function obj = InitialRiseTimeLagCalculationProcess(owner,varargin)
%             obj = obj@DataProcessingProcess(owner, InitialRiseTimeLagCalculationProcess.getName);
%             obj.funName_ = @calculateInitialRiseTimeLagFromTracks; % This should be variation from colocalizationAdhesionWithTFM
%             obj.funParams_ = InitialRiseTimeLagCalculationProcess.getDefaultParams(owner,varargin{1});
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
                super_args{2} = InitialRiseTimeLagCalculationProcess.getName;
                super_args{3} = @calculateInitialRiseTimeLagFromTracks;
                
                if isempty(funParams)
                    funParams = InitialRiseTimeLagCalculationProcess.getDefaultParams(owner,outputDir);
                end
                
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
        end
        
        function [output, output2] = loadChannelOutput(obj, iChan, varargin)
            outputList = {'tracksNA','idClass'};
            nOutput = length(outputList);

            ip = inputParser;
            ip.addRequired('iChan',@(x) obj.checkChanNum(x));
            ip.addOptional('iOutput',1,@(x) ismember(x,1:nOutput));
            ip.addParamValue('output','',@(x) all(ismember(x,outputList)));
            ip.addParamValue('useCache',false,@islogical);
            ip.addParamValue('idSelected',[]);
            ip.parse(iChan,varargin{:})
    
            outputRequested=ip.Results.output;
            p = obj.funParams_;
            %% Reading tracks from master channel
            % Use the previous analysis folder structure
            % It might be good to save which process the tracksNA was obtained.
            MD=obj.owner_;
            if strcmp(outputRequested,'tracksNA')
                disp('Reading tracksNA ...')
                tic
                try
                    iAdhProc = MD.getProcessIndex('AdhesionAnalysisProcess');
                    adhAnalProc = MD.getProcess(iAdhProc);
                catch
                    iFAPack =  MD.getPackageIndex('FocalAdhesionPackage');
                    faPack = MD.getPackage(iFAPack);
                    adhAnalProc = faPack.processes_{7};
                end
                idSelected=ip.Results.idSelected;
                if isempty(idSelected)
                    tracksNA=adhAnalProc.loadChannelOutput(p.ChannelIndex,'output','tracksNA','wantFullTrack',true);
                else
                    tracksNA=adhAnalProc.loadChannelOutput(p.ChannelIndex,'output','tracksNA','idSelected',idSelected,'wantFullTrack',true);
                end
                % numChans = numel(p.ChannelIndex);

                % Now we have to combine this with readings from step 9 and 10
                iFAPack = MD.getPackageIndex('FocalAdhesionPackage');
                FAPack=MD.packages_{iFAPack}; iTheOtherProc=9; iForceRead=10;
                theOtherReadProc=FAPack.processes_{iTheOtherProc};
                forceReadProc=FAPack.processes_{iForceRead};

                if ~isempty(theOtherReadProc)
                    ampObj = load(theOtherReadProc.outFilePaths_{1,p.ChannelIndex},'tracksAmpTotal'); % the later channel has the most information.
                    tracksAmpTotal = ampObj.tracksAmpTotal;
                    if ~isempty(idSelected)
                        tracksAmpTotal = tracksAmpTotal(idSelected);
                    end
                %     tracksAmpTotal = theOtherReadProc.loadChannelOutput(p.ChannelIndex);

                    if isfield(tracksAmpTotal,'ampTotal2')
                        [tracksNA(:).ampTotal2] = tracksAmpTotal.ampTotal2;
                    end
                    if isfield(tracksAmpTotal,'amp2')
                        [tracksNA(:).amp2] = tracksAmpTotal.amp2;
                    end
                    if isfield(tracksAmpTotal,'ampTotal3')
                        [tracksNA(:).ampTotal3] = tracksAmpTotal.ampTotal3;
                    end
                end

                if ~isempty(forceReadProc)
                    forceReadObj = load(forceReadProc.outFilePaths_{1,p.ChannelIndex},'tracksForceMag'); % the later channel has the most information.
                    
                    tracksForceMag = forceReadObj.tracksForceMag;
                    idxTracksObj = load(forceReadProc.outFilePaths_{2,p.ChannelIndex},'idxTracks');
                    if ~isfield(idxTracksObj,'idxTracks')
                        idxTracksObj = load(forceReadProc.outFilePaths_{6,p.ChannelIndex},'idxTracks');
                    end
                    idxTracks = idxTracksObj.idxTracks;
                    %Apply if idSelected is applied
                    if ~isempty(idSelected)
                        idxTracks = find(idxTracks(idSelected));
                        tracksForceMag = tracksForceMag(idxTracks);
                    end
                
                    tracksNA = tracksNA(idxTracks);

                    if isfield(tracksForceMag,'forceMag')
                        [tracksNA(:).forceMag] = tracksForceMag.forceMag;
                    end
                    
                    if any(idxTracks)
                        output2 = idxTracks;
                    else
                        output2 = [];
                    end
                else
                    output2 = [];
                end
                toc
                output = tracksNA;          
                
            elseif strcmp(outputRequested,'idClass')
                %% Reading classes
                disp('Reading idsClassified ...')
                iFAPack = MD.getPackageIndex('FocalAdhesionPackage');
                FAPack=MD.packages_{iFAPack}; iForceRead=10;

                forceReadProc=FAPack.processes_{iForceRead};
                if ~isempty(forceReadProc)
                    idxTracksObj = load(forceReadProc.outFilePaths_{6,p.ChannelIndex},'idxTracks');
                    if ~isfield(idxTracksObj,'idxTracks')
                        idxTracksObj = load(forceReadProc.outFilePaths_{2,p.ChannelIndex},'idxTracks');
                        disp('Reading idxTracks from other source. OK!')
                    end
                    idxTracks = idxTracksObj.idxTracks;
                end
                try
                    try
                        iClaProc = MD.getProcessIndex('AdhesionClassificationProcess');
                        classProc = MD.getProcess(iClaProc);
                    catch
                        classProc = faPack.processes_{8};
                    end
                    % numChans = numel(p.ChannelIndex);
                    idsClassified = load(classProc.outFilePaths_{4,p.ChannelIndex});
                    idGroup1 = idsClassified.idGroup1;
                    idGroup2 = idsClassified.idGroup2;
                    idGroup3 = idsClassified.idGroup3;
                    idGroup4 = idsClassified.idGroup4;
                    idGroup5 = idsClassified.idGroup5;
                    idGroup6 = idsClassified.idGroup6;
                    idGroup7 = idsClassified.idGroup7;
                    idGroup8 = idsClassified.idGroup8;
                    idGroup9 = idsClassified.idGroup9;
                catch
                    disp('No Classified groups. Using no classification...')
                    % Potentially we can use the dedactically chosen classes here (e.g.
                    % shown in the code used for Tristan's movie analysis.
                    idGroup1 = []; idGroup2 = []; idGroup3 = []; idGroup4 = [];
                    idGroup5 = []; idGroup6 = []; idGroup7 = []; idGroup8 = []; idGroup9 = [];
                end    
                if ~isempty(forceReadProc)
                    idsClassified.idGroup1 = idGroup1(idxTracks);
                    idsClassified.idGroup2 = idGroup2(idxTracks);
                    idsClassified.idGroup3 = idGroup3(idxTracks);
                    idsClassified.idGroup4 = idGroup4(idxTracks);
                    idsClassified.idGroup5 = idGroup5(idxTracks);
                    idsClassified.idGroup6 = idGroup6(idxTracks);
                    idsClassified.idGroup7 = idGroup7(idxTracks);
                    idsClassified.idGroup8 = idGroup8(idxTracks);
                    idsClassified.idGroup9 = idGroup9(idxTracks);
                else
                    disp('Traction reading was not done. No further filtering...')
                end
                output = idsClassified;          
            end
        end
    end
    methods (Static)
        function name = getName()
            name = 'Initial-Rise Time Lag Calculation';
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieObject'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            adhAnalProc = owner.getProcess(owner.getProcessIndex('AdhesionAnalysisProcess'));
            pAnal=adhAnalProc.funParams_;
            
            ip.addOptional('ChannelIndex',pAnal.ChannelIndex,...
               @(x) all(owner.checkChanNum(x)));
            ip.parse(owner,varargin{:})
            
            % Set default parameters
            funParams.OutputDirectory = [ip.Results.outputDir filesep 'InitialRiseTimeLagCalculation'];
            funParams.ChannelIndex = ip.Results.ChannelIndex;
            funParams.mainSlave = 1;
        end
        
        function h = GUI()
            h = @initialRiseTimeLagCalculationProcessGUI;
        end
    end
end
