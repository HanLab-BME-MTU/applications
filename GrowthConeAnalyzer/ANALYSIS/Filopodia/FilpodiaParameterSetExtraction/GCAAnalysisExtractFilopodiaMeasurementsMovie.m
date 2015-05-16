function [ output_args ] = GCAAnalysisExtractFilopodiaMeasurementsMovie(movieData,varargin)
%GCAAnalysisExtractFilopodiaParamsMovie
%   This function makes the default filopodia filter sets and extracts all
%   default filopodia parameters
%% INPUT
%
%   movieData - The MovieData object describing the movie, as created using
%   movieSelectorGUI.m
%
%   Parameter Structure Field Names:
%
% Generic Fields: (Input/Output Fields Needed for Wrapper)
%       ('OutputDirectory' -> Optional. A character
%       string specifying the directory to save the filopodia measurements to.
%       If not input, the output will be saved in the same directory
%       as the movieData, in a sub-directory called
%       'GrowthConeAnalyzer' filesep 'GCAMeasurementExtraction' filesep
%       'WholeNeurite' 
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) on which to perform the filopodia
%       reconstruct
%       This index corresponds to the channel's location in the array
%       movieData.channels_. If not input, all channels will be analyzed
%
%       ('ProcessIndex' -> Positive integer scalar or vector)
%       This parameter specifies the output process to use for performing the
%       estimation
%       This allows a previously processed image (ie shade corrected or
%       background subtracted to potentially be input). If not input, the
%       backbone information will be calculated from the channels (ie raw
%       images)
%% Check input
% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
defaultOutDir = [movieData.outputDirectory_ filesep...
   'GCAMeasurementExtraction' filesep 'WholeNeurite' ];

ip.addParameter('OutputDirectory',defaultOutDir,@(x) ischar(x));
ip.addParameter('ChannelIndex',1);
ip.addParameter('ProcessIndex',0);

ip.parse(varargin{:});
p = ip.Results;
%% Init:
nFrames = movieData.nFrames_;
nChan = numel(p.ChannelIndex);
channels = p.ChannelIndex;
%% Loop for each channel
for iCh = 1:nChan
    
    %% Load the segmenatation analysis
    display('Please Be Patient This File Takes a While to Load...');
    load([movieData.outputDirectory_ filesep 'filopodia_fits' filesep 'Filopodia_Fits_Channel_' num2str(channels(iCh)) filesep 'analInfoTestSave.mat']);
    display('Finished Loading');
    
    
    
    % Check to make sure everything run completely eventually will just look at
    % the movieData process.
    for iFrame = 1:numel(analInfo) -1
        filoInfo = analInfo(iFrame).filoInfo;
        % arrayfun(@(x)
        % test if the field associated with endpoint coordinate is empty
        % 1 if missing fits 0 if not
        missingFits(iFrame,1)= sum(arrayfun(@(x) isempty(x.Ext_endpointCoordFitPix),filoInfo));
    end
    if sum(missingFits) ==0  % continue all is ok
        
        
        
        %% Set Up 'ConnectToVeil_LengthInt' Filopodia Filter Analysis Structure
        % Each analInput Structure specifies the filo filter type,
        % the param functions that need to be called, and the descriptive parameter name
        
        analInput(1).filterType = 'ConnectToVeil_LengthInt'; %
        
        % Length To Veil
        analInput(1).paramFunc{1} = 'filoLength'; % % function ID
        analInput(1).paramName{1} = 'filoLengthToVeil'; % paramName-
        analInput(1).paramInput{1} = 'Ext_';
        
        % Length Embedded
        analInput(1).paramFunc{2} = 'filoLength';
        analInput(1).paramName{2} = 'filoLengthEmbedded';
        analInput(1).paramInput{2} = 'Int_';
        
        % Total Length Actin Bundle
        analInput(1).paramFunc{3} = 'filoLength';
        analInput(1).paramName{3} = 'filoLengthFullActinBundle';
        analInput(1).paramInput{3} = 'Tot';
        
        % Intensity To Veil
        analInput(1).paramFunc{4} = 'filoAvgIntensity'; % function ID
        analInput(1).paramName{4} = 'filoIntensityToVeil'; % paramName for output
        analInput(1).paramInput{4} = 'Ext'; % other information for function
        
        % Intensity Embed
        analInput(1).paramFunc{5} = 'filoAvgIntensity';
        analInput(1).paramName{5} = 'filoIntensityEmbedded';
        analInput(1).paramInput{5}  = 'Int';
        
        
        
        %% Set up ConnectToVeil_DensityOrient Filter Type
        analInput(2).filterType = 'ConnectToVeil_DensityOrient';
        %
        % % Orientation
        analInput(2).paramFunc{1} = 'filoOrient';
        analInput(2).paramName{1} = 'filoOrient';
        analInput(2).paramInput{1} = [];
        %
        % % Density
        analInput(2).paramFunc{2} = 'filoDensityAlongVeil';
        analInput(2).paramName{2} = 'filoDensityAlongVeil';
        analInput(2).paramInput{2} = [];
        %% Wrap through for each analysis type
        for iAnalType = 1:length(analInput);
            
            % get the filopodia filter for analInput
            [filoFilterSet,filterParams] = GCACreateFilopodiaFilterSet(analInfo,analInput(iAnalType).filterType);
            
            
            % Make the Respective Folders
            newFiloDir  = [ ip.Results.OutputDirectory filesep 'Channel' num2str(channels(iCh))...
                filesep 'Descriptor' filesep 'Filopodia' filesep analInput(iAnalType).filterType];
            if isdir(newFiloDir)
                display([newFiloDir ' Found: SKIPPING']) ; % for now skip but can make an option to re-write
            else
                mkdir(newFiloDir);
            end
            
            % save the filter set used
            save([newFiloDir filesep 'filoFilterSet.mat'],'filoFilterSet','filterParams');
            
            % perform the various assocatied parameter extraction by calling the
            % associated function
            nParams  = numel(analInput(iAnalType).paramFunc);
            
            for iParamExtract = 1:nParams
                
                % Call the extraction function : all of these extraction functions have names starting with
                % GCAAnalysisExtract_ '
                % some are very short however I found it useful to keep each
                % extraction as a separate function- all output will be a cell of
                % iFrames holding the distribution of values
                % for that parameter
                paramFuncC = str2func(['GCAAnalysisExtract_' analInput(iAnalType).paramFunc{iParamExtract}]);
                inputC =  analInput(iAnalType).paramInput{iParamExtract};
                if ~isempty(inputC)
                    paramC =  paramFuncC(analInfo,filoFilterSet, inputC);
                else
                    paramC = paramFuncC(analInfo,filoFilterSet);
                end
                % Make the output directory: Each extraction has its own folder for now so have a place to
                % save associated plots and movies- will collect in a later step
                analOutputDir = [newFiloDir filesep  analInput(iAnalType).paramFunc{iParamExtract}];
                if ~isdir(analOutputDir)
                    mkdir(analOutputDir)
                end
                
                % save the parameter extraction
                nameParam = ['param_'  analInput(iAnalType).paramName{iParamExtract}];
                timeStamp = clock;
                save([analOutputDir filesep nameParam '.mat' ],'paramC','timeStamp');
                
            end % iParamExtract % iterate over all params using the same filter set
            
            
        end % iAnal % iterate over all analysis groups with the same filter set
        
    else % some of the fits have not been run
        startProblemFrame = find(missingFits);
        display(['Fits for ' movieData.outputDirectory_ 'were not fully run (starting at frame '  num2str(startProblemFrame) ...
            ' ): SKIPPING ' ]);
        
        
    end % sum missing
end  % for iCh
end









