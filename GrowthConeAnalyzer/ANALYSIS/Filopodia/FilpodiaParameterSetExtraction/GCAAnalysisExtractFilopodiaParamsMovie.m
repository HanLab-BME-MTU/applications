function [ output_args ] = GCAAnalysisExtractFilopodiaParamsMovie(MD)
%GCAAnalysisExtractFilopodiaParamsMovie
%   
%
 

%% Load the segmenatation analysis
display('Please Be Patient This File Takes a While to Load...');
load([MD.outputDirectory_ filesep 'filopodia_fits' filesep 'Filopodia_Fits_Channel_1' filesep 'analInfoTestSave.mat']);
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
    newFiloDir  = [MD.outputDirectory_ filesep 'PARAMETER_EXTRACTION' ...
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
    display(['Fits for ' MD.outputDirectory_ 'were not fully run (starting at frame '  num2str(startProblemFrame) ...
        ' ): SKIPPING ' ]); 
       
    
end % sum missing 
end 









