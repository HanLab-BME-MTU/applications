function [ output_args ] = GCAfindVeilStemLongestPathMovie( movieData,varargin)
%%% GCAfindVeilStemLongestPathMovie
% STEP IV in GCA Segmentation
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%
%   movieData (REQUIRED)  - The MovieData object describing the movie, as created using
%   movieSelectorGUI.m
%
% Generic Fields:
%       ('Input/Output Fields Needed for Wrapper' -> Possible Values)
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the backboneInfo structure to.
%       If not input, the backboneInfo will be saved in the same directory
%       as the movieData, in a sub-directory called
%       "neurite_orientation"
%       % PERSONAL NOTE : You might need
%       to truncate these names for the windows paths. %
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) on which to perform the neurite orientation
%       estimation.
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
%
%       ('Restart'  -> structure with fields = 'auto';   % 'auto'
%       paramsIn.restart.endFrame = 'auto'; % 'auto' or number, default auto.
%       paramsIn.plots = 1;
%
%  SPECIFIC INPUT
%  See GCAfindVeilStemLongestPath.m for details.
%% %% INPUTPARSER
% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
defaultOutDir = [movieData.outputDirectory_ filesep...
    'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'IV_veilStem_length'];

ip.addParameter('OutputDirectory',defaultOutDir,@(x) ischar(x));
ip.addParameter('ChannelIndex',1);
ip.addParameter('ProcessIndex',0);
ip.addParameter('StartFrame','auto');
ip.addParameter('EndFrame','auto');

% Specific
ip.addParameter('TSOverlays',true,@(x) islogical(x));
%ip.addParameter('TSMovie',false,@(x) islogical(x));
ip.addParameter('neuriteElongTS_medFiltWindSize',10,@(x) isscalar(x)); % see gcaFindOutliersFromMedFilt.m
ip.addParameter('neuriteElongTSOutlier_outlierDef_k',3,@(x) isscalar(x)); % see gcaFindOutliersFromMedFilt.m

ip.parse(varargin{:});
p = ip.Results;
%% Initiate
nFrames = movieData.nFrames_;
nChan = numel(p.ChannelIndex);
imSize = movieData.imSize_;

%% Start Channel Wrapper

for iCh = 1:nChan
    
    display(['Finding VeilStem Longest Path for Channel ' num2str(iCh)]);
    %% Get Start and End Frames Based on Restart Choice
    
    % make final output dir where backboneInfo will be saved
    outDirC =  [ip.Results.OutputDirectory filesep 'Channel_' num2str(iCh)];
    
    if ~isdir(outDirC)
        mkdir(outDirC);
    end
    
    if p.ProcessIndex == 0
        imgDir =  movieData.channels_(p.ChannelIndex).channelPath_; % currently mainly needed if you do the local thresholding/ otherwise just overlay
    else
        imgDir = movieData.proccesses_(p.ProcessIndex).outFilePaths_{p.ChannelIndex};
    end
    
    % collect images and initiate
    [listOfImages] = searchFiles('.tif',[],imgDir ,0);
    
    veilStemFile = [outDirC filesep 'veilStem.mat'];
    
    %% Get Restart Information
    % longPathFile = [outDirC filesep 'veilStemInfo.mat'];
    
    % If veilStem.mat  already exists load it
    if  exist(veilStemFile,'file')==2;
        display('Previous Measurments Detected'); 
        load(veilStemFile) % load the file
        load([ip.Results.OutputDirectory filesep 'Channel_' num2str(iCh) filesep 'neuriteLength.mat']);
        display('Loading Previously Run Longest Path VeilStem Measurements');
        if strcmpi(ip.Results.StartFrame,'auto')
            startFrame = find(arrayfun(@(x) isempty(x.neuriteLongPathIndices),veilStem),1,'first')-1;
            if isempty(startFrame)
                startFrame = length(veilStem);
            end
            
            
            if startFrame == 0
                startFrame = 1; % reset to 1;
            end
            display(['Auto Start: Starting Veil/Stem Reconstruction at Frame ' num2str(startFrame)]);
        else
            startFrame = ip.Results.StartFrame; % use user input
            display(['Manual Start: Starting Veil/Stem Reconstruction at Frame ' num2str(startFrame)]);
        end
    else % if doesn't exist
        
        startFrame = 1;
        
        display('No Veil/Stem Folder Found: Creating and Starting at Frame 1');
        %%    Load Veil Stem Information from Step III 
         load([ [movieData.outputDirectory_ filesep...
        'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'III_veilStem_reconstruction'] filesep 'Channel_' num2str(iCh)...
        filesep 'veilStem.mat']);
        
    end % exist(orientFile,'file') == 2
    
    % if restarting add saved parameters
    %     if startFrame ~= 1
    %         load([outDirC filesep 'params.mat']);
    %     end
    
    
    if strcmpi(ip.Results.EndFrame,'auto');
        endFrame = nFrames;
        display(['Auto End: Finding Veil/Stem From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    else
        endFrame = ip.Results.EndFrame;
        display(['Manual End: Finding Veil/Stem From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    end
   
    %% Main Function:
    for iFrame = startFrame:endFrame
        
        veilStemMaskC = veilStem(iFrame).finalMask;
        idxEnterNeuriteC = veilStem(iFrame).idxEnterNeurite;
        
        [neuriteLengthC,longPathLinIndC,EPLongPathC,~] = GCAfindVeilStemLongestPath(veilStemMaskC,idxEnterNeuriteC);
       
        neuriteLength(iFrame,1) = neuriteLengthC;
        
        veilStem(iFrame).neuriteLongPathIndices = longPathLinIndC;
        veilStem(iFrame).endPointLeadingProt = EPLongPathC;
        
        save([ip.Results.OutputDirectory filesep 'Channel_' num2str(iCh) filesep 'veilStem.mat'],'veilStem','-v7.3') ;
        save([ip.Results.OutputDirectory filesep 'Channel_' num2str(iCh) filesep 'neuriteLength.mat'], 'neuriteLength');
        display(['Longest Path Found for Frame ' num2str(iFrame)]);
    end
    
    
    % After complete put into the neurite outgrowth measurement file
    measDir =  [ movieData.outputDirectory_ filesep 'MEASUREMENT_EXTRACTION' filesep 'GlobalFunctional' filesep ...
        'neurite_outgrowth_measurements'];
    if ~isdir(measDir)
        mkdir(measDir) ;
    end
    timeStamp = clock;
    save([measDir filesep 'neuriteLengthOutput.mat'],'neuriteLength','timeStamp');
    %% Peform the median filtering to test for veil/stem outliers
    %gcaPerformMedianFiltering(
    % for each outlier put flag
    if isfield(veilStem,'flagOutlier');
        veilStem = rmfield(veilStem,'flagOutlier');
        display('A .flagOutlier field was previously found and was overwritten');
    end
    
    
    neuriteLength = neuriteLength - neuriteLength(1);
    
    [~, outlierIdx, TSFig] = gcaFindOutliersFromMedFilt(neuriteLength,ip.Results.neuriteElongTS_medFiltWindSize,ip.Results.neuriteElongTSOutlier_outlierDef_k); % add these parameters
    
    fileType{1} = '.fig';
    fileType{2} = '.png';
    for ifile = 1:numel(fileType)
        if ip.Results.TSOverlays   == true;
            saveas(TSFig.h,[ip.Results.OutputDirectory filesep 'Channel_' num2str(iCh) filesep TSFig.name fileType{ifile}]);
            
            saveas(TSFig.h,[measDir filesep TSFig.name fileType{ifile} ]);
        end
    end
    
    
    for i = 1:length(outlierIdx)
        frame = outlierIdx(i);
        veilStem(frame).flagOutlier = true;
    end
    
    save([ip.Results.OutputDirectory filesep 'Channel_' num2str(iCh) filesep 'veilStem.mat'],'veilStem','-v7.3') ;

    save([ip.Results.OutputDirectory filesep 'Channel_' num2str(iCh) filesep 'params.mat'],'p');  
 %% Extra: Extract Thickness (may do this eventually in a different step  
 GCAAnalysisExtract_veilStemThicknessMovie(movieData,veilStem); 
    
end % for iCh
end % function
