function [] = GCAVisualsMakeOverlaysWindows2Movie(movieData,varargin)
%GCAVisualsMakeOverlaysWindowMovie :

%% INPUT
%
%   movieData - The MovieData object describing the movie, as created using
%   movieSelectorGUI.m
%
%   Parameter Structure Field Names:
%
% Generic Fields: (Input/Output Fields Needed for Wrapper)
%       ('OutputDirectory' -> Optional. A character
%       string specifying the directory to save the Visualization Output
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) on which to the windows will be overlaid
%
%       This index corresponds to the channel's location in the array
%       movieData.channels_. If not input, all channels will be analyzed
%
%       ('ProcessIndex' -> Positive integer scalar or vector)
%       This parameter specifies the output process to use for the overlay
%       This allows a previously processed image (ie shade corrected or
%       background subtracted to potentially be input as the base of the overaly).
%       If not input, the raw images will be used for the overlay
%
%% Check Input
% REQUIRED
% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%% Check Input


if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
defaultOutDir = [movieData.outputDirectory_ filesep ...
    'VisualizationOverlays' filesep 'WholeNeurite' filesep ...
    'VeilWindows'];



ip.addParameter('OutputDirectory',defaultOutDir,@(x) ischar(x));
ip.addParameter('ChannelIndex',1);
ip.addParameter('ProcessIndex',0);
ip.addParameter('StartFrame','auto');
ip.addParameter('EndFrame','auto');

ip.parse(varargin{:});
p = ip.Results;

%% Init:
nFrames = movieData.nFrames_;
nChan = numel(p.ChannelIndex);
channels = p.ChannelIndex;
imSize = movieData.imSize_;
ny = imSize(2);
nx = imSize(1);

%% Loop for each channel
for iCh = 1:nChan
    
    display(['Making Window Overlays For ' num2str(channels(iCh))]);
    
    %% Make the OutputDirectory
    saveDir = [ip.Results.OutputDirectory filesep 'Channel_' num2str(iCh)];
    if ~isdir(saveDir) ;
        mkdir(saveDir);
    end
    %% Load the Windows
    idxWind = find(cellfun(@(x) sum(strcmpi(x.name_,'Windowing')),movieData.processes_));
    
    if ~isempty(idxWind)
        % load the three protrusion process associated cell structures,
        % 'normals','protrusion','smoothedEdge'
        windowFiles = searchFiles('.mat',[],movieData.processes_{idxWind(end)}.outFilePaths_,0);
        % If for some reason there is more than one protrusion process
        % associated with the data, tell the user that your are using the
        % most recently run
        if length(idxWind) >1
            
            display(['There is more than one veil windowing process associated with ' ...
                'this movie: the most recent will be plotted']);
        end % if length(protS)
        
        
        check(1)= 1;
    else % no veil protrusion vectors were found
        warning(['Windowing not found for this movie :'...
            'Run Windowing Package and Try Again']);
        check(1) = 0 ;
    end
    %% Load Protrusion Vectors if user wants to plot
    idxProt = find(cellfun(@(x) sum(strcmpi(x.name_,'Protrusion')),movieData.processes_));
    if ~isempty(idxProt)
        % load the three protrusion process associated cell structures,
        % 'normals','protrusion','smoothedEdge'
        load(movieData.processes_{idxProt(end)}.outFilePaths_);
        % If for some reason there is more than one protrusion process
        % associated with the data, tell the user that your are using the
        % most recently run
        if length(idxProt) >1
            display(['There was more than one veil protrusion process associated with ' ...
                'this movie: '...
                'will be performed using the most recently run protrusion process']);
            % if length(protS)
        else % no veil protrusion vectors were found
            display(['The veil protrusion process was not found for this movie :'...
                'Please Run the Protrusion Vectors']);
        end
    end % idxProt
    %% Load the edge velocity analysis
    load([movieData.outputDirectory_ filesep 'EdgeVelocityQuantification' filesep 'EdgeMotion.mat']);
    
    windowInfoProt = analysisResults.protrusionAnalysis.windows;
    %% MOVIE WRAPPER
    for iFrame = 1:nFrames-1
        % load img
        img = double(imread([movieData.getChannelPaths{1} filesep movieData.getImageFileNames{iCh}{iFrame}]));
        
        % load windows
        load([windowFiles{iFrame,2} filesep windowFiles{iFrame,1}]);
        
        normalsC = normals{iFrame};
        protrusionC = protrusion{iFrame};
        edgeC = smoothedEdge{iFrame};
        
        setFigure(nx,ny,'on');
        imshow(-img,[]) ;
        hold on
        
        if ~isempty(windows{1}) ; % not sure why but it appears hunter leaves this empty sometimes ?
            
            windowFirst = windows{1}(1) ;
        else
            windowFirst = windows{2}(1);
        end
        gcaPlotWindows(windowFirst,{'y','FaceAlpha',1},0,'bandMax',1);
        
       
        %display(['Plotting Frame' num2str(iFrame)]);
        
        gcaPlotWindows(windows,{'k','FaceAlpha',0},5,'bandMax',1);
        % find all windows that were protrusion events that passed the
        % time series analysis criteria in the current frame
        nWinds = length(analysisResults.protrusionAnalysis.windows(:));
        nWindsC = numel(windows); 
        idxBlackOut = arrayfun(@(x) ~isempty(find(find(isnan(analysisResults.data.procEdgeMotion(x,:)))==iFrame)),1:nWindsC); 
        windowsOut = windows(idxBlackOut); 
        if ~isempty(windowsOut) 
         gcaPlotWindows(windowsOut,{'k','FaceAlpha',1},0,'bandMax',1); 
         end
        idxProtC =  arrayfun(@(x) ~isempty(find(vertcat(analysisResults.protrusionAnalysis.windows(x).blockOut{:})==iFrame)),1:nWinds);
        windowsProt =  windows(find(idxProtC));
        
        %         % for each window
        %         for iWind = 1:nWinds
        %         %windowInfoProt = analysisResults.protrusionAnalysis.windows(1).blockOut;
        %         blockTimes = vertcat(analysisResults.protrusionAnalysis.windows(iWind).blockOut{:});
        %
        gcaPlotWindows(windowsProt,{'r','FaceAlpha',1},0,'bandMax',1);
        %GCAVisualsMakeOverlaysVeilStemWindows(windowsProt,normalsC,protrusionC,edgeC);
        hold on
        idxRetC = arrayfun(@(x) ~isempty(find(vertcat(analysisResults.retractionAnalysis.windows(x).blockOut{:}) == iFrame)),1:nWinds);
        windowsRet = windows(find(idxRetC));
        gcaPlotWindows(windowsRet,{'b','FaceAlpha',1},0,'bandMax',1);
        %GCAVisualsMakeOverlaysVeilStemWindows(windowsRet,normalsC,protrusionC,edgeC);
        saveas(gcf,[saveDir filesep   num2str(iFrame,'%03d') '.png']);
        saveas(gcf,[saveDir filesep   num2str(iFrame,'%03d') '.eps'],'psc2');
        saveas(gcf,[saveDir filesep   num2str(iFrame,'%03d') '.fig']);
        close gcf
    end % iFrame
    
end % nCh

%% Create Movies if user desires.












%%


