function [maxNMSInt,img] = GCAGetSteerableFilterScaleOverlaysMovie(movieData,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%
%   movieData - The MovieData object describing the movie, as created using
%   movieSelectorGUI.m
%
%   Parameter Structure Field Names:
%
% Generic Fields: (Input/Output Fields Needed for Wrapper)
%       ('OutputDirectory' -> Optional. A character
%       string specifying the directory to save the filoBranchInfo structure to.
%       If not input, the filoBranch will be saved in the same directory
%       as the movieData, in a sub-directory called
%       "IIIa_veilStem_refine_ActiveContours"
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
%
%  GCAveilStemRefineWithActiveContours specific functions
%
% OUTPUT: (see main function GCAveilStemRefineWithActiveContours

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
    'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'IIIa_veilStem_refine_ActiveContour'];

ip.addParameter('OutputDirectory',defaultOutDir,@(x) ischar(x));
ip.addParameter('ChannelIndex',1);
ip.addParameter('ProcessIndex',0);
ip.addParameter('frame',1);


% Specific
ip.addParameter('TSOverlays',true,@(x) islogical(x));


ip.addParameter('RidgeScalesAC',1:10,@(x) isvector(x)); % you can use the previous scale integration (0.214 to 2 um)
ip.addParameter('ThreshNMSResponseAC',75,@(x) isscalar(x));
ip.addParameter('MinRidgeScale',3, @(x) isscalar(x)); 
ip.addParameter('Iterations',100,@(x) isscalar(x)); 
ip.addParameter('SmoothingFactor',1,@(x) isscalar(x));                                          

ip.parse(varargin{:});
p = ip.Results;

%% Init:
nFrames = movieData.nFrames_;
nChan = numel(p.ChannelIndex);
channels = p.ChannelIndex;
imSize = movieData.imSize_;
ny = imSize(1);
nx = imSize(2);

%%

for iCh = 1:nChan
  
    display(['Refining Thin Veil/Stems Using Active Contours ' num2str(channels(iCh))]);
    % make final output dir where backboneInfo will be saved
    outDirC =  [ip.Results.OutputDirectory filesep 'Channel_' num2str(channels(iCh))];
    
    if ~isdir(outDirC)
        mkdir(outDirC);
    end
    
      
%% Get Restart Information 
    veilStemFile = [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep 'StepsToReconstruct'... 
        filesep 'III_veilStem_reconstruction' filesep 'Channel_' num2str(channels(iCh)) ... 
        filesep 'veilStem.mat'];

        load(veilStemFile) % load the file
       
       
        % check to make sure veilStemComplete 
%% Maybe have a check for outliers step ...

%% 

    veilStemMaskC = veilStem(ip.Results.frame).finalMask;    
    img = double(imread([movieData.getChannelPaths{channels(iCh)} filesep movieData.getImageFileNames{channels(iCh)}{ip.Results.frame} ])); 
    [ TSFigs,maxNMSInt,scaleMapInt] =  GCAGetSteerableFilterScaleOverlays(img,veilStemMaskC,p);    
    %veilStem(iFrame).finalMask = veilStemMaskRefineC;  
%     veilStem(iFrame).maxNMSInt  = maxNMSInt; 
%     veilStem(iFrame).scaleMapInt = scaleMapInt; 
    
 % make the directories for the figures if required. 
        for iFig = 1:length(TSFigs)
           if ~isdir([outDirC filesep TSFigs(iFig).name]); 
               mkdir([outDirC filesep TSFigs(iFig).name]); 
           end  
        end 
            type{1} = '.fig'; 
            type{2} = '.tif'; 
            
        if ~isempty(TSFigs)
            for iType = 1:numel(type)
            arrayfun(@(x) saveas(x.h,...
                [outDirC filesep x.name filesep num2str(ip.Results.frame,'%03d') type{iType}]),TSFigs);   
            end 
        end 
    
    
  save([outDirC filesep 'maxNMSInt.mat'],'maxNMSInt') ;  


end 
end 








