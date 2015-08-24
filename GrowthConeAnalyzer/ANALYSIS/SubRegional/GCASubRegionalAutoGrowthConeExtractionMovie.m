function [ output_args ] = GCAGetGrowthConeSubregionsMovie(movieData,varargin)
%GCASubregionalAnalysisMovie : Function to partition the GC automatically
%into subregions- module first introduced is currently based on a pre-set 
% value from the tip of the leading protrusion 
% Requires that the skeleton for the original GC has already been run. 
% NOTE: % Currently assumes you want to make a mask from the smoothed Edge after protrusion! 
%% 
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
%       If not input, 
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) on which to perform the filopodia
%       reconstruct
%       This index corresponds to the channel's location in the array
%       movieData.channels_. If not input, all channels will be analyzed
%
%
%
%  GCAGetGrowthConeSubregions specific functions
%
% OUTPUT: (see main function GCAGetGrowthConeSubregions- for details):
%      
%% InputParser
%% INPUTPARSER
% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
defaultOutDir = [movieData.outputDirectory_ filesep...
    'SegmentationPackage'  filesep 'GCASubRegions'];

 defaultInDir = [movieData.outputDirectory_ filesep ... 
    'SegmentationPackage' filesep 'StepsToReconstruct'  filesep 'IV_veilStem_length' ]; 



ip.addParameter('OutputDirectory',defaultOutDir,@(x) ischar(x));
ip.addParameter('InputDirectory', defaultInDir,@(x) ischar(x)); 
ip.addParameter('ChannelIndex',1);
ip.addParameter('ProcessIndex',0);
ip.addParameter('StartFrame','auto');
ip.addParameter('EndFrame','auto');


% Specific
% PARAMETERS
ip.addParameter('TSOverlays',true); 

% Distance from leading edge protrusion (in um) 
ip.addParameter('maskFromSmoothedEdge',true); 
ip.addParameter('distFromLeadProt',20); 
ip.addParameter('GCFinder',[]); 
ip.addParameter('angle',90); % in degrees- relative to the local direction of the skeleton at the idxPt. 

ip.parse(varargin{:});
params = ip.Results; 
%% Initiate 
%% Init:
nFrames = movieData.nFrames_;
nChan = numel(params.ChannelIndex);
channels = params.ChannelIndex;
imSize = movieData.imSize_;
ny = imSize(1);
nx = imSize(2);
% load the veilStem structure : note this is a very large structure and
% contains the data for the entire movies - can in the future potentially 
% partition data per frame 

if ip.Results.maskFromSmoothedEdge == true 
    % get protrusion process 
    idxProt =  cellfun (@(x) strcmpi(x.name_,'Protrusion'),movieData.processes_); 
    load(movieData.processes_{idxProt}.outFilePaths_);    
end 
 load([ip.Results.InputDirectory filesep 'Channel_' num2str(channels(1)) filesep  'veilStem.mat']);

 %% Make the directories 
    subNames{1} = 'GC';
    subNames{2} = 'Stem';
    
    subRegDir = ip.Results.OutputDirectory; 
    % make the directories for masks
    arrayfun(@(i) mkClrDir([subRegDir filesep subNames{i} filesep 'masks']),1:2);
 
     if ip.Results.TSOverlays == true 
     mkClrDir([subRegDir filesep 'Overlays']); 
     filenames = movieData.getImageFileNames{1}; 
     end
    
%% Loop 
for iFrame = 1:nFrames
    % load the longest path linear coordinates for the given frame
    nLPC = veilStem(iFrame).neuriteLongPathIndices; 
    if ip.Results.maskFromSmoothedEdge == true 
        idx = sub2ind([ny,nx],round(smoothedEdge{iFrame}(:,2)),round(smoothedEdge{iFrame}(:,1))); 
        veilStemC = zeros([ny,nx]); 
        veilStemC(idx) = 1;
        veilStemC = logical(veilStemC); 
        veilStemC = bwmorph(veilStemC,'bridge'); 
        veilStemC = imfill(veilStemC,'holes');  
    else 
    % load the veilStem mask 
    veilStemC = veilStem(iFrame).finalMask; 
    end 
    % Cut the longest path indices by the user specified value
    [~,measIndices] = calculateDistance(nLPC,[ny,nx],'distCutOff',ip.Results.distFromLeadProt); 
    
    [ subRois,xVect,yVect,pixelInfoGC,defineGCPlot ] =  GCASubRegionalAutoGrowthConeExtraction(measIndices,veilStemC,'angle',ip.Results.angle); 
    
    % Write the masks to the appropriate directory.
           arrayfun(@(i) imwrite(subRois(:,:,i),[subRegDir filesep subNames{i} filesep 'masks' filesep 'mask' num2str(iFrame,'%03d') '.tif']),1:2);
if ip.Results.TSOverlays == true 
    img = double(imread([movieData.getChannelPaths{1} filesep filenames{iFrame}]));
    setFigure(nx,ny,'off'); 
    imshow(-img,[]) 
    hold on 
    roiYX = arrayfun(@(i) bwboundaries(subRois(:,:,i)),1:2,'uniformoutput',0);
    color(1) = 'm'; 
    color(2) = 'b'; 
    
    measIndicesMask = zeros(ny,nx); 
    measIndicesMask(measIndices) = 1;
    arrayfun(@(i) cellfun(@(x) plot(x(:,2),x(:,1),color(i)),roiYX{i}),1:2); 
    hold on 
    spy(measIndicesMask,'k'); 
    text(5,5,' Within 20 um from tip of leading protrusion','color','m'); 
    text(5,20,'> 20 um from tip ','color','b'); 
    saveas(gcf, [ subRegDir filesep 'Overlays' filesep num2str(iFrame,'%03d') '.png']); 
    close gcf 
end 
    

    
end 
    
end 


