function [ output_args ] =GCASubRegionalAutoGrowthConeExtractionMovie(movieData,varargin)
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
ip.addParameter('overwrite',false);

% Overlay Options 
ip.addParameter('TSOverlays',true);
ip.addParameter('ShowFig','off');  % flag to have figure pop up on screen
ip.addParameter('CirclesBeginEnd',false);
ip.addParameter('plotScaleBar',false)

% Specific Movie Function 
ip.addParameter('maskFromSmoothedEdge',true);
ip.addParameter('distFromLeadProt',[]); % Distance from leading edge protrusion (in um)
ip.addParameter('maxCMapValue',10); 

% Specific : GCASubRegionalAutoGrowthConeExtraction
ip.addParameter('GCFinder',true);
ip.addParameter('GCFindNeckCutOff',2,@isscalar); % the thicknes cut-off for the veil stem. 
ip.addParameter('GCFindMinLength',5,@isscalar); % the smallest allowed length of a growth cone
ip.addParameter('GCFindThickPt',2.5,@isscalar); % assumes the local max is at least this thick

ip.addParameter('angle',90); % in degrees- relative to the local direction of the skeleton at the idxPt.




ip.parse(varargin{:});
params = ip.Results;
%% Initiate

run = 0;
nFrames = movieData.nFrames_;
nChan = numel(params.ChannelIndex);
channels = params.ChannelIndex;
imSize = movieData.imSize_;
ny = imSize(1);
nx = imSize(2);
% FIX 
if strcmpi(ip.Results.StartFrame,'auto');
    startFrame = 1;
else
    startFrame = ip.Results.StartFrame;
end
if strcmpi(ip.Results.EndFrame,'auto');
    endFrame = nFrames;
else
    endFrame = ip.Results.EndFrame;
end

% load the veilStem structure : note this is a very large structure and
% contains the data for the entire movies - can in the future potentially
% partition data per frame

if ip.Results.maskFromSmoothedEdge == true
    % get protrusion process
    idxProt =  cellfun (@(x) strcmpi(x.name_,'Protrusion'),movieData.processes_);
    load(movieData.processes_{idxProt}.outFilePaths_);
end
load([ip.Results.InputDirectory filesep 'Channel_' num2str(channels(1)) filesep  'veilStem.mat']);

%% Check distFromLeadProt to make sure it is not 
% Check if in any frame the neurite length profile was sampled at less than 
  % distFromLeadProt: if so give an error message 
  
   % get all of the neurite length metrics for each frame 
         allDist =  arrayfun(@(x) calculateDistance(veilStem(x).neuriteLongPathIndices,[ny,nx]),1:nFrames);
         allDist = allDist'; 
        
         if ~isempty(ip.Results.distFromLeadProt)
              toTest = ip.Results.distFromLeadProt; 
             if sum(allDist < toTest ) >0;
                 problemFrames = find(allDist<toTest);
                 shortestDist = min(allDist);
                 %problemFrames = cellfun(@(x) num2str(x),problemFrames,'uniformoutput',0);
                 errorMessage =['The distFromLeadingProt parameter needs '  ...
                     'to be truncated less than or equal to ' num2str(floor(shortestDist))  ...
                     ': or the movie needs to be recropped '] ;
                 save([ip.Results.OutputDirectory filesep 'ErrorReport.mat'],'errorMessage','problemFrames');
                 error(errorMessage);      
             end
         end

%% Make the directories
subNames{1} = 'GC';
subNames{2} = 'Stem';

% Tentative portion for file names for testing - include if GC Finder on
% and params 
if ~ip.Results.GCFinder
forName = ''; 
else 

  pn{1} = ['_NeckW_' num2str(ip.Results.GCFindNeckCutOff)]  ;   
  pn{2} = ['_MinL_' num2str(ip.Results.GCFindMinLength) ] ; 
  pn{3} = ['_ThickPt_' num2str(ip.Results.GCFindThickPt) ] ; 
  
  pn = cellfun(@(x)  strrep(x,'.','pt'),pn,'uniformoutput',0); 
    
forName = ['GCFinder' horzcat(pn{:}) ]; 
end 

subRegDir = [ip.Results.OutputDirectory forName];

if ~isdir(subRegDir)
    run = 1 ;
end

if ip.Results.overwrite
    run = 1;
end






if run == 1
    % make the directories for masks
    arrayfun(@(i) mkClrDir([subRegDir filesep subNames{i} filesep 'masks']),1:2);
    
    if ip.Results.TSOverlays == true
        visDir = [subRegDir filesep 'Visualization'];
        mkClrDir(visDir);
        filenames = movieData.getImageFileNames{1};
    end
     
   
     if isempty(ip.Results.distFromLeadProt)
      
         maxDistProfile = ceil(max(allDist)); 

     else 
         maxDistProfile = ip.Results.distFromLeadProt; 
     end
    
    %% Loop
    for iFrame =startFrame:endFrame
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
        % Cut the longest path indices by the user specified value if desired
        if ~isempty(ip.Results.distFromLeadProt)
            [~,measIndices] = calculateDistance(nLPC,[ny,nx],'distCutOff',ip.Results.distFromLeadProt);
        else
            measIndices = nLPC;
        end
        [ subRois,xVect,yVect,pixelInfoGC,defineGCPlot,idxCMap,cmap,distTransInMicGC,GCDetect ] =  GCASubRegionalAutoGrowthConeExtraction(measIndices,veilStemC,... 
            'maxDistProfile', maxDistProfile, params);
        
        
        % Write the masks to the appropriate directory.
        arrayfun(@(i) imwrite(subRois(:,:,i),[subRegDir filesep subNames{i} filesep 'masks' filesep 'mask' num2str(iFrame,'%03d') '.tif']),1:2);
        if ip.Results.TSOverlays == true
            img = double(imread([movieData.getChannelPaths{1} filesep filenames{iFrame}]));
            
            % Make the main overlays (Currently only for the GC Detector need to
            % remake for the vector partitioning)
            cDir = [visDir filesep 'OverlayPartitions'];
            if ~isdir(cDir)
                mkdir(cDir)
            end
            setFigure(nx,ny,ip.Results.ShowFig);
            imshow(-img,[])
            hold on
            roiYX = arrayfun(@(i) bwboundaries(subRois(:,:,i)),1:2,'uniformoutput',0);
            color(1) = 'r';
            color(2) = 'k';
            %     measIndicesMask(measIndices) = 1;
            arrayfun(@(i) cellfun(@(x) plot(x(:,2),x(:,1),color(i)),roiYX{i}),1:2);
            saveas(gcf, [ cDir filesep num2str(iFrame,'%03d') '.png']);
            saveas(gcf, [ cDir filesep num2str(iFrame,'%03d') '.fig']);
            close gcf
            %% Save the Plot for the distance transform definition.
            if ~isempty(defineGCPlot)
                if ~isdir([visDir filesep 'GCDefinition']);
                    mkdir([visDir filesep 'GCDefinition']);
                end
                saveas(gcf,[visDir filesep 'GCDefinition' filesep num2str(iFrame,'%03d') '.png']);
                saveas(gcf,[visDir filesep 'GCDefinition' filesep num2str(iFrame,'%03d') '.fig']);
                close gcf
                
                
                %% Make and Save the plot for the
                
                cDir = [visDir filesep 'DistTransOverlay'];
                if ~isdir(cDir)
                    mkdir(cDir)
                end
                
                setFigure(nx,ny,ip.Results.ShowFig);
                imshow(-img,[])
                hold on
                roiYX = arrayfun(@(i) bwboundaries(subRois(:,:,i)),1:2,'uniformoutput',0);
                color(1) = 'k';
                color(2) = 'k';
                [plotY,plotX] = ind2sub([ny,nx],measIndices);
                
                for iColor = 1:length(idxCMap)
                    if ~isempty(measIndices(idxCMap(:,1)==iColor));
                        scatter(plotX(idxCMap(:,1)==iColor),plotY(idxCMap(:,1)==iColor),50,cmap(iColor,:),'filled');
                    end
                end
                idxMax =  find(distTransInMicGC==max(distTransInMicGC));
                
                distAtMax = distTransInMicGC(idxMax);
                distAtMaxPix = distAtMax./0.216; % convert back to pix
                [yMax,xMax] =  ind2sub([ny,nx],pixelInfoGC(idxMax));
                yMax = yMax(1);
                xMax = xMax(1);
                distAtMaxPix = distAtMaxPix(1);
                idxMax = idxMax(1);
                
                gcaCircles(xMax,yMax,distAtMaxPix,'facecolor','none','edgeColor',cmap(idxCMap(idxMax,1),:));
                
                if ip.Results.CirclesBeginEnd
                    dist1 = distTransInMicGC(1)./0.216;
                    gcaCircles(plotX(1),plotY(1),dist1,'facecolor','none','edgeColor',cmap(idxCMap(1,1),:));
                   
                    distEnd = distTransInMicGC(end)./0.216;
                    gcaCircles(plotX(end),plotY(end),distEnd,'facecolor','none','edgeColor',cmap(idxCMap(end,1),:));
                end
                
                scatter(xMax,yMax,50,'k');
                
                if ip.Results.plotScaleBar 
                      pixSizeMic = movieData.pixelSize_/1000; 
                      width  = ip.Results.maxCMapValue./pixSizeMic; 
                      plotScaleBar(width,2,'Location','SouthWest','Color',[0 0 0]);  
                end 
                
                %     measIndicesMask = zeros(ny,nx);
                %     measIndicesMask(measIndices) = 1;
                arrayfun(@(i) cellfun(@(x) plot(x(:,2),x(:,1),color(i)),roiYX{i}),1:2);
                hold on
                %spy(measIndicesMask,'k');
                %     text(5,5,' Within 20 um from tip of leading protrusion','color','m');
                %     text(5,20,'> 20 um from tip ','color','b');
                saveas(gcf, [ cDir filesep num2str(iFrame,'%03d') '.png']);
                saveas(gcf,[cDir filesep num2str(iFrame,'%03d') '.fig']); 
                close gcf
            end
        end
        % add the timestamp and save the parameters 
        params.timeStamp = clock; 
        p(iFrame) = params; 
        save([subRegDir filesep 'params.mat'],'p'); 
    end
    
    
    
    
else
    display(['GC SubRois Found for ' movieData.outputDirectory_ ': Skipping' ]) ;
end
