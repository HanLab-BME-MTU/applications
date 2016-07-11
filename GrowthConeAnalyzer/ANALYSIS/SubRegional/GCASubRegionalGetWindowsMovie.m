function [ idxWindFinal ] = GCASubRegionalGetWindowsMovie( movieData,varargin)
% GCASubRegionalGetWindowsMovie: finds windows in a subroi that last
% throughout the time interval chosen 
%
%%
%
% OUTPUT: 
% idxWindFinal 
% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;
%ip.addParameter('TSOverlays',true);
defaultInDir = [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep 'GCASubRegions' filesep  'GC' filesep 'masks'];
ip.addParameter('subRoiMaskDirIn',defaultInDir); % where the subRoimasks are stored.
ip.addParameter('windMethod','ConstantNumber'); % what windowing method to use - currently assumes you will use one of the
% windowing methods associated with the MD.
defaultOutDir = [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep 'GCASubRegions' filesep 'GC' filesep 'subRegion_windows'];
ip.addParameter('OutputDirectory_', defaultOutDir);
ip.addParameter('EndFrame',movieData.nFrames_-1);
ip.addParameter('StartFrame',1);
ip.addParameter('FillSpatialHoles',true); 
% ip.addParameter('OverlayColorByVel',true ); % this option will find the first and last window
% ip.addParameter('OverlayWindTracker',true); 
% ip.addParameter('OverlayOutlierFilter',true);  % outlier  
% ip.addParameter('OverlaySignalDetectType',{'P','R'}); % will only plot protrusion or retraction or quiescent events 

ip.addParameter('ReInit',61); 
% number to be defined in the subregion throughout the whole movie to
% define which windows are included. This allows for temporal gaps in the
% windows (ie where a window disappears and reappears) to be included in
% the subRoi.
ip.parse(varargin{:});
%% Initiate
% find the correct window output file.
idxWindProc =  cellfun (@(x) strcmpi(x.name_,'Windowing'),movieData.processes_);
windProcesses =  movieData.processes_(idxWindProc);
toLoad = cellfun(@(x) strcmpi(x.funParams_.MethodName, (ip.Results.windMethod)) & x.funParams_.ReInit==ip.Results.ReInit ,windProcesses);
windDir = windProcesses{toLoad}.outFilePaths_;
%windDir = movieData.processes_{toLoad}.outFilePaths_;
listOfWindFiles = searchFiles('.mat',[],windDir,0,'all',1);
if ~isdir(ip.Results.OutputDirectory_)
    mkdir(ip.Results.OutputDirectory_);
end

% if (ip.Results.OverlayColorByVel == true || ip.Results.OverlayWindTracker == true)
%     forName{1,1} = 'ColorByVelocity'; 
%     forName{2,1} = 'WindowTracker';
%     forName = forName([ip.Results.OverlayColorByVel;ip.Results.OverlayWindTracker],:); 
%     
%     if ~isempty(ip.Results.OverlayOutlierFilter) || ~isempty(ip.Results.OverlaySignalDetectType);  
%     idxSampProc =  cellfun (@(x) strcmpi(x.name_,'Protrusion Sampling'),movieData.processes_);
%     sampProcesses =  movieData.processes_(idxSampProc);
%     
%     load(sampProcesses{toLoad}.outFilePaths_{1}) ;
%     end    
% end
   
%%% Make OutputDirectories %%%  

% for now include the windMethod in outputdirectory name for organization 
outputDir_CMethodType = [ip.Results.OutputDirectory_ filesep ip.Results.windMethod]; 
if ~isdir(outputDir_CMethodType) 
    mkdir(outputDir_CMethodType) 
end 

%% % Get SubRegional Windows 

if ~isempty(ip.Results.subRoiMaskDirIn)
    
for iFrame = ip.Results.StartFrame:ip.Results.EndFrame
    % load the windows for the current frame
    load(listOfWindFiles{iFrame});
    
    roiMask  = logical(imread(([ip.Results.subRoiMaskDirIn filesep 'mask' num2str(iFrame,'%03d') '.tif'])));
    
    % load the subRoiMask
    [idxWindTestLog] = GCASubRegionalGetWindows(windows,roiMask);
    idxWindTest = find(idxWindTestLog);
    
    % get only winds that likewise existed in the preceding frame
    if iFrame>ip.Results.StartFrame
        % rewrite
        idxWindFinal=  intersect(idxWindTest , idxWindFinal);
    else
        idxWindFinal = idxWindTest;
    end
    
end

if ip.Results.FillSpatialHoles
    idxWindFinal = idxWindFinal(1):idxWindFinal(end);
end

% find start and end points

save([outputDir_CMethodType filesep 'idxWindFinal_' num2str(ip.Results.StartFrame) '_' num2str(ip.Results.EndFrame) '.mat' ],'idxWindFinal');

% if plot make a folder to test the windows.
% find the window numbers common to all
% chec
else 
    idxWindFinal = []; 
end  
  %% end subRegionalSpecific   
  
  
  
  
%% Make Overlays if Applicable 
% 
% if ~isempty(forName)
%      filenames=movieData.getImageFileNames;
%      
% for iPlot = 1:numel(forName)    
%     % mkClrDir([ip.Results.OutputDirectory_ filesep 'TSOverlays']);
%     % mkClrDir([ip.Results.OutputDirectory_ filesep 'WindowsAll']);
%     if ~isdir([outputDir_CMethodType  filesep 'TSOverlays_'  forName{iPlot}]);
%         mkdir([outputDir_CMethodType filesep 'TSOverlays_'  forName{iPlot}]);
%     end
% 
% 
% %%  If label by signal detect load analysis results %%% 
%             
%             if (ip.Results.OverlayOutlierFilter || ~isempty(ip.Results.OverlaySignalDetectType))
%                 
%                 sigDetectDir = upDirectory(sampProcesses{toLoad}.outFilePaths_{1},1);
%                 for iA = 1:2
%                     load([sigDetectDir filesep 'EdgeVelocityQuantification_CutMovie_' num2str(iA) filesep 'EdgeMotion.mat']);
%                     analysis{iA} = analysisResults;
%                     clear analysisResults;
%                 end
%             end           
%            %% Start Overlay Wrapper
%     for iFrame = ip.Results.StartFrame:ip.Results.EndFrame
%        
%         
%         
%         
%         
%         
%         img = double(imread([movieData.getChannelPaths{1} filesep filenames{1}{iFrame}]));
%         
%         if ~isempty(ip.Results.subRoiMaskDirIn)
%             
%             roiMask  = logical(imread(([ip.Results.subRoiMaskDirIn filesep 'mask' num2str(iFrame,'%03d') '.tif'])));
%             roiYX = bwboundaries(roiMask);
%         end
%         
%         load(listOfWindFiles{iFrame});
%         [ny,nx] = size(img);
%         
%          nWinds = numel(windows); 
%           
%         if strcmpi(forName{iPlot},'ColorByVelocity');
%             cmap = brewermap(128,'RdBu');
%         else
%             cmap = lines(nWinds);
%         end
%          
%          if isempty(idxWindFinal)
%              idxWindFinal = 1:nWinds;   
%          end 
%  %% Refine by Signal Detect Process if necessary 
% if iFrame<ip.Results.ReInit;
%     analysisResults = analysis{1};
%     iFrameFind = iFrame; 
% else
%     analysisResults = analysis{2};
%     iFrameFind = iFrame - ip.Results.ReInit+1; % the values for each of the analysisResults
% end
%  
% if  ~isempty(ip.Results.OverlaySignalDetectType);  
%     
%     idxProtC = []; 
%     idxRetC = [] ; 
%     idxQuies = []; 
%     
%     if sum(strcmpi(ip.Results.OverlaySignalDetectType,'prot')) >0 
%           idxProtC =  arrayfun(@(x) ~isempty(find(vertcat(analysisResults.protrusionAnalysis.windows(x).blockOut{:})==iFrameFind)),1:nWinds);
%           idxProtC= find(idxProtC); 
%           
%     end
%     
%     if sum(strcmpi(ip.Results.OverlaySignalDetectType,'retract'))>0 
%          idxRetC = arrayfun(@(x) ~isempty(find(vertcat(analysisResults.retractionAnalysis.windows(x).blockOut{:}) == iFrameFind)),1:nWinds);
%          idxRetC = find(idxRetC); 
%     end 
%     
%      
%     if sum(strcmpi(ip.Results.OverlaySignalDetectType,'quies'))>0 
%       
%           idxProtC =  arrayfun(@(x) ~isempty(find(vertcat(analysisResults.protrusionAnalysis.windows(x).blockOut{:})==iFrameFind)),1:nWinds);
%          
%           idxRetC = arrayfun(@(x) ~isempty(find(vertcat(analysisResults.retractionAnalysis.windows(x).blockOut{:}) == iFrameFind)),1:nWinds);
%       
%           idxQuies = ~idxProtC & ~idxRetC; 
%        
%           idxQuies = find(idxQuies); 
%     end 
%     
%     
%     idxKeep = [idxProtC,idxRetC,idxQuies]; 
%     idxWindFinalC = intersect(idxKeep,idxWindFinal); 
% else 
%     idxWindFinalC = idxWindFinal;
% 
% end % isempty(ip.Results.OverlaySignalDetectType
%                 
%  %%
%  
%         if  strcmpi(forName{iPlot},'ColorByVelocity')
%          
%             setFigure(nx,ny,'on');
%             
%             imshow(-img,[]) ;
%             hold on
%             if ~isempty(ip.Results.subRoiMaskDirIn)
%                 
%                 roiMask  = logical(imread(([ip.Results.subRoiMaskDirIn filesep 'mask' num2str(iFrame,'%03d') '.tif'])));
%                 roiYX = bwboundaries(roiMask);
%                 cellfun(@(x) plot(x(:,2),x(:,1),'color','k'),roiYX);
%             end
%            
%             % plot with no color to get the window numbers
%             gcaPlotWindows(windows,{'k','FaceAlpha',0},'bandMax',1,'showNum',10);
%             
% 
%     plotValues = protSamples.avgNormal(idxWindFinalC,iFrame);
%     plotValues = plotValues*movieData.pixelSize_/movieData.timeInterval_; 
%             
%     mapper=linspace(100,-100,128)'; % lower values are red in the brewermap
%     D=createDistanceMatrix(plotValues,mapper);
%     [sD,idxCMap]=sort(abs(D),2);
% 
%             for iColor = 1:length(cmap) 
%                 windToPlot = idxWindFinalC(idxCMap(:,1)==iColor);         
%                 if ~isempty(windToPlot)
%                     for iWind = 1:length(windToPlot) % I know it is silly but do each individually for now - just because it is the way I can
%                         % get the boxes overlaid -
%                         gcaPlotWindows(windows(windToPlot(iWind)),{cmap(iColor,:),'FaceAlpha',1},windToPlot(iWind),'bandMax',1,'colorWind',cmap(iColor,:));
%                         
%                     end
%                 end
%             end 
%            
%         % find all windows that were protrusion events that passed the
%         % time series analysis criteria in the current frame
%         
%        % nWinds = length(analysisResults.protrusionAnalysis.windows(:));
%         nWindsC = numel(windows);
%         idxBlackOut = arrayfun(@(x) ~isempty(find(find(isnan(analysisResults.data.procEdgeMotion(x,:)))==iFrame)),1:nWindsC);
%         idxBlackOut = find(idxBlackOut); 
%         %idxBlackOut == idxWindFinalC; 
%         idxBlackOut = intersect(idxBlackOut,idxWindFinalC); 
%           if ~isempty(idxBlackOut) 
%               gcaPlotWindows(windows(idxBlackOut),{'k','FaceAlpha',1},'bandMax',1); 
%           end 
%          clear idxBlackOut
%          
%          
%        
%          
%          
%         else  % color by number (default)
%            
%             setFigure(nx,ny,'off');
%             
%             imshow(-img,[]) ;
%             hold on
%             if ~isempty(ip.Results.subRoiMaskDirIn)                
%                 cellfun(@(x) plot(x(:,2),x(:,1),'color','k'),roiYX);
%             end
%             % plot with no color to get the window numbers
%             gcaPlotWindows(windows,{'k','FaceAlpha',0},'bandMax',1,'showNum',10);     
%            
%             for iWind = 1:length(idxWindFinalC)
%                 
%                 gcaPlotWindows(windows(idxWindFinalC(iWind)),{cmap(iWind,:),'FaceAlpha',1},idxWindFinalC(iWind),'bandMax',1,'colorWind',cmap(iWind,:));   
%             end      
%         end
%         
%         saveas(gcf,[outputDir_CMethodType  filesep 'TSOverlays_' forName{iPlot} filesep num2str(iFrame,'%03d') '.png']);
%         close gcf
%         %        setFigure(ny,nx,'on');
%         %        imshow(-img,[]);
%         %        hold on
%         %        % plot all the windows as well as a check so can see the number
%         %        gcaPlotWindows(windows,'bandMax',1,'showNum',5);
%         %        saveas(gcf, [ip.Results.OutputDirectory_ filesep 'WindowsAll' filesep num2str(iFrame,'%03d') '.png']);
%         
%         
%         %        close gcf
%     end
% end
% 
% 
% 
% 
% end


