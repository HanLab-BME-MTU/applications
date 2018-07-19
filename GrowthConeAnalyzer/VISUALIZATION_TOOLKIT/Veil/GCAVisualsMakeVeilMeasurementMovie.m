function [ output_args ] = GCAVisualsMakeVeilMeasurementMovie(movieData,varargin)
% GCAVisualsMakeVeilMeasurementMovie (FINISH ME!!- most of this functionality has been implemented in GCASubRegionalGetWindowsMovie
% just needs to be split  
% for now keep this completely separete from the filopodia (eventually
% would like to merge
% Features to include: 
%  SubRegions, Signal Processing overlays, color by velocity and, veil tracker capabilities 
% 
% 
% OUTPUT: 

if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('OutputDirectory_', []);
ip.addParameter('OverlayType','ColorByVel'); 
ip.addParameter('ScaleBar',false); 
ip.addParameter('ShowWindowNum',true);

% SubRoiOption
ip.addParameter('subRoiMaskDirIn',[]); % where the subRoimasks are stored.
ip.addParameter('OutlierFilter',false); 
ip.addParameter('StartFrame',1); 
ip.addParameter('EndFrame',119);

ip.addParameter('cLims',[-100,100]); % currently in nm/sec, set's the scale limits for the colormap

ip.addParameter('SignalType',[]); % {'P'},{'R'},{'Q'} if empty will plot all signals 

% Protrusion Map Settings 
ip.addParameter('MakeProtMap',true); 
ip.addParameter('SmoothActivityMap',true); 

ip.addParameter('windMethod','ConstantNumber'); 
ip.addParameter('ReInit',61); 

ip.addParameter('firstFrameLimits',false);

ip.addParameter('overlayLongPath',false); 
ip.addParameter('neuriteElongDir',[]); 

ip.addParameter('screen2png',false); 
ip.addParameter('highResImage',false); 

ip.addParameter('collectDir',[]); 

ip.parse(varargin{:});

%% Initiate
% find the correct window output file.
idxWindProc =  cellfun (@(x) strcmpi(x.name_,'Windowing'),movieData.processes_);
windProcesses =  movieData.processes_(idxWindProc);

if ~isempty(ip.Results.ReInit)
    toLoad = cellfun(@(x) strcmpi(x.funParams_.MethodName, (ip.Results.windMethod)) & x.funParams_.ReInit==ip.Results.ReInit ,windProcesses);
else
    toLoad = cellfun(@(x) strcmpi(x.funParams_.MethodName, (ip.Results.windMethod))  ,windProcesses);
end
windDir = windProcesses{toLoad}.outFilePaths_;

%windDir = movieData.processes_{toLoad}.outFilePaths_;
listOfWindFiles = searchFiles('.mat',[],windDir,0,'all',1);

if (isempty(ip.Results.subRoiMaskDirIn) && isempty(ip.Results.OutputDirectory_))
    outDir =  [movieData.outputDirectory_ filesep 'VisualizationOverlays' filesep 'WholeNeurite' filesep 'VeilWindows'];
elseif   ~isempty(ip.Results.subRoiMaskDirIn) && isempty(ip.Results.OutputDirectory_)
    outDir = [ip.Results.subRoiMaskDirIn filesep 'subRegion_windows']; 
else 
    outDir = ip.Results.OutputDirectory_; 
end 

 outDir = [outDir filesep ip.Results.windMethod filesep ip.Results.OverlayType];

% get image filenames 
filenames=movieData.getImageFileNames;



if strcmpi(ip.Results.OverlayType,'ColorByVel');  
    % load the signal detection if applicable
    
    idxSampProc =  cellfun (@(x) strcmpi(x.name_,'Protrusion Sampling'),movieData.processes_);
    sampProcesses =  movieData.processes_(idxSampProc);
    
    load(sampProcesses{toLoad}.outFilePaths_{1}) ;
end
    
if (ip.Results.OutlierFilter || ~isempty(ip.Results.SignalType))   
    sigDetectDir = upDirectory(sampProcesses{toLoad}.outFilePaths_{1},1);
    for iA = 1:2
        load([sigDetectDir filesep 'EdgeVelocityQuantification_CutMovie_' num2str(iA) filesep 'EdgeMotion.mat']);
        analysis{iA} = analysisResults;
        clear analysisResults;
    end
    
end  % if signal processing information required

if ~isempty(ip.Results.subRoiMaskDirIn)
    load([ip.Results.subRoiMaskDirIn  filesep 'subRegion_windows' filesep ip.Results.windMethod filesep 'idxWindFinal_' ... 
        num2str(ip.Results.StartFrame) '_' num2str(ip.Results.EndFrame) '.mat']);  
end 

%% Make output directory 
% 

if ~isempty(ip.Results.SignalType)
    detectType = horzcat(ip.Results.SignalType{:});
    outDir = [outDir '_SigDetect_' detectType];
end

if ip.Results.OutlierFilter
    outDir = [outDir '_OutlierDetect'];
end


if ip.Results.overlayLongPath

    if isempty(ip.Results.neuriteElongDir);
        neuriteElongDir  = [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep 'StepsToReconstruct' filesep ...
            'IV_veilStem_length' filesep 'Channel_1'];
    else 
        neuriteElongDir =  ip.Results.neuriteElongDir; 
    end
  
    load([neuriteElongDir filesep 'veilStem.mat']); 
end 

if ~isdir(outDir)
    mkdir(outDir)
end

if ~isempty(ip.Results.StartFrame)
%% Start Overlay loop 
  for iFrame = ip.Results.StartFrame:ip.Results.EndFrame
            
      img = double(imread([movieData.getChannelPaths{1} filesep filenames{1}{iFrame}]));
      
      
        if  iFrame ==1 
            if ip.Results.firstFrameLimits
            lims = [min(-img(:)) max(-img(:))]; 
            else 
                lims = []; 
            end 
        
        end 
      if ~isempty(ip.Results.subRoiMaskDirIn)         
          % currently the input is such that 
          roiMask  = logical(imread(([ip.Results.subRoiMaskDirIn filesep 'masks' filesep  'mask' num2str(iFrame,'%03d') '.tif'])));
          roiYX = bwboundaries(roiMask);
      end
      
     
      load(listOfWindFiles{iFrame});
     
      
      [ny,nx] = size(img);
      
      nWinds = numel(windows);
      
      setFigure(nx,ny,'off');
      
      imshow(-img,lims) ;
      hold on
      
      if ~isempty(ip.Results.subRoiMaskDirIn)
          cellfun(@(x) plot(x(:,2),x(:,1),'color','k'),roiYX);
         idxWindFinalC = idxWindFinal; 
         hold on
      else 
          idxWindFinalC = 1:length(windows);
          
      end
      
      if ip.Results.overlayLongPath
          % plot the veil/stem 
          l = veilStem(iFrame).neuriteLongPathIndices;
         
          lmask = zeros(size(img)); 
          lmask(veilStem(iFrame).neuriteLongPathIndices)= 1;
          spy(lmask,'k');
      end 
      
      
      
      switch ip.Results.OverlayType
      %% Start Window Tracker  
          case 'WindTrack'  % window tracker
                 
           
                    
              % color by number
              cmap = lines(nWinds);
              
              if ip.Results.ShowWindowNum
                  % plot with no color to get the window numbers
                  gcaPlotWindows(windows,{'k','FaceAlpha',0},'bandMax',2,'bandMax',2,'showNum',10);
              end    
                  for iWind = 1:length(idxWindFinalC)
                      gcaPlotWindows(windows(idxWindFinalC(iWind)),{cmap(iWind,:),'FaceAlpha',1},idxWindFinalC(iWind),'bandMin',2,'bandMax',2,'colorWind',cmap(iWind,:));
                  end
             
          %% start ColorByVelocity 
          case 'ColorByVel'
                        
              % quick fix for my workflow - I store the Re-Int timeseries
              % in two separate analysis (1:60) typically and (61-120)
              if (ip.Results.OutlierFilter || ~isempty(ip.Results.SignalType))
                  if iFrame<ip.Results.ReInit;
                      analysisResults = analysis{1};
                      iFrameFind = iFrame;
                  else
                      analysisResults = analysis{2};
                      iFrameFind = iFrame - ip.Results.ReInit+1; % the values for each of the analysisResults
                  end
              end
              
              %% Filter by Signal Type Assignment if required
              if  ~isempty(ip.Results.SignalType);
                 
                  idxProtC = [];
                  idxRetC = [] ;
                  idxQuies = [];
                  
                  if sum(strcmpi(ip.Results.SignalType,'P')) >0
                      idxProtC =  arrayfun(@(x) ~isempty(find(vertcat(analysisResults.protrusionAnalysis.windows(x).blockOut{:})==iFrameFind)),1:nWinds);
                      idxProtC= find(idxProtC);
                      
                  end
                  
                  if sum(strcmpi(ip.Results.SignalType,'R'))>0
                      idxRetC = arrayfun(@(x) ~isempty(find(vertcat(analysisResults.retractionAnalysis.windows(x).blockOut{:}) == iFrameFind)),1:nWinds);
                      idxRetC = find(idxRetC);
                  end
                  
                  
                  if sum(strcmpi(ip.Results.SignalType,'Q'))>0 
                      
                      idxProtC =  arrayfun(@(x) ~isempty(find(vertcat(analysisResults.protrusionAnalysis.windows(x).blockOut{:})==iFrameFind)),1:nWinds);
                      
                      idxRetC = arrayfun(@(x) ~isempty(find(vertcat(analysisResults.retractionAnalysis.windows(x).blockOut{:}) == iFrameFind)),1:nWinds);
                      
                      idxQuies = ~idxProtC & ~idxRetC;
                      
                      idxQuies = find(idxQuies);
                  end
                  
                  
                  idxKeep = [idxProtC,idxRetC,idxQuies];
                  idxWindFinalC = intersect(idxKeep,idxWindFinalC);
              else
                  idxWindFinalC = idxWindFinalC;
                  
              end % isempty(ip.Results.OverlaySignalDetectType
              
                       
              %% Create Mapper %%
              % initiate the colormap (% make default)
              cmap = brewermap(128,'RdBu');
              
              % get the average velocity values for all the windows in the current frame
              % and assign a color based on the the mapper.
              plotValues = protSamples.avgNormal(idxWindFinalC,iFrame);
              plotValues = plotValues*movieData.pixelSize_/movieData.timeInterval_;
              
              mapper=linspace(ip.Results.cLims(2),ip.Results.cLims(1),128)'; % lower values are red in the brewermap
              D=createDistanceMatrix(plotValues,mapper);
              [sD,idxCMap]=sort(abs(D),2);
              
              
              %% Plot the windows
              
              
                %NOTE 20160607
                % for movies you want this to be plotted first and the
                % patches 2nd - this will look smooth
                % 
%               if ip.Results.ShowWindowNum
%                   % plot with no color to get the window numbers
%                   gcaPlotWindows(windows,{'k','FaceAlpha',0},'bandMax',1,'showNum',10);
%               end
%               
%               
          
   


              % Plot Events Selected 
              for iColor = 1:length(cmap)
                  windToPlot = idxWindFinalC(idxCMap(:,1)==iColor);
                  if ~isempty(windToPlot)
                      for iWind = 1:length(windToPlot) % I know it is silly but do each individually for now - just because it is the way I can
                          % get the boxes overlaid -
                          gcaPlotWindows(windows(windToPlot(iWind)),{cmap(iColor,:),'FaceAlpha',1},windToPlot(iWind),'bandMax',1,'colorWind',cmap(iColor,:));
                          
                      end
                  end
              end
            
           % Note 20160607: note for vectorization (ie .eps files and figures) want this to be
              % plotted second - somehow having the lines plotted second signals it to vectorize and not
              % smooth
              if ip.Results.ShowWindowNum
                  % plot with no color to get the window numbers
                  gcaPlotWindows(windows,{'k','FaceAlpha',0},'bandMax',1,'showNum',10);
              else 
                  %gcaPlotWindows(windows,{'k','FaceAlpha',0},'bandMax',1,'lineWidth',.01); 
              end
              

              
              
              if ip.Results.ScaleBar == true
                  pixSizeMic = movieData.pixelSize_./1000;
                  plotScaleBar(10/pixSizeMic,1,'Color',[0,0,0],'Location','SouthEast');
              end
              
              
              %% Black Out Outlier Windows if Neccessary 
              if  ip.Results.OutlierFilter;     % filter out the outliers based on the criteria
                
                  c = brewermap(9,'Set1'); 
                  outlierC = c(end,:); 
                  
                  % nWinds = length(analysisResults.protrusionAnalysis.windows(:));
                  nWindsC = numel(windows);
                  idxBlackOut = arrayfun(@(x) ~isempty(find(find(isnan(analysisResults.data.procEdgeMotion(x,:)))==iFrameFind)),1:nWindsC);
                  idxBlackOut = find(idxBlackOut);
                  %idxBlackOut == idxWindFinalC;
                  %idxBlackOut = intersect(idxBlackOut,idxWindFinalC);
                  if ~isempty(idxBlackOut)
                      
                     % if ~isempty(windToPlot)
                      for iWind = 1:length(idxBlackOut) % I know it is silly but do each individually for now - just because it is the way I can
                          % get the boxes overlaid -
                          gcaPlotWindows(windows(idxBlackOut(iWind)),{outlierC,'FaceAlpha',1},idxBlackOut(iWind),'bandMax',1,'colorWind',outlierC);
                          
                      end
                  end
                      
                      
%                         clear idxBlackOut
                      
                    %  gcaPlotWindows(windows(idxBlackOut),{'y','FaceAlpha',1},0,'bandMax',1);
                
                
              end % if output
                          
      end  % switch
      
       %% save the results 
       if ip.Results.screen2png
          helperScreen2png([outDir filesep num2str(iFrame,'%03d') '.png']); 
       else 
       saveas(gcf,[outDir filesep num2str(iFrame,'%03d') '.png']);
       end
       
       if ip.Results.highResImage
           hiResFolder = [outDir filesep '600dpi' ];
           if ~isdir(hiResFolder)
               mkdir(hiResFolder);
           end
           print(gcf,'-dpng',[hiResFolder filesep num2str(iFrame,'%03d') '.png'], '-r600');
           
           
       end
       
       saveas(gcf,[outDir filesep num2str(iFrame,'%03d') '.fig']);
       saveas(gcf,[outDir filesep num2str(iFrame,'%03d'),'.eps']);
       
    
       
       
       close gcf
           
  end % iframe  
  
  if ~isempty(ip.Results.collectDir)
      copyfile([outDir,collectDir]);
  end
end  
if ip.Results.MakeProtMap
     cmap = brewermap(128,'RdBu');
    % make the directory for this.
    mapDir = [outDir filesep 'protrusionMap']; 
    if ~isdir(mapDir)
        mkdir(mapDir); 
    end 
    
    mapValues =  protSamples.avgNormal;
    mapValues = mapValues*movieData.pixelSize_/movieData.timeInterval_;
    
    
    % Always make one protrusion map with the original map values
    cmapC = flip(cmap,1);
    add = [0,0,0]; % plot NaN in black 
    cmapFinal = [add; cmapC]; 
    
    if ip.Results.SmoothActivityMap 
     processed = smoothActivityMap(mapValues,'upSample',1,'SmoothParam',0.75); 
    end 
      
     GCAVisualsProtrusionMap(mapValues,'colorMap',cmapFinal,'CLim',[-100,100]);  
    
     saveas(gcf,[mapDir filesep 'protrusionMapUnSmoothed.png']);
     saveas(gcf,[mapDir filesep 'protrusionMapUnSmoothed.fig']);
     
     
     close gcf 
     GCAVisualsProtrusionMap(processed,'colorMap',cmapFinal,'CLim',[-100,100]); 
     saveas(gcf,[mapDir filesep 'protrusionMapSmoothed.tif']);
     saveas(gcf,[mapDir filesep 'protrusionMapSmoothed.fig']);
     clear cmapFinal
     close gcf
end 

    
if (ip.Results.MakeProtMap && ip.Results.OutlierFilter)
    cmapC = flip(cmap,1);
    % pad 
%     forPad = cmapC(1,:); 
%     forPad2 = cmapC(end,:); 
    
  
%     add = [0,0,0];
%     add2 = outlierC;
%     cmapFinal = [add ; cmapC ; add2];
    
    raw = analysisResults.data.rawTimeSeries;
    raw(isnan(raw)) = -inf; % set NaN to lowest value ;
    %         raw(isnan(raw)) = -inf; % set NaN to lowest value ;
    %idxExclude = vertcat(analysis{1}.data.excludedWin{:});
    % REMEMBER : Marco does the outlier detection on the ENTIRE time series
    % not the partition- so it doesn't ma
    outlierOut = analysis{1}.data.procEdgeMotion;
    %        processed1 = analysis{1}.data.procTimeSeries; % I am not sure if marco processed the whole time
    outlierOut(raw==-inf) = -inf;
    outlierOut(isnan(outlierOut)) = inf;
    
    data = raw.*analysisResults.data.scaling; 
    
    if ip.Results.SmoothActivityMap
        data = smoothActivityMap(data,'upSample',1,'SmoothParam',0.75);
    end
    
    
    %        processed2 = analysis{2}.data.procTimeSeries;
    %        processed2(raw==-inf) = -inf;
    %        processed2(isnan(processed2)) = inf;
    %        isequal(processed,processed2);
    %
    %
    %
    %       processed = [processed1 processed2];
    %
    %GCAVisualsProtrusionMap(raw.*analysisResults.data.scaling,'colorMap',cmapC,'CLim',[-100,100]);
    GCAVisualsProtrusionMap(data,'colorMap',cmapC,'CLim',ip.Results.cLims,'ShowColorBar',false,'FontSize',50);
    
    
    hold on 
    notANum = find(outlierOut==-Inf);
    if ~isempty(notANum)
    [y,x] = ind2sub(size(outlierOut),notANum); 
    scatter(x,y,'k','filled'); 
    end 
    
    outliers = find(outlierOut == Inf);
    if ~isempty(outliers)
        [y,x] = ind2sub(size(outlierOut),outliers); 
        scatter(x,y,'k','x'); 
    end 
    
    
    saveas(gcf,[mapDir filesep 'protrusionMapRemoveSignalDetectOutlier.tif']);
    saveas(gcf,[mapDir filesep 'protrusionMapRemoveSignalDetectOutlier.fig']);
    saveas(gcf,[mapDir filesep 'protrusionMapRemoveSignalDetectOutlier.eps'],'psc2'); 
    % series for outliers or not.
    
    
    
    %processed2 = analysis{2}.data.procTimeSeries;
    
    %cmap =  cmap(end:
    
    close gcf 
end
 
end

