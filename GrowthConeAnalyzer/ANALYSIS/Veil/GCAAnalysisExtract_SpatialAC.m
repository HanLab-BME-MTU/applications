function [ output_args ] = GCAAnalysisExtract_SpatialAC(movieData,varargin)
% GCAAnalysisExtract_SpatialAC 

%% CheckInput
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end



%%Input check
ip = inputParser;
ip.CaseSensitive = false;

ip.addParameter('OutputDirectory',[]); 

ip.addParameter('subRegionFlag',false);
ip.addParameter('subRoiMaskDirIn',[]); 
ip.addParameter('subRoiTSOverlays',false); 
ip.parse(varargin{:});

%%
iProtProc = movieData.getProcessIndex('ProtrusionProcess',1);

if ~movieData.processes_{iProtProc}.checkChannelOutput
    error('Movie must have valid protrusion vectors calculated first!')
end


%% if subRegional load the masks 
 if ip.Results.subRegionFlag
     if isempty(ip.Results.subRoiMaskDirIn) ; 
         veilDir = [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep ... 
             'GCASubRegions' filesep 'GC' filesep 'masks']; 
     else 
         veilDir = ip.Results.subRoiMaskDirIn;  
     end 
 end 


tmp = movieData.processes_{iProtProc}.loadChannelOutput;

protrusion = tmp.protrusion;
smoothedEdge = tmp.smoothedEdge;
normals = tmp.normals;

nPtsTot = 0;
nFrames = movieData.nFrames_;

%outDir = [movieData.outputDirectory_  filesep ...
%   'PARAMETER_EXTRACTION' filesep 'Spatial_AutoCorrelation_Test'];

if isempty(ip.Results.OutputDirectory)
    if ip.Results.subRegionFlag
        outDir = [movieData.outputDirectory_ filesep  'GCAMeasurementExtraction' filesep 'GCASubRegions' ...
            filesep 'GC' filesep 'Descriptor' filesep 'Veil' filesep 'SpatialACF'];
    else
        outDir = [movieData.outputDirectory_ filesep ...
            'MEASUREMENT_EXTRACTION' filesep 'Descriptor' filesep 'Veil' filesep ...
            'SpatialACF'];
    end
else
    outDir = ip.Results.OutputDirectory;
end



if ~isempty(outDir)
    mkdir(outDir)
end

for iFrame = 1:(nFrames-1)
    
      if ip.Results.subRegionFlag
           roiMask  = logical(imread(([veilDir filesep 'mask' num2str(iFrame,'%03d') '.tif'])));
           roiMask(1,1)= 0; 
           smoothedEdgeC = smoothedEdge{iFrame};
           [ny,nx] = size(roiMask); 
           idxAll = sub2ind([ny,nx],round(smoothedEdge{iFrame}(:,2)),round(smoothedEdge{iFrame}(:,1)));
           IN = roiMask(idxAll); 
           idxKeep = find(IN,1,'first'):find(IN,1,'last') ; 
           
           if ip.Results.subRoiTSOverlays
               setFigure(nx,ny,'off');
               
               smoothedEdgeTrunc = smoothedEdgeC(idxKeep,:);
               
               
               
               test = bwboundaries(roiMask);
               
               scatter(test{1}(:,2),test{1}(:,1),10,'r','filled');
               hold on
               scatter(round(smoothedEdge{iFrame}(:,1)),round(smoothedEdge{iFrame}(:,2)) ,50,'k');  
               
          
               scatter(round(smoothedEdgeTrunc(:,1)),round(smoothedEdgeTrunc(:,2)),50,'b');
               
               if ~isdir([outDir filesep 'TSOverlays']); 
                   mkdir([outDir filesep 'TSOverlays']); 
               end 
               
               saveas(gcf,[outDir filesep 'TSOverlays' filesep num2str(iFrame,'%03d') '.png']); 
               close gcf
               %scatter(x(:,1),x(:,2),'g');
           end
                  
           
      end  
    
    %Normalize the normals - sam's function does not return them with unit
    %length
    normals{iFrame} = normals{iFrame} ./ ...
        repmat(sqrt(dot(normals{iFrame}',normals{iFrame}')'),1,2);
    
    
    %Get the normal component of the protrusion at each point
    ncCurr = dot(normals{iFrame}',protrusion{iFrame}')';
    %Get the magnitude of the protrusion at each point
    mCurr = sqrt(dot(protrusion{iFrame}',protrusion{iFrame}'))';
    
    if ip.Results.subRegionFlag
        ncCurr = ncCurr(idxKeep);
        mCurr = mCurr(idxKeep); 
        
    end 
    
    %     %use Khulouds splitting trick - divide more than once??? Why only two??
    %     %Middle justified? Or something fancier?
    %     nPtsCurr = numel(mCurr);
    %     nHalf = round(nPtsCurr/2);
    %     normalComponent(iFrame).observations = ncCurr(1:nHalf);
    %     normalComponent(iFrame+nFrames-1).observations = ncCurr(nHalf+1:end);
    %     magnitude(iFrame).observations = mCurr(1:nHalf);
    %     magnitude(iFrame+nFrames-1).observations = mCurr(nHalf+1:end);
    %     nPtsTot = nPtsTot + nPtsCurr;
    
    nPtsCurr = numel(mCurr);
    nPtsAll(iFrame) = nPtsCurr;
    normalComponents{iFrame} =  ncCurr;
    magnitude{iFrame} = mCurr;
    nPtsTot = nPtsTot + nPtsCurr;
    
end

% calculate the maximum lag based on the number of windows
maxLag1 = round(min(nPtsAll)/4);
% get the minimum number of prot calcs for a given window and divide by 4
maxLag = round(min(cellfun(@(x) length(x), normalComponents))/4);

% make an rxc double array where r is the number of observations and c is
% the time frame
NCpad = reformatDataCell(normalComponents); % r  = spatial metric, c = frame

NCpadInvert = NCpad; % r = frame (to average) and c = spatial metric (to correlate)

% get average correlation for spatial metric
[muCC,muCI,lags,ACF] =  getAverageCorrelation(NCpad','maxLag',maxLag);

muCCPos = muCC(lags>0);
% for now just set tMax at the value where muCC goes to 0.2.
maxSigLag = find(muCCPos<0.2,1,'first');


save([outDir filesep 'AvgAC'],'muCC','muCI','lags','ACF','maxSigLag');

% reformat into the parameter version (cell for each frame) - for now as
% just replicate the values as this

measC = arrayfun(@(x) maxSigLag,1:nFrames,'uniformoutput',0)';
save([outDir filesep 'meas_maxACFLagSpatial'],'measC');


%% Plot
makePlot = 1;

if makePlot == 1;
    setAxis('off'); 
    stem(lags,muCC,'fill');
    hold on
    arrayfun(@(i) plot(lags,muCI(i,:),'color','r','lineStyle','--'),1:2);
    xlabel('Lag (Pixels)','FontSize',16,'FontName','Arial');
    ylabel('AutoCorrelation','FontSize',16,'FontName','Arial');
    
    saveas(gcf,[outDir filesep 'SpatialAutoCorrelation.fig']);
    saveas(gcf,[outDir filesep 'SpatialAutoCorrelation.png']);
    %plotTransparent(lags,muCI);
    %h = plot(lags,muCC,'r','linewidth',2);
    close gcf
    
end

