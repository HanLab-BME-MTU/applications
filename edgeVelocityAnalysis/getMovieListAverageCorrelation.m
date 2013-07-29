function [cellData,dataSet] = getMovieListAverageCorrelation(movieObj,varargin)
% This function estimates the global and local cross-correlation between the sampled signal and edge velocity
%It also returns the autocorrelation for the edge and sampled signal
%
%USAGE:
%       [cellData,dataSet] = getMovieListAverageCorrelation(movieObj,varargin)
%
%Input:
%       movieObj      - movie list or movie data object
%       includeWin    - cell array with indexes of the windows to be included for the analysis. Each element of the array corresponds to a cell from the movie list
%       outLevel      - scalar for outlier detection(see detectOutlier)
%       trendType     - scalar - algorihtm used to detrend the data(See getTimeSeriesTrend)
%       minLength     - scalar - mininum time serie length
%       maxLag        - scalar - maximum cross/auto-correlation lag
%       layer         - scalar - layer for the sampled signal
%       signalChannel - scalar
%
%Output:
%       cellData - structure array with all parameters for each cell
%                                    cellData(iCell).meanValue
%                                                   .CI
%                                                   .total
%                                                         .edgeAutoCorr
%                                                         .signalAutoCorr
%                                                         .crossCorr
%                                                         .lag
%
%       dataSet - structure array with all parameters for the whole data set. Same fields as in cellData
%
%
%Marco Vilela, 2013

ip = inputParser;
ip.addRequired('movieObj',@(x) isa(x,'MovieList') || isa(x,'MovieData'));


if isa(movieObj,'MovieData')
    
    ML = movieData2movieList(movieObj);
    
else
    
    ML = movieObj;
    
end

nCell = numel(ML.movies_);

ip.addParamValue('includeWin', cell(1,nCell),@iscell);
ip.addParamValue('outLevel',  zeros(1,nCell),@isvector);
ip.addParamValue('trendType',-ones(1,nCell),@isvector);
ip.addParamValue('minLength', 30*ones(1,nCell),@isvector);
ip.addParamValue('gapSize',   zeros(1,nCell),@isvector);
ip.addParamValue('scale',false,@islogical);
ip.addParamValue('maxLag',0,@isscalar);
ip.addParamValue('layer',1,@isscalar);
ip.addParamValue('signalChannel',1,@isscalar);
ip.addParamValue('interval',{[]},@iscell);

ip.parse(movieObj,varargin{:});
scale            = ip.Results.scale;
includeWin       = ip.Results.includeWin;
outLevel         = ip.Results.outLevel;
minLen           = ip.Results.minLength;
trend            = ip.Results.trendType;
maxLag           = ip.Results.maxLag;
layer            = ip.Results.layer;
channel          = ip.Results.signalChannel;
interval         = ip.Results.interval;

cellData{1,nCell} = [];
edgeInputParam    = {'outLevel',outLevel,'minLength',minLen,'trendType',trend,'includeWin',includeWin,'gapSize',ones(1,nCell),'scale',scale,'outputPath','correlationEstimation'};
signalInputParam  = {'outLevel',outLevel,'minLength',minLen,'trendType',trend,'includeWin',includeWin,'gapSize',ones(1,nCell),'outputPath','correlationEstimation'};
edge              = edgeVelocityQuantification(ML,edgeInputParam{:});
signal            = sampledSignalQuantification(ML,channel,signalInputParam{:});
nInterval         = cellfun(@(x) 1:size(x.data.procEdgeMotion,2),edge,'Unif',0);



for iInt = 1:numel(interval)
    if ~isempty(interval{iInt})
        nInterval{iInt} = interval{iInt};
    end
end

totalEdgeACF     = [];
totalSignalACF   = [];
totalCCF         = [];
protrusionCCF    = [];
retractionCCF    = [];

for iCell = 1:nCell
    
    
    windows         = intersect(edge{iCell}.data.includedWin,signal{iCell}.data.includedWin{layer});
    motionMap       = edge{iCell}.data.procEdgeMotion(windows,nInterval{iCell});
    activity        = squeeze(signal{iCell}.data.procSignal(windows,nInterval{iCell},layer));
    
    cellData{iCell} = internalGetCorrelation(motionMap,activity,maxLag,1);
    
    totalEdgeACF    = cat(2,cellData{iCell}.total.edgeAutoCorr,totalEdgeACF);
    totalSignalACF  = cat(2,cellData{iCell}.total.signalAutoCorr,totalSignalACF);
    totalCCF        = cat(2,cellData{iCell}.total.crossCorr,totalCCF);
    
    protIdx         = arrayfun(@(x) cell2mat(x.blockOut(:)),edge{iCell}.protrusionAnalysis.windows(windows),'Unif',0);
    retrIdx         = arrayfun(@(x) cell2mat(x.blockOut(:)),edge{iCell}.retractionAnalysis.windows(windows),'Unif',0);
    
    protMask        = nan(size(motionMap));
    retrMask        = nan(size(motionMap));
    
    for iWin = 1:numel(windows)
        
        protMask(iWin,protIdx{iWin}') = 1;
        retrMask(iWin,retrIdx{iWin}') = 1;
        
    end
    
    protrusion    = motionMap.*protMask;
    retraction    = motionMap.*retrMask;
    sigProtrusion = activity.*protMask;
    sigRetraction = activity.*retrMask;
    
    cellData{iCell}.total.protrusionSignal = nanmean(sigProtrusion(:));
    cellData{iCell}.total.protrusionSignal = nanmean(sigRetraction(:));

    cellData{iCell}.protrusionCorrelation  = internalGetCorrelation(protrusion,activity,maxLag,0);
    cellData{iCell}.retractionCorrelation  = internalGetCorrelation(retraction,activity,maxLag,0);
    protrusionCCF                          = cat(2,cellData{iCell}.protrusionCorrelation.total.crossCorr,protrusionCCF);
    retractionCCF                          = cat(2,cellData{iCell}.retractionCorrelation.total.crossCorr,retractionCCF);
    
end

[dataSet.meanValue.crossCorr,dataSet.CI.crossCorr]                       = correlationBootstrap(totalCCF, repmat( 2/sqrt(size(motionMap,2)), 1, size(totalCCF,2) ) );
[dataSet.meanValue.edgeAutoCorr,dataSet.CI.edgeAutoCorr]                 = correlationBootstrap(totalEdgeACF, repmat( 2/sqrt(size(motionMap,2)), 1, size(totalCCF,2) ) );
[dataSet.meanValue.signalAutoCorr,dataSet.CI.signalAutoCorr]             = correlationBootstrap(totalSignalACF, repmat( 2/sqrt(size(motionMap,2)), 1, size(totalCCF,2) ) );
[dataSet.protrusion.meanValue.crossCorr,dataSet.protrusion.CI.crossCorr] = correlationBootstrap(protrusionCCF, repmat( 2/sqrt(size(motionMap,2)), 1, size(totalCCF,2) ) );
[dataSet.retraction.meanValue.crossCorr,dataSet.retraction.CI.crossCorr] = correlationBootstrap(retractionCCF, repmat( 2/sqrt(size(motionMap,2)), 1, size(totalCCF,2) ) );

savingMovieResultsPerCell(ML,cellData,'correlationEstimation','correlation')
savingMovieDataSetResults(ML,dataSet,'correlationEstimation','correlation')

end


function out = internalGetCorrelation(protrusion,activity,maxLag,acFlag)

[muCcf,muCCci,lags,xCorr] = getAverageCorrelation(protrusion,activity,'maxLag',maxLag);

out.total.crossCorr       = xCorr;
out.total.lag             = lags;

out.meanValue.crossCorr   = muCcf;
out.meanValue.lag         = lags;

out.CI.crossCorr          = muCCci;
out.CI.lag                = lags;

if acFlag
    
    [muProtAcf,protCI,~,protAcf] = getAverageCorrelation(protrusion,'maxLag',maxLag);
    [muSignAcf,signCI,~,signAcf] = getAverageCorrelation(activity,'maxLag',maxLag);
    
    out.total.edgeAutoCorr       = protAcf;
    out.total.signalAutoCorr     = signAcf;
    
    out.meanValue.edgeAutoCorr   = muProtAcf;
    out.meanValue.signalAutoCorr = muSignAcf;
    
    out.CI.edgeAutoCorr          = protCI;
    out.CI.signalAutoCorr        = signCI;
    
end

end