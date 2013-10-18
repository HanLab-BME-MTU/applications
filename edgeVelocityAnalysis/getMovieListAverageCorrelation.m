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
ip.addParamValue('winInterval',num2cell(cell(1,nCell)),@iscell);
ip.addParamValue('signalOutLevel',  zeros(1,nCell),@isvector);
ip.addParamValue('edgeOutLevel',  zeros(1,nCell),@isvector);
ip.addParamValue('trendType',-ones(1,nCell),@isvector);
ip.addParamValue('minLength', 30*ones(1,nCell),@isvector);
ip.addParamValue('gapSize',   zeros(1,nCell),@isvector);
ip.addParamValue('scale',false,@islogical);
ip.addParamValue('maxLag',0,@isscalar);
ip.addParamValue('layer',1,@isscalar);
ip.addParamValue('signalChannel',1,@isscalar);
ip.addParamValue('interval',{[]},@iscell);
ip.addParamValue('robust',false,@islogical);

ip.parse(movieObj,varargin{:});
scale            = ip.Results.scale;
includeWin       = ip.Results.includeWin;
winInterval      = ip.Results.winInterval;
signalOutLevel   = ip.Results.signalOutLevel;
edgeOutLevel     = ip.Results.edgeOutLevel;
minLen           = ip.Results.minLength;
trend            = ip.Results.trendType;
maxLag           = ip.Results.maxLag;
layer            = ip.Results.layer;
channel          = ip.Results.signalChannel;
interval         = ip.Results.interval;
robust           = ip.Results.robust;

cellData{1,nCell} = [];
edgeInputParam    = {'winInterval',winInterval,'outLevel',edgeOutLevel,'minLength',minLen,'trendType',trend,'includeWin',includeWin,'gapSize',ones(1,nCell),'scale',scale,'outputPath','correlationEstimation'};
signalInputParam  = {'winInterval',winInterval,'outLevel',signalOutLevel,'minLength',minLen,'trendType',trend,'includeWin',includeWin,'gapSize',ones(1,nCell),'outputPath','correlationEstimation'};
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
    
    windows         = intersect(edge{iCell}.data.includedWin{1},signal{iCell}.data.includedWin{layer});
    tsLengths       = cell2mat( cellfun(@(x) numel(x),edge{iCell}.data.winInterval,'Unif',0) );
    
    nWin            = numel(windows);
    motionMap       = nan(nWin,max(tsLengths));
    activity        = motionMap;
    
    for iWin = 1:nWin
        motionMap(iWin,1:tsLengths(iWin)) = edge{iCell}.data.procExcEdgeMotion{iWin};
        activity(iWin,1:tsLengths(iWin))  = signal{iCell}.data.procExcSignal{layer}{iWin};
    end
    
    cellData{iCell} = internalGetCorrelation(motionMap,activity,maxLag,1,robust);
    
    totalEdgeACF    = cat(2,cellData{iCell}.total.edgeAutoCorr,totalEdgeACF);
    totalSignalACF  = cat(2,cellData{iCell}.total.signalAutoCorr,totalSignalACF);
    totalCCF        = cat(2,cellData{iCell}.total.crossCorr,totalCCF);
    
    protIdx         = arrayfun(@(x) cell2mat(x.blockOut(:)),edge{iCell}.protrusionAnalysis.windows,'Unif',0);
    retrIdx         = arrayfun(@(x) cell2mat(x.blockOut(:)),edge{iCell}.retractionAnalysis.windows,'Unif',0);
    
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
    cellData{iCell}.total.retractionSignal = nanmean(sigRetraction(:));

    cellData{iCell}.protrusionCorrelation  = internalGetCorrelation(protrusion,activity,maxLag,0,robust);
    cellData{iCell}.retractionCorrelation  = internalGetCorrelation(retraction,activity,maxLag,0,robust);
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


function out = internalGetCorrelation(protrusion,activity,maxLag,acFlag,robust)

[muCcf,muCCci,lags,xCorr] = getAverageCorrelation(protrusion,activity,'maxLag',maxLag,'robust',robust);

out.total.crossCorr       = xCorr;
out.total.lag             = lags;

out.meanValue.crossCorr   = muCcf;
out.meanValue.lag         = lags;

out.CI.crossCorr          = muCCci;
out.CI.lag                = lags;

if acFlag
    
    [muProtAcf,protCI,~,protAcf] = getAverageCorrelation(protrusion,'maxLag',maxLag,'robust',robust);
    [muSignAcf,signCI,~,signAcf] = getAverageCorrelation(activity,'maxLag',maxLag,'robust',robust);
    
    out.total.edgeAutoCorr       = protAcf;
    out.total.signalAutoCorr     = signAcf;
    
    out.meanValue.edgeAutoCorr   = muProtAcf;
    out.meanValue.signalAutoCorr = muSignAcf;
    
    out.CI.edgeAutoCorr          = protCI;
    out.CI.signalAutoCorr        = signCI;
    
end

end