function [cellData,dataSet] = edgeVelocityQuantification(movieObj,varargin)
% This function quantifies:
%                           mean protrusion/retraction instantaneous values and confidence interval;
%                           mean protrusion/retraction persistence time and confidence interval;
%                           cluster velocities and persistence time
%                           for each window: persistence time;protrusion/retraction block;Upper and lower noise threshold
%
% All the above quantifications are done for each cell of the movieObj (movieList or movieData)
% If the input 'interval' is set, all the described quantifications are calculated for each time interval
%
% IMPORTANT: This function requires to format the time series first. To do that, run the function "formatEdgeVelocity".
%
%Usage: [cellData,dataSet] = edgeVelocityQuantification(movieObj,varargin)
%
% Input:
%       nBoot - # of boostrap samples to be used (default value 1000)
%
%       alpha - alpha used to generate the bootstrap confidence intervals
%         (default value 0.05)
%
%       interval - cell array with the time intervals in frames
%                   Ex : interval = {[1:30],[20:40]} - Quantification will be done with velocities calculated at each element of the cell array
%
%       scale         - convert velocity from pixel/frame into nm/sec.Default:false
%
%       includeWin    - cell array. Each element is a vector with the window indexes to be used in the analysis. Enter [] to analize all windows.
%
%       outLevel      - scalar for outlier detection(see detectOutlier)
%
%       trendType     - scalar - algorihtm used to detrend the data(See getTimeSeriesTrend)
%
%       minLength     - scalar - mininum time serie length
%
%       gapSize       - scalar - length of the nan gap to be closed. Default:0
%
%
%       outputPath    - string. Name of the output folder. Default: 'EdgeVelocityQuantification'
%
%       fileName      - string. Name of the output file. Default: 'edgeVelocity'
%
% Output:
%       cellData - this is a long structure array cellData(for each cell) with the following fields:
%                .data.excludeWin - indexes of windows that were excluded from analysis . Ex: border windows
%                     .rawEdgeMotion - raw edge velocity time series. Protrusion map. This is redundant.
%                     .procEdgeMotion - pre-processed edge velocity time series. Ex: trend, mean and NaN removed.
%
%                .protrusionAnalysis(for each interval).meanValue.persTime
%                                                                .maxVeloc
%                                                                .minVeloc
%                                                                .meanVeloc
%                                                                .mednVeloc
%                                                                .meanValue
%                                                                .cluster.persTime
%                                                                        .maxVeloc
%                                                                        .minVeloc
%                                                                        .meanVeloc
%                                                                        .mednVeloc
%                                                                        .meanValue
%
%           The structure continues with the confidence interval for each of the measurement above
%
%                                                      CI.persTimeCI
%                                                        .maxVelocCI
%                                                        .minVelocCI
%                                                        .meanVelocCI
%                                                        .mednVelocCI
%                                                        .meanValueCI
%                                                        .cluster.persTime
%                                                                .maxVeloc
%                                                                .minVeloc
%                                                                .meanVeloc
%                                                                .mednVeloc
%                                                                .meanValue
%
%                                                       Analysis at the individual time series level
%
%                                                      .windows(for each window).limit
%                                                                               .PersTime
%                                                                               .BlockOut
%                                                                               .MaxVeloc
%                                                                               .MeanVeloc
%                                                                               .MinVeloc
%                                                                               .MednVeloc
%
%Marco Vilela, 2012


ip = inputParser;
ip.addRequired('movieObj',@(x) isa(x,'MovieList') || isa(x,'MovieData'));
ip.addParamValue('nBoot',1e3,@isscalar);
ip.addParamValue('alpha',.05,@isscalar);
ip.addParamValue('cluster',false,@isscalar);
ip.addParamValue('nCluster',2,@isscalar);

if isa(movieObj,'MovieData')
    
    ML = movieData2movieList(movieObj);
    
else
    
    ML = movieObj;
    
end

nCell = numel(ML.movies_);

%% Time Series Pre-Processing operations
ip.addParamValue('includeWin',cell(1,nCell),@iscell);
ip.addParamValue('winInterval',num2cell(cell(1,nCell)),@iscell);
ip.addParamValue('outLevel',  zeros(1,nCell),@isvector);
ip.addParamValue('trendType',-ones(1,nCell),@isvector);
ip.addParamValue('minLength', 30*ones(1,nCell),@isvector);
ip.addParamValue('gapSize',   zeros(1,nCell),@isvector);
ip.addParamValue('scale',     false,@islogical);
ip.addParamValue('outputPath','EdgeVelocityQuantification',@isstr);
ip.addParamValue('fileName','edgeVelocity',@isstr);
ip.addParamValue('interval',num2cell(cell(1,nCell)),@iscell);

ip.parse(movieObj,varargin{:});
nBoot       = ip.Results.nBoot;
alpha       = ip.Results.alpha;
cluster     = ip.Results.cluster;
nCluster    = ip.Results.nCluster;
interval    = ip.Results.interval;
scale       = ip.Results.scale;
outputPath  = ip.Results.outputPath;
fileName    = ip.Results.fileName;

includeWin  = ip.Results.includeWin;
outLevel    = ip.Results.outLevel;
minLen      = ip.Results.minLength;
trend       = ip.Results.trendType;
gapSize     = ip.Results.gapSize;
winInterval = ip.Results.winInterval;

if numel(interval) == 1
   interval = repmat(interval,1,nCell);
end

%% Formatting Time Series
operations = {'includeWin',includeWin,'winInterval',winInterval,'outLevel',outLevel,'minLength',minLen,'trendType',trend,'gapSize',gapSize,'saveOn',true,'outputPath',outputPath,'fileName',fileName};
cellData   = formatMovieListTimeSeriesProcess(ML,'ProtrusionSamplingProcess',operations{:});
scaling    = 1;

for iCell = 1:nCell
    
    
    currMD  = ML.movies_{iCell};
     
    if ~isfield(cellData{iCell}.data,'scaling')
        cellData{iCell}.data.scaling = false;
    end
        
        
    if xor(scale,cellData{iCell}.data.scaling)
        
        if scale
            
            if isempty(currMD.pixelSize_) || isempty(currMD.timeInterval_)
                error(['Movie ' num2str(iCell) 'does not have the pixel size and/or time interval setup'])
            end
            
            scaling = (currMD.pixelSize_/currMD.timeInterval_);
        
        else
            
            scaling = (currMD.timeInterval_/currMD.pixelSize_);
            
        end
        
    end
    
    cellData{iCell}.data.scaling           = scale;
    cellData{iCell}.data.pixelSize         = currMD.pixelSize_;
    cellData{iCell}.data.frameRate         = currMD.timeInterval_;
    cellData{iCell}.data.nFrames           = currMD.nFrames_;
    cellData{iCell}.data.rawEdgeMotion     = cellData{iCell}.data.rawTimeSeries.*scaling;
    cellData{iCell}.data.procEdgeMotion    = cellData{iCell}.data.procTimeSeries.*scaling;
    cellData{iCell}.data.procExcEdgeMotion = cellfun(@(x) x*scaling, cellData{iCell}.data.procExcTimeSeries{1},'Unif',0);% 1 means first layer
    
    if ~isfield(cellData{iCell}.data,'interval')
        
        cellData{iCell}.data.interval = {[]};
        
    end
    
    if isfield(cellData{iCell}.data,'rawTimeSeries')
        
        cellData{iCell}.data = rmfield(cellData{iCell}.data,{'rawTimeSeries','procTimeSeries','procExcTimeSeries'});
        
    end
    
end


%% Getting Average Velocities and Persistence Time per Cell

%indexes for non-processed cells
flag1 = cellfun(@(x) isfield(x,'protrusionAnalysis'),cellData);

caseFlag = true;
if isempty(winInterval{1}{1})
    
    %If the interval is empty
    empInterval           = cell2mat( cellfun(@(x) isempty(x{1}),interval,'Unif',0) );
    interval(empInterval) = cellfun(@(x) {1:x.data.nFrames},cellData(empInterval),'Unif',0);
  
    %If number of the interval is different
    diffInter             = true(1,nCell);
    diffInter(flag1)      = cell2mat(cellfun(@(x,y) abs(numel(x.data.interval)-numel(y)),cellData(flag1),interval(flag1),'Unif',0)) ~= 0;
    %Same number of interval but the interval is different
    dWinInter             = true(1,nCell);
    dWinInter(~diffInter) = cell2mat(cellfun(@(x,y)  all(cell2mat(cellfun(@(w,z) isequaln(w,z),x.data.interval,y,'Unif',0))),cellData(~diffInter),interval(~diffInter),'Unif',0));
    flag2                 = diffInter | ~dWinInter;
    finalIdx              = ~flag1 | flag2;
    caseFlag              = false;
    
else
    
    finalIdx = ~flag1;
    
end


commonGround = @(x,z,y) mergingEdgeResults(x,'cluster',cluster,'nCluster',nCluster,'alpha',alpha,'nBoot',nBoot,'deltaT',z,'winInterval',y);
protrusion   = [];
retraction   = [];

if sum(finalIdx) ~= 0
    
    
    if caseFlag
        
        [protrusion,retraction] ...
                   = cellfun(@(x,y) commonGround(x.data.procExcEdgeMotion,x.data.frameRate,y),cellData(finalIdx),winInterval(finalIdx),'Unif',0);
        protrusion = num2cell(protrusion);
        retraction = num2cell(retraction);
        
    else
        mIdx           = flag1 & flag2;
        cellData(mIdx) = cellfun(@(x) rmfield(x,{'protrusionAnalysis','retractionAnalysis'}),cellData(mIdx),'Unif',0);
        firstLevel     = @(x,y,z) commonGround( cellfun(@(w) w(x),y,'Unif',0), z, {[]});
        secondLevel    = @(x,y,z) cellfun(@(w) firstLevel(w,y,z),x,'Unif',0);
        
        [protrusion,retraction] ...
                    = cellfun(@(x,y) secondLevel(y,x.data.procExcEdgeMotion,x.data.frameRate),cellData(finalIdx),interval(finalIdx),'Unif',0);
        
    end
    
end


[cellData,dataSet] = getDataSetAverage(cellData,protrusion,retraction,interval,alpha,nBoot,finalIdx);

for iCell = find(finalIdx)
    
    if isempty(interval{1})
        cellData{iCell}.data.interval = {1:cellData{iCell}.data.nFrames};
    else
        cellData{iCell}.data.interval =interval{iCell};
    end
    
end

%% Saving results

savingMovieResultsPerCell(ML,cellData,outputPath,fileName)
savingMovieDataSetResults(ML,dataSet,outputPath,fileName)

end%End of main function

function [cellData,dataSet] = getDataSetAverage(cellData,protrusion,retraction,interval,alpha,nBoot,idx)
%This function pull all the data from individual cells and calculates the dataSet mean value and CI
%It also formats the data structure for plotting
%Input:
%       cellData   - structure created by the function "formatEdgeVelocity.m".
%       interval   - Cell array containing the time intervals where the analysis is performed
%       alpha      - confidence interval Ex: - 0.05 = 95 percent
%       nBoot      - number of bootstrap samples
%
%Output:
%       dataSet    - averages and CI for the data set
%
nCell = numel(cellData);
cc    = 1;

for iCell = find(idx)
    
            
        for iiInt = 1:numel(interval{iCell})
            cellData{iCell}.protrusionAnalysis(iiInt) = protrusion{cc}{iiInt};
            cellData{iCell}.retractionAnalysis(iiInt) = retraction{cc}{iiInt};
        end
        cc = cc + 1;
        
end

dataSet = [];
nInter  = cellfun(@(x) numel(x),interval);

if sum(rem(nInter,nInter(1))) == 0
    
    for iInt = 1:numel(interval{1})
        
        total = struct('ProtPersTime',[],'ProtMaxVeloc',[],'ProtMinVeloc',[],'ProtMeanVeloc',[],'ProtMednVeloc',[],'RetrPersTime',[],'RetrMaxVeloc',[],'RetrMinVeloc',[],'RetrMeanVeloc',[],'RetrMednVeloc',[]);
        
        for iCell = 1:nCell
            
            total.ProtPersTime    = [total.ProtPersTime;cellData{iCell}.protrusionAnalysis(iInt).total.persTime];
            total.ProtMaxVeloc    = [total.ProtMaxVeloc;cellData{iCell}.protrusionAnalysis(iInt).total.maxVeloc];
            total.ProtMinVeloc    = [total.ProtMinVeloc;cellData{iCell}.protrusionAnalysis(iInt).total.minVeloc];
            total.ProtMeanVeloc   = [total.ProtMeanVeloc;cellData{iCell}.protrusionAnalysis(iInt).total.Veloc];
            total.ProtMednVeloc   = [total.ProtMednVeloc;cellData{iCell}.protrusionAnalysis(iInt).total.mednVeloc];
            
            total.RetrPersTime    = [total.RetrPersTime;cellData{iCell}.retractionAnalysis(iInt).total.persTime];
            total.RetrMaxVeloc    = [total.RetrMaxVeloc;cellData{iCell}.retractionAnalysis(iInt).total.maxVeloc];
            total.RetrMinVeloc    = [total.RetrMinVeloc;cellData{iCell}.retractionAnalysis(iInt).total.minVeloc];
            total.RetrMeanVeloc   = [total.RetrMeanVeloc;cellData{iCell}.retractionAnalysis(iInt).total.Veloc];
            total.RetrMednVeloc   = [total.RetrMednVeloc;cellData{iCell}.retractionAnalysis(iInt).total.mednVeloc];
            
        end
        
        [dataSet.CI.interval(iInt),dataSet.meanValue.interval(iInt)] = structfun(@(x) bootStrapMean(x,alpha,nBoot),total,'Unif',0);
        dataSet.total.interval(iInt) = total;
        
    end
    
end

end