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
%       cluster - scalar "1" to perform cluster analysis; "0" otherwise       
%
%       nCluster - number of cluster 
%
%       interval - cell array with the time intervals in frames
%                   Ex : interval = {[1:30],[20:40]} - Quantification will be done with velocities calculated at each element of the cell array
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
ip.addParamValue('interval',{[]},@iscell);

ip.parse(movieObj,varargin{:});
nBoot    = ip.Results.nBoot;
alpha    = ip.Results.alpha;
cluster  = ip.Results.cluster;
nCluster = ip.Results.nCluster;
interval = ip.Results.interval;

if isa(movieObj,'MovieData')
    
    ML = movieData2movieList(movieObj);
    
else
    
    ML = movieObj;
    
end


cellData = loadingMovieResultsPerCell(ML);

%% Getting Average Velocities and Persistence Time per Cell

commonGround = @(x,z) mergingEdgeResults(x,'cluster',cluster,'nCluster',nCluster,'alpha',alpha,'nBoot',nBoot,'deltaT',z);
if isempty(interval{1})
    
    [protrusionA,retractionA] ...
                = arrayfun(@(x) commonGround(x.data.excProcEdgeMotion,x.data.frameRate),cellData,'Unif',0);
     protrusion = cellfun(@(x) {x},protrusionA,'Unif',0);
     retraction = cellfun(@(x) {x},retractionA,'Unif',0);
else
    
    firstLevel  = @(x,y,z) commonGround( cellfun(@(w) w(x),y,'Unif',0), z);
    secondLevel = @(x,y,z) cellfun(@(w) firstLevel(w,y,z),x,'Unif',0);
    
    [protrusion,retraction] ...
                = arrayfun(@(x) secondLevel(interval,x.data.excProcEdgeMotion,x.data.frameRate),cellData,'Unif',0);
    
end

[cellData,dataSet] = getDataSetAverage(cellData,protrusion,retraction,interval,alpha,nBoot);

%% Saving results
savingMovieResultsPerCell(ML,cellData)
savingMovieDataSetResults(ML,dataSet)

end%End of main function

function [cellData,dataSet] = getDataSetAverage(cellData,protrusion,retraction,interval,alpha,nBoot)
%This function pull all the data from individual cells and calculates the dataSet mean value and CI
%It also formats the data structure for plotting
%Input:
%       cellData   - structure created by the function "formatEdgeVelocity.m". 
%       interval   - Cell array containing the time intervals where the analysis is performed
%       alpha      - confidence interval Ex: - 0.05 = 95 percent
%       nBoot      - number of bootstrap samples
%
%Output:
%       dataSet         - averages and CI for the data set
%
nCell = numel(cellData);
total = struct('ProtPersTime',[],'ProtMaxVeloc',[],'ProtMinVeloc',[],'ProtMeanVeloc',[],'ProtMednVeloc',[],'RetrPersTime',[],'RetrMaxVeloc',[],'RetrMinVeloc',[],'RetrMeanVeloc',[],'RetrMednVeloc',[]);

for iInt = 1:numel(interval)
    
    for iCell = 1:nCell
        
        cellData(iCell).protrusionAnalysis(iInt) = protrusion{iCell}{iInt};
        cellData(iCell).retractionAnalysis(iInt) = retraction{iCell}{iInt};
        
        timeScale             = cellData(iCell).data.frameRate;
        
        total.ProtPersTime    = [total.ProtPersTime;cellData(iCell).protrusionAnalysis(iInt).total.persTime*timeScale];
        total.ProtMaxVeloc    = [total.ProtMaxVeloc;cellData(iCell).protrusionAnalysis(iInt).total.maxVeloc];
        total.ProtMinVeloc    = [total.ProtMinVeloc;cellData(iCell).protrusionAnalysis(iInt).total.minVeloc];
        total.ProtMeanVeloc   = [total.ProtMeanVeloc;cellData(iCell).protrusionAnalysis(iInt).total.Veloc];
        total.ProtMednVeloc   = [total.ProtMednVeloc;cellData(iCell).protrusionAnalysis(iInt).total.mednVeloc];
        
        total.RetrPersTime    = [total.RetrPersTime;cellData(iCell).retractionAnalysis(iInt).total.persTime*timeScale];
        total.RetrMaxVeloc    = [total.RetrMaxVeloc;cellData(iCell).retractionAnalysis(iInt).total.maxVeloc];
        total.RetrMinVeloc    = [total.RetrMinVeloc;cellData(iCell).retractionAnalysis(iInt).total.minVeloc];
        total.RetrMeanVeloc   = [total.RetrMeanVeloc;cellData(iCell).retractionAnalysis(iInt).total.Veloc];
        total.RetrMednVeloc   = [total.RetrMednVeloc;cellData(iCell).retractionAnalysis(iInt).total.mednVeloc];
        
    end
    
    
    [dataSet.CI.cond(iInt),dataSet.meanValue.cond(iInt)] = structfun(@(x) bootStrapMean(x,alpha,nBoot),total,'Unif',0);
    dataSet.total.cond(iInt) = total;
    
end
    
end