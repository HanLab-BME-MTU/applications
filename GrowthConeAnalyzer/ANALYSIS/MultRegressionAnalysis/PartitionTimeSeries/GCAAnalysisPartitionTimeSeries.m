function [ localMeasGrouped,predVarMat,responseMat ] = GCAAnalysisPartitionTimeSeries(localMeas,globalMeas,varargin)
%GCAAnalysisPartitionTimeSeries
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%  localMeas:    REQUIRED:      structure with parameter fields, number of
%                                 fields will be dependent upon the what
%                                 the user has run / chose to analyze in
%                                 the previous step. Each field is
%                                 contains a rxc double array where
%                                 N is the number of observations (padded
%                                 with NaN and M is the number of time
%                                 frames). Note I chose this format for now
%                                 instead of a cell format as it allows for
%                                 easy feeding into both a .csv/excel file
%                                 and/or the matlab boxplot function which
%                                 facilitates plotting different data
%                                 groupings.
% grouping:        REQUIRED:      nx1 double array providing the group ID
%                                 where n is the number of time points of
%                                 the movie. So far this grouping is
%                                 currently the output of
%                                 GCAfindPausingInNeuriteTrajectory
%
% writePredMat:      PARAM:       logical: flag to output an rxc predictor
%                                 matrix
%                                 where r (row) is the number of groups in the time
%                                 series and c (col) is the number of params (or predictor variables)
%                                 DEFAULT: True
%
% predFunc PARAM:     function:  The function used to
%                                to obtain the final predictor variables
%                                opertation is performed on each group
%                                cluster for each time series (still trying
%                                to decide if breaking up the time series
%                                in this manner is smart) maybe could
%                                simply used all the points of the time
%                                series and see the difference?- for now we
%                                will keep as it
%                                DEFAULT: @nanmean
%
%  secPerFrame          PARAM:    scalar: time interval per frame in sec
%                                 DEFAULT: 5 sec
%
%  makePlot:            PARAM:    logical: flag for sanity plots
%                                 DEFAULT: true
%
%
%  outPath:             PARAM:    character or empty: if is empty do nothing
%                                 else save results (including any figures)
%                                 in directory specfied by outPath
%                                 DEFAULT: [];
%
% OUTPUT:
% localMeasGrouped             Structure with N Fields where N is the number ... 
%                                  of local measurements collected from the descriptor folders 
%                                  .measurementName 
%                                 where each field holds an 1xG number of cells 
%                                 where G is the number of groups defined in globalMeas.   
%                                 Each cell contains a rxc double array where
%                                 N is the number of observations (padded
%                                 with NaN and M is the number of time
%                                 frames in that partitioned group.
%
%
%
%
%
%% Check Input
ip = inputParser;
% REQUIRED
ip.addRequired('localMeas',@(x) isstruct(x));
ip.addRequired('globalMeas',@(x) isstruct(x));

% PARAMS
ip.addParameter('writePredMat',true,@(x) islogical(x));
ip.addParameter('predFunc',@nanmedian, @(x) isa(x, 'function_handle')); % so this is the tricky
% bit about how I should format these predictors... currently pooling over
% all the frames of the trajectory group and then just taking a mean/or a
% median value. 
ip.parse(localMeas,globalMeas,varargin{:});

writePredMat = ip.Results.writePredMat;
predFunc = ip.Results.predFunc;

%% Initiate
grouping = globalMeas.outgrowth.grouping;
params = fieldnames(localMeas);
localParamGroups.grouping = grouping;
nParams = numel(params);
nGroups = length(unique(grouping));

if writePredMat == true;
    % initiate the predVarMat
    predVarMat = nan(nGroups,nParams);
   
end

%% Cluster

for iParam = 1:numel(params)
    groupingF = grouping(1:end-1); 
    values = localMeas.(params{iParam}); % remember each set is of the form rxc
    % such that r is the number of obs in a given frame
    numGroups = max(unique(grouping));
    % matrices of values for each distribution maintaining frame info per
    % column 
    
    % test to make sure the grouping var and the descriptor per frame match
    % (they should always but I was sloppy with the number of frames I
    % used)
    if size(values,2) > length(groupingF) % measured an extra frame at the end
        % trunc values
        values = values(:,1:length(groupingF));
        display(['Number frames calculated for ' params{iParam} 'is larger than the '...
            'neurite length calculation!: Assuming Truncation at End of Movie']);
    elseif size(values,2) < length(groupingF) 
        % trunc the grouping var
         groupingF = groupingF(1:size(values,2));
        
      
        display(['Number of frames calculated for ' params{iParam} 'is shorter than the '...
            'neurite length calculation! : Assuming Truncation of the End of Movie']); 
        
    end
    
    
    valsClust = arrayfun(@(i) values(:,groupingF==i),1:numGroups,'uniformoutput',0);
      
    % add field
    localMeasGrouped.(params{iParam}) = valsClust;
    if writePredMat == true; % if write the predictor
        % put into a matrix where r (row) is the number of observations -
        % here the number of groups in the time series
        % and c (col) is the predictor variable
        
        % pool the distribution over the relevent frames 
        valsClustPooled = cellfun(@(x) x(:),valsClust,'uniformoutput',0); % for now just pool the data over the frames in the cluster based on outgrowth.
        % remove all extra NaNs (just in case user specifies a function
        % that requires these be removed. 
        valsClustPooled = cellfun(@(x) x(~isnan(x)),valsClustPooled,'uniformoutput',0); 
        
        
        % peform the predictor operation (for example averaging)
        predVarMat(:,iParam) = cellfun(predFunc,valsClustPooled);
        
    end
   
    
end
 responseMat = cellfun(predFunc,globalMeas.outgrowth.groupedVelSmoothed)';

