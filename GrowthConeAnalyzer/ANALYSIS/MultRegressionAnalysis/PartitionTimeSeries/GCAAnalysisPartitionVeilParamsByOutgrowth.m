function [ localParamsGrouped,predVarMat ] = GCAAnalysisPartitionVeilParamsByOutgrowth( analysisResults1,analysisResults2,globalParams,varargin)
%GCAAnalysisPartitionVeilParamsByOutgrowth
% a small function to reformat the analysisResults from Marcos veil
% protrusion/retraction identifier and cluster by
% INPUT:
%       analysisResults:   REQUIRED:   very complex structure that is the
%                                      output of Marcos edge
%                                      protrusion/retraction event detector
%                                      within it contains
%                                      protrusion/retraction measurements
%                                      for a number of individual
%                                      protrusion/retraction events
%
%      grpFrames:          REQUIRED:   a cell array of the frame numbers
%                                      for each outgrowth group
%
%OUTPUT:
% localParamsGrouped                     :  a structure of the parameters
%%
%% Check Input
ip = inputParser;
% REQUIRED
ip.addRequired('analysisResults1',@(x) isstruct(x));
ip.addRequired('analysisResults2',@(x) isstruct(x));
ip.addRequired('globalParams',@(x) isstruct(x));

% PARAMS
ip.addParamValue('writePredMat',true,@(x) islogical(x));
ip.addParamValue('predFunc',@nanmedian, @(x) isa(x, 'function_handle')); % so this is the tricky
% bit about how I should format these predictors... currently pooling over
% all the frames of the trajectory group and then just taking a mean/or a
% median value.
ip.parse(analysisResults1,analysisResults2,globalParams,varargin{:});



writePredMat = ip.Results.writePredMat;
predFunc = ip.Results.predFunc;

%% Initiate
grpFrames = globalParams.outgrowth.groupedFrames;
nGroups = numel(grpFrames); 
% Define what parameters you would like to extract from marcos data (can
% make input)
params{1} = 'Veloc';
params{2} = 'persTime';

analysisC{1} = 'protrusionAnalysis';
analysisC{2} = 'retractionAnalysis';

if writePredMat == true;
    % initiate the predVarMat
    predVarMat = nan(nGroups,4);
end
count = 1;
%% Initiate Loop to collect clusters
for iAnal = 1:numel(analysisC)
    for iParam = 1: numel(params)
    
        paramName = [analysisC{iAnal} '_' params{iParam}];
        
        % initiate valuesClust to hold the veil measurements per outgrowth
        % group
        valuesClust = cell(numel(grpFrames),1);
        
        for iGroup = 1:numel(grpFrames);
            
            grpC= grpFrames{iGroup};
            
            % compile the time frames for all the protrusion or retraction events from all windows
            framesOfVeilEvent1 = horzcat(analysisResults1.(analysisC{iAnal}).windows(:).blockOut);
            
            % compile the time frames for all the protrusion or retraction events from all windows
            framesOfVeilEvent2 = horzcat(analysisResults2.(analysisC{iAnal}).windows(:).blockOut);
            
            % note it looks like the block out is not in absolute number of
            % frames but idx convert the second half of hte movie so that
            % it is frames. 
            framesOfVeilEvent2 = cellfun(@(x) (x + 61), framesOfVeilEvent2,'uniformoutput',0); 
            framesOfVeilEvent = [framesOfVeilEvent1 framesOfVeilEvent2]; 
            
            % collect all the measurements for the current protrusion/retraction
            % event
            valuesC1 = vertcat(analysisResults1.(analysisC{iAnal}).windows(:).(params{iParam}));
            valuesC2 = vertcat(analysisResults2.(analysisC{iAnal}).windows(:).(params{iParam}));
            valuesC = [valuesC1;valuesC2]; 
            
            % find all the measurements that lie completely within the
            % outgrowth event timescale.
            idxC  = cellfun(@(x) length(intersect(grpC,x)) == length(x),framesOfVeilEvent);
            valuesClust{iGroup} = valuesC(idxC);
            
        end % for iGroup
        
        % save the parameter
        localParamsGrouped.(paramName) = valuesClust;
        
        if writePredMat == true; % if write the predictor
            % put into a matrix where r (row) is the number of observations -
            % here the number of groups in the time series
            % and c (col) is the predictor variable
            
            % peform the predictor operation (for example averaging)
            % rows observations in time series, columns are predictor vars
            predVarMat(:,count) = cellfun(predFunc,valuesClust);
            
        end
        
        
        count = 1+count;
    end %for iParam
    
    
    
    
end % for iAnal
end

