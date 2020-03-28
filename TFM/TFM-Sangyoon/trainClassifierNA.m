function [trainedClassifier, validationAccuracy,C,order,validationPredictions, validationScores] = trainClassifierNA(datasetTable)
% Extract predictors and response
predictorNames = {'decayingIntensityNAs', 'edgeAdvanceSpeedNAs', 'advanceSpeedNAs', 'lifeTimeNAs', 'meanIntensityNAs', 'distToEdgeFirstNAs', 'startingIntensityNAs',...
    'distToEdgeChangeNAs', 'distToEdgeLastNAs', 'edgeAdvanceDistFirstChangeNAs', 'edgeAdvanceDistLastChangeNAs', 'maxEdgeAdvanceDistChangeNAs',...
    'maxIntensityNAs', 'timeToMaxInten', 'edgeVariation', 'Area', 'FAfinishing', 'ampSlopeNAs', 'earlyAmpSlopeNAs', 'lateAmpSlopeNAs',...
    'startingAsNA','homogeneityAll'}; %, 'earlyAmpSlopeNAs', 'lateAmpSlopeNAs'};

missingVarNames=setdiff(predictorNames,datasetTable.Properties.VariableNames);
if ~isempty(missingVarNames)
    predictorNames(strcmp(predictorNames,missingVarNames))=[];
end
predictors = datasetTable(:,predictorNames);
predictors = table2array(varfun(@double, predictors));
response = datasetTable.Group;
% Get the unique resonses
[totalGroups, ~, ic] = unique(response);
% Filter out small training data group
arrayIc = unique(ic);
numEntitiesInGroup=arrayfun(@(x) sum(ic==x), arrayIc);
minNumGroup=5; bigEnoughGroups=[];
while isempty(bigEnoughGroups)
    bigEnoughGroups=find(numEntitiesInGroup>minNumGroup);
    minNumGroup=minNumGroup-1;
    if minNumGroup<1
        error('There is not enough training data! Cannot proceed with training')
    end
end

bigEnoughGroupsIcCellArray=arrayfun(@(x) (ic==x), bigEnoughGroups,'UniformOutput',false);
bigEnoughGroupsIc = bigEnoughGroupsIcCellArray{1};
for kk=2:numel(bigEnoughGroupsIcCellArray)
    bigEnoughGroupsIc = bigEnoughGroupsIc | bigEnoughGroupsIcCellArray{kk};
end
predictors = predictors(bigEnoughGroupsIc,:);
response = response(bigEnoughGroupsIc,:);
totalGroups = totalGroups(bigEnoughGroups);
% Train a classifier
template = templateSVM('KernelFunction', 'polynomial', 'PolynomialOrder', 2, 'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1);
trainedClassifier = fitcecoc(predictors, response,'FitPosterior',1, 'Learners', template, 'Coding', 'onevsone', 'PredictorNames', ...
    predictorNames, 'ResponseName', 'Group', 'ClassNames', totalGroups'); %, 'earlyAmpSlopeNAs', 'lateAmpSlopeNAs'

% Perform cross-validation
partitionedModel = crossval(trainedClassifier, 'KFold', 5);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

% confusion matrix
predictedLabels = trainedClassifier.predict(predictors);
[C,order] = confusionmat(response,predictedLabels);

%% Uncomment this section to compute validation predictions and scores:
% Compute validation predictions and scores
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);
end
