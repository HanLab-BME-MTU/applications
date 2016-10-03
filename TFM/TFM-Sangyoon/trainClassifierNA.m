function [trainedClassifier, validationAccuracy,C,order,validationPredictions, validationScores] = trainClassifierNA(datasetTable)
% Extract predictors and response
predictorNames = {'decayingIntensityNAs', 'edgeAdvanceSpeedNAs', 'advanceSpeedNAs', 'lifeTimeNAs', 'meanIntensityNAs', 'distToEdgeFirstNAs', 'startingIntensityNAs', 'distToEdgeChangeNAs', 'distToEdgeLastNAs', 'edgeAdvanceDistFirstChangeNAs', 'edgeAdvanceDistLastChangeNAs', 'maxEdgeAdvanceDistChangeNAs'};
predictors = datasetTable(:,predictorNames);
predictors = table2array(varfun(@double, predictors));
response = datasetTable.Group;
% Get the unique resonses
totalGroups = unique(response);
% Train a classifier
template = templateSVM('KernelFunction', 'polynomial', 'PolynomialOrder', 2, 'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1);
trainedClassifier = fitcecoc(predictors, response,'FitPosterior',1, 'Learners', template, 'Coding', 'onevsone', 'PredictorNames', {'decayingIntensityNAs' 'edgeAdvanceSpeedNAs' 'advanceSpeedNAs' 'lifeTimeNAs' 'meanIntensityNAs' 'distToEdgeFirstNAs' 'startingIntensityNAs' 'distToEdgeChangeNAs' 'distToEdgeLastNAs' 'edgeAdvanceDistFirstChangeNAs' 'edgeAdvanceDistLastChangeNAs' 'maxEdgeAdvanceDistChangeNAs'}, 'ResponseName', 'Group', 'ClassNames', totalGroups');

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
