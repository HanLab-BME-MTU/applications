function [trainedClassifier, validationAccuracy, C, order, validationPredictions, validationScores] = trainClassifierKNN(datasetTable)
% Extract predictors and response
predictorNames = {'decayingIntensityNAs', 'edgeAdvanceSpeedNAs', 'advanceSpeedNAs', 'lifeTimeNAs', 'meanIntensityNAs', 'distToEdgeFirstNAs', 'startingIntensityNAs', 'distToEdgeChangeNAs', 'distToEdgeLastNAs', 'edgeAdvanceDistFirstChangeNAs', 'edgeAdvanceDistLastChangeNAs', 'maxEdgeAdvanceDistChangeNAs'};
predictors = datasetTable(:,predictorNames);
predictors = table2array(varfun(@double, predictors));
response = datasetTable.Group;
% Train a classifier
trainedClassifier = fitcknn(predictors, response, 'PredictorNames', {'decayingIntensityNAs' 'edgeAdvanceSpeedNAs' 'advanceSpeedNAs' 'lifeTimeNAs' 'meanIntensityNAs' 'distToEdgeFirstNAs' 'startingIntensityNAs' 'distToEdgeChangeNAs' 'distToEdgeLastNAs' 'edgeAdvanceDistFirstChangeNAs' 'edgeAdvanceDistLastChangeNAs' 'maxEdgeAdvanceDistChangeNAs'}, 'ResponseName', 'Group', 'ClassNames', {'Group1' 'Group2' 'Group3' 'Group4' 'Group5' 'Group7' 'Group8' 'Group9'}, 'Distance', 'Euclidean', 'Exponent', '', 'NumNeighbors', 10, 'DistanceWeight', 'SquaredInverse', 'StandardizeData', 1);

% Perform cross-validation
partitionedModel = crossval(trainedClassifier, 'KFold', 5);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

% confusion matrix
predictedLabels = trainedClassifier.predict(predictors);
[C,order] = confusionmat(response,predictedLabels);
%% Uncomment this section to compute validation predictions and scores:
% % Compute validation predictions and scores
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);