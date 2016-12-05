function  [validationAccuracy,C,order] = validateClassifier(trainedClassifier,datasetTable)
% Extract predictors and response
predictorNames = {'decayingIntensityNAs', 'edgeAdvanceSpeedNAs', 'advanceSpeedNAs', 'lifeTimeNAs', 'meanIntensityNAs', 'distToEdgeFirstNAs',...
    'startingIntensityNAs', 'distToEdgeChangeNAs', 'distToEdgeLastNAs', 'edgeAdvanceDistFirstChangeNAs', 'edgeAdvanceDistLastChangeNAs', 'maxEdgeAdvanceDistChangeNAs',...
    'maxIntensityNAs', 'timeToMaxInten', 'edgeVariation'};
predictors = datasetTable(:,predictorNames);
predictors = table2array(varfun(@double, predictors));
response = datasetTable.Group;
% Perform cross-validation
% figure; imagesc(predictors);
% predictedLabels = predict(trainedClassifier,predictors);
% [predictedLabels,NegLoss,PBScore] = trainedClassifier.predict(predictors);
predictedLabels = trainedClassifier.predict(predictors);
% if the number of response is different (e.g. no Group6), merge into
% Group9.
results = nan(1,numel(predictedLabels));
for i = 1 : numel(predictedLabels)
    results(i) = strcmp(predictedLabels{i},response{i});
end
validationAccuracy=sum(results)/length(results);
% confusion matrix
ascendingOrder = unique(sort([response; predictedLabels]));
[C,order] = confusionmat(response,predictedLabels,'order',ascendingOrder);
end