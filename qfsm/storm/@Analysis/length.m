function avLen = length(obj)
% Compute average model length
avLen = mean(obj.data.modelLength);

% Display histogram
nBins = round(sqrt(size(obj.data.modelLength,1)));
hist(obj.data.modelLength,nBins);
xlabel('Model length [nm]');
ylabel('Number of models');
end