function meanIntensity = intensity(obj)
meanIntensity = mean(obj.data.intensity);
fprintf('Analyzer: Mean intensity: %f\n',meanIntensity);

% Display intensity histogram
hist(obj.data.intensity,round(sqrt(numel(obj.data.intensity))));
end