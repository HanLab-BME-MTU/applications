function [maxout, minout, meanout, medianout, stdout] = summary(array)

maxout    = max(array);
minout    = min(array);
meanout   = mean(array);
medianout = median(array);
stdout    = std(array);