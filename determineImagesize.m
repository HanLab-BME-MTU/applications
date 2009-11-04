function [data] = determineImagesize(data)

for i = 1:length(data)
    tifFiles = dir([data(i).source filesep '*.tif']);
    data(i).imagesize = size(imread([data(i).source filesep tifFiles(1).name]));
end