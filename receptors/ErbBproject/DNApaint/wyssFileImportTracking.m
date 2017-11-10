function [ success ] = wyssFileImportTracking(dir)
%Reads in wyss data and then applies global tracking
%   dir is the base directory containing .mat files for one image from the
%   wyss.

success =1;

ip = inputParser;

ip.addRequired('dir',@ischar);
ip.parse(dir);

PointList = wyssFileImport(dir);
wyssFileTracking(PointList);


end

