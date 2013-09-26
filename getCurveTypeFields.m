function [curvTypes,curvNames,curvUnits,curvConv] = getCurveTypeFields(pixelSizeNM,useMicrons)
%GETCURVETYPEFIELDS returns all the curvature measures and associated parameters
%
% [curvTypes,curvNames,curvUnits,curvConv] = getCurveTypeFields(pixelSizeNM,useMicrons)
%
% Figure it out.
%

% Hunter Elliott
% 3/2013

%Quick and dirty way to have all the curv fields available to multiple
%functions.

curvTypes = {'gaussCurvSampMean',...             
             'meanCurvSampMean',...             
             'PC1CurvSampMean',...             
             'PC2CurvSampMean',...                         
             'MaxAbsPCCurvSampMean'};

%And intuitive names for them...
curvNames = {'Gaussian Curvature',...             
             'Mean Curvature',...             
             'k2',...             %Due to sign switch, we change the names here.
             'k1',...             
             'Max Absolute Curvature'};
         
if nargin < 2 || isempty(useMicrons)
    useMicrons = false;
end

if useMicrons
    uStr = 'um';
    pixelSize = pixelSizeNM / 1e3;
else
    uStr = 'nm';
    pixelSize = pixelSizeNM;
end

if nargin > 0 && ~isempty(pixelSize)         
    curvUnits = cell(numel(curvNames),1);
    curvUnits{1} = ['1/' uStr '^2'];
    [curvUnits{2:end}] = deal(['1/' uStr '']);            
    curvConv = nan(numel(curvNames),1);
    curvConv(1) = 1/pixelSize^2;
    curvConv(2:4) = -1/pixelSize;
    curvConv(5) = 1/pixelSize;
else
    curvUnits = cell(numel(curvNames),1);
    curvUnits{1} = '1/pixels^2';
    [curvUnits{2:end}] = deal('1/pixels');
    curvConv = nan(numel(curvNames),1);
    curvConv(1) = 1;
    curvConv(2:4) = -1;    
    curvConv(5) = 1;
    
end