function [ figureData, commonInfo ] = getFigureData( matFile )
%getFigureData Obtain figure data from matFile

if(nargin < 1 || isempty(matFile))
    matFile = [pwd filesep 'figureData.mat'];
elseif(exist(matFile,'dir'))
    matFile = [matFile filesep 'figureData.mat'];
end

try
    S = load(matFile);
    figureData = S.figureData;
    commonInfo = S.commonInfo;
catch
    error('timeCourseAnalysis:plot:getFigureData:MatFileLoadFailed', ...
        'Could not load matFile: %s',matFile);
end


end

