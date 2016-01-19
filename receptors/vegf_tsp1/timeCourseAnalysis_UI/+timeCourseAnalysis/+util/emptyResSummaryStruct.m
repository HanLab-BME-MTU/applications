function [ s ] = emptyResSummaryStruct( )
%emptyResSummaryStruct Returns an empty res summary struct

s = struct('diffSummary',[],'diffCoefMeanPerClass',[],...
        'confRadMeanPerClass',[],'ampMeanPerClass',[],'ampStatsF20',[],...
        'ampStatsL20',[],'statsMS',[],'msTimeInfo',[],'cellArea',[]);

end

