function [data] = fillStructSurvivalFunction_segment(data)
% fillStructSurvivalFunction fills in the values for the survival function
% based on the lifetime dat (lftHist_censored)
%
% SYNOPSIS [data] = fillStructLifetimeHist_censored(data)
%
% INPUT     data:   experiment structure, which has to contain the fields
%                   .source
%                   .lftHist_InRegion
%                   .lftHist_OutRegion
%                   source is the path to the data location; at this
%                   location, the function reads the lftInfo
%                   from a folder called LifetimeInfo
%
% OUTPUT    data:   creates new fields
%                   .survivalFunction_InRegion
%                   .survivalFunction_OutRegion
% REMARKS 
%
%
% last modified DATE: 21-May-2008 (Dinah)



% determine survival function for each movie 

for k=1:length(data)
    
    if isfield(data,'lftHist_InRegion')
        currHistIn = data(k).lftHist_InRegion;
        currHistOut = data(k).lftHist_OutRegion;
    else
        error('function requires existence of a structure field .lftHist_InRegion');
    end
          
    % lifetime function
    pdefIn = find(isfinite(currHistIn));
    hfIn = currHistIn(pdefIn);
    cfIn = cumsum(hfIn);
    currLifetimeFunctionIn = currHistIn;
    currLifetimeFunctionIn(pdefIn) = cfIn;
    
    % lifetime function
    pdefOut = find(isfinite(currHistOut));
    hfOut = currHistOut(pdefOut);
    cfOut = cumsum(hfOut);
    currLifetimeFunctionOut = currHistOut;
    currLifetimeFunctionOut(pdefOut) = cfOut;
    
    
    % survival function
    currMaxIn = max(currLifetimeFunctionIn);
    currSurvivalFunctionIn = currMaxIn - currLifetimeFunctionIn;
    
    % survival function
    currMaxOut = max(currLifetimeFunctionOut);
    currSurvivalFunctionOut = currMaxOut - currLifetimeFunctionOut;
    
    data(k).survivalFunction_InRegion = currSurvivalFunctionIn;
    data(k).survivalFunction_OutRegion = currSurvivalFunctionOut;
end



end % of function







    
