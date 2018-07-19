function [data] = fillStructSurvivalFunction(data, fieldname)
% fillStructSurvivalFunction fills in the values for the survival function
% based on the lifetime dat (lftHist_censored)
%
% SYNOPSIS [data] = fillStructLifetimeHist_censored(data)
%
% INPUT     data:   experiment structure, which has to contain the fields
%                   .source
%                   .lftHist_censored  
%                   source is the path to the data location; at this
%                   location, the function reads the lftInfo
%                   from a folder called LifetimeInfo
%
% OUTPUT    data:   creates a new field
%                   .lftHist_censored  
% REMARKS 
%
%
% last modified DATE: 21-May-2008 (Dinah)
% Francois Aguet 01/21/2010


% determine survival function for each movie 
if (nargin == 1)
    fieldname = 'lftHist_censored';
end

for k = 1:length(data)
    
    if isfield(data, fieldname)
        currHist = data(k).(fieldname);
    else
        error(['Function requires the field ' fieldname ' in ' data(k) '.']);
    end
          
    % lifetime function
    pdef = find(isfinite(currHist));
    hf = currHist(pdef);
    cf = cumsum(hf);
    currLifetimeFunction = currHist;
    currLifetimeFunction(pdef) = cf;
    
    % survival function    
    data(k).survivalFunction = max(currLifetimeFunction) - currLifetimeFunction;
end