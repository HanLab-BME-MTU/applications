function [dataPool] = poolExperiments(varargin);
% this function pools different experiment structures into one, e.g.
% different control datasets of different origins
% INPUT:    various input experiments, e.g. data1, data2, data3
%           where each experiment has multiple fields (e.g. .source,
%           .framerate, .movieLength)
% OUTPUT:   merged experiment file
%
% last modified: Dinah Loerke, March 15, 2009


for i=1:length(varargin)
    
    data_curr = varargin{i};
    
    if i==1
        names = fieldnames(data_curr);
        data_all = data_curr;
    else

        clen = length(data_all);
        
        for n=1:length(data_curr)
            curr_struct = data_curr(n);
            for k=1:length(names)
                curr_field  = getfield(curr_struct, names{k}); 
                data_all    = setfield(data_all, {1,clen+n}, names{k}, curr_field);
            end
        end
    end % of if
end % of for-loop

dataPool = data_all;

end % of function
            
            
    
    