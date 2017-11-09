function [ C ] = structMerge( A , B)
% structMerge merges struct A and B together with values in A taking
% precedence
%
% The output C is a struct that contains all the fields in A and B with the
% values of A if the fields exist in both A and B
%
% See also mergestruct
%
% Mark Kittisopikul
% Jaqaman Lab
% UT Southwestern
% November 19th, 2014

C = A;

if(isscalar(A))
    
    if(~isscalar(B))
        return;
    end

    % Iterature through B's fields
    fields = fieldnames(B);
    for f = 1:length(fields)
        fieldName = fields{f};
        B_value = B.(fieldName);
        % If A has that field
        if(isfield(A,fieldName))
            A_value = A.(fieldName);
            % and the value is a struct, then recurse
            if(isstruct(A_value) && isstruct(B_value))
                C.(fields{f}) = structMerge(A_value,B_value);
            end
            % if value is not a struct, do nothing
        else
            % if A does not have that field, then add it
            C.(fieldName) = B_value;
        end
    end
    
elseif(length(A) == length(B))
    % iterate over equally sized arrays
    C = arrayfun(@structMerge,A,B);
else
    % otherwise quit
    return;    
end

end