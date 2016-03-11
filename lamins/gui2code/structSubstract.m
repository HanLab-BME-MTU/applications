function [ C ] = structSubstract( A, B )
%structSubstract struct substracts struct B from A to obtain struct C
% C contains all the fields in struct A that are not in B and do not have
% the same value in B
%
% Mark Kittisopikul
% Jaqaman Lab
% UT Southwestern
% November 19th, 2014

C = A;

if(isscalar(A))
    
    if(~isscalar(B))
        % quit is A is scalar and B is not
        return;
    end

    % Iterate through A's fields
    fields = fieldnames(A);
    for f = 1:length(fields)
        fieldName = fields{f};
        A_value = A.(fieldName);
        % If B has the field...
        if(isfield(B,fieldName))
            B_value = B.(fieldName);
            if( isequal( A_value , B_value) )
                % and the value in A and B are the same, remove it from A
                C = rmfield(C,fieldName);
            elseif(isstruct(A_value) && isstruct(B_value))
                % otherwise if the value is a struct, recurse
                C.(fields{f}) = structSubstract(A_value,B_value);
            end
        end
    end
    
elseif(length(A) == length(B))
    % if the length of A and B are the same compare the elements pairwise
    C = arrayfun(@structSubstract,A,B);
else
    % otherwise if A and B are of different size, quit, return A
    return;    
end

end

