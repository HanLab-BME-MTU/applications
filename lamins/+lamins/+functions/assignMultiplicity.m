function matrix = assignMultiplicity(matrix,mult)
% Replace values in a matrix by their multiplicity
%

if(nargin < 2)
    [mult] = getMultiplicityInt(matrix);
end

[uniq,ia,ic] = unique(matrix);
matrix = reshape(mult(ic),size(matrix));

end
