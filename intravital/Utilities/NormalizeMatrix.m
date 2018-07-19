function B = NormalizeMatrix( A )
% Normalizes the input matrix to have values between [0,1]

    B = A(:) - min(A(:)) / (max(A(:)) - min(A(:)));

end