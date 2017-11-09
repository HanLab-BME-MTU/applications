function [ shift, C ] = registerInt( A, B, maxShift )
%registerInt Determine shift by correlation

if(nargin < 3)
    maxShift = 3;
end

A = double(A) - mean(A(:));
sizeA = size(A);
B = double(B) - mean(B(:));
sizeB = size(B);

shifts = -maxShift:maxShift;
nShifts = length(shifts);

C = zeros(nShifts,nShifts,nShifts);

for i = 1:nShifts
    for j = 1:nShifts
        for k = 1:nShifts
            aS = max(1 + shifts([i j k]),1);
            aE = min(sizeA + shifts([i j k]),sizeA);
            bS = max(1 - shifts([i j k]),1);
            bE = min(sizeB - shifts([i j k]),sizeB);
            C(i,j,k) = corr2(joinColumns([],A(aS(1):aE(1),aS(2):aE(2),aS(3):aE(3))), ...
                             joinColumns([],B(bS(1):bE(1),bS(2):bE(2),bS(3):bE(3))));
        end
    end
end

[~,maxidx] = max(C(:));
[shift(1),shift(2),shift(3)] = ind2sub(size(C),maxidx);
shift = shifts(shift);


end

