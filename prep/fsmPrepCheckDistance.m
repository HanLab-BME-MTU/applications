function cands2=fsmPrepCheckDistance(cands2,cands1)

% fsmPrepCheckDistance uses a look-up-table to verify if a new (secondary or tertiary) speckle are not closer than 
% theoretically possible to a primary (for tertiary also secondary) one; if a candidate is closer - it is discarded
%
% SYNOPSIS    cands2=fsmPrepCheckDistance(cands2,cands1)
%
% INPUT      cands2      :   cands for the new/secondary speckles
%            cands1      :   cands for the old/primary speckles
%
% OUTPUT     cands2      :   updated cands for the new/secondary speckles
%
%
% DEPENDENCES   fsmPrepCheckDistance uses { createDistanceMatrix }
%               fsmPrepCheckDistance is used by { fsmPrepMainSecondarySpeckles }
%
% Alexandre Matov, November 7th, 2002
% Modified by Sylvain Berlemont, 2010

validIdx1 = [cands1(:).status] == 1;
validIdx2 = [cands2(:).status] == 1;

R1 = vertcat(cands1(validIdx1).Lmax);
R2 = vertcat(cands2(validIdx2).Lmax);

deltaI1 = [cands1(validIdx1).deltaI];
deltaI2 = [cands2(validIdx2).deltaI];

D=createSparseDistanceMatrix(R2,R1,5.21);

toBeRemoved=false(size(D,1),1); % initialization

for i=1:size(D,1)
    rowI = D(i,:);
    minDistRow=find(rowI);
    
    if ~isempty(minDistRow)
        % SB: this criteria doesn't make sense to me:
        RelInt = deltaI1(minDistRow) / deltaI2(i);
        r=5.2-3.5.*exp(.5-RelInt(:));
        d = nonzeros(rowI) + eps;
            
        toBeRemoved(i) = any(d <= r);
    end
end

if any(toBeRemoved)
    ind2 = find(validIdx2);
    cands2(ind2(toBeRemoved)) = [];
end

 