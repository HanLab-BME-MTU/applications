function candsTot=fsmPrepCheckInfo(candsTot)

% fsmPrepCheckInfo checks for repetitions amongst the lists of the primary, 
% secondary and tertiary speckles and removes them from the cands structure
%
%
% SYNOPSIS   candsTot=fsmPrepCheckInfo(candsTot)
%
% INPUT      candsTot    :   cands structure
%
% OUTPUT     candsTot    :   updated cands
%
%
% DEPENDENCES   fsmPrepCheckInfo uses { createSparseDistanceMatrix }
%               fsmPrepCheckInfo is used by { fsmPrepMainSecondarySpeckles }
%
% Alexandre Matov, November 7th, 2002


c=0; % counter
rep=[]; % repetitions

Lmax=[candsTot.Lmax]; % coordinates
candsTotPos=reshape(Lmax,2,length(Lmax)/2)'; % insignificant loc max positions

nula=1e-10; % a very small number

D=createSparseDistanceMatrix(candsTotPos,candsTotPos,.9,nula); % finds the repetitions

[y x]=find(D==nula);

if length(y)>length(candsTot)
    for i=1:length(y)
        if y(i)~=x(i)
            c=c+1;
            if candsTot(y(i)).status==0
                rep(c)=y(i);
            elseif candsTot(x(i)).status==0
                rep(c)=x(i);
            else
                error('unusual speckles repetition: both speckles (at the same position) are significant');
            end
        end
    end
end

candsTot(rep)=[];


    