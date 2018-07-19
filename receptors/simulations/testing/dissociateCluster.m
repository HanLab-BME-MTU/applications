function dissociateFlag = dissociateCluster(dissociationProb,clusterMembers)
%DISSOCIATECLUSTER determines receptors to dissociate based on the total 
%number of receptors.  
%   Allows for multiple receptors to dissociate instantly. The
%   probability of removing r receptors from a cluster of size N is
%   calculated as p(1)^r, if r = N-1 or (N choose r)*p(1)^r, if r < N-1.
%   See 072213_summary for details.
%
%   INPUT:      
%           dissociationProb - dissociation probability
%           clusterMembers   - a vector of receptors forming the cluster
%
%   OUTPUT:
%           dissociateFlag   - a vector of the same size as clusterMembers
%                              with boolean values indicating receptors
%                              that will or will not dissociate (1/0)
%
%   Robel Yirdaw, 07/22/13
%

    %Current dissociating cluster size
    clusterSize = length(clusterMembers);
    %Initialize dissociation flag
    dissociateFlag = zeros(clusterSize,1);
    %The calculated probabilities
    allDissProb = zeros(clusterSize-1,1);
    %The maximum number of receptors that can be removed
    maxRecep2Remove = clusterSize - 1;
        
    %Calculate probabilities as described above.
    for recep2remove=1:maxRecep2Remove
        term1 = power(dissociationProb,recep2remove);        
        if (recep2remove < (maxRecep2Remove))
            term2 = factorial(clusterSize)/(factorial(recep2remove)*factorial(clusterSize - recep2remove));
        else
            term2 = 1;
        end
        allDissProb(recep2remove,1) = (term1*term2);        
    end
  
    %Determine which "event" will occur
    dissEvent = rand(1) < allDissProb;
    
    %One to maximum number of receptors can dissociate, starting from the
    %first.
    if (any(dissEvent))
        dissociateFlag(1:sum(dissEvent),1) = 1;
    end
    
end %function
    