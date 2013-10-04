function receptor2cluster = fixClusterLabels(receptor2clusterNew,receptor2clusterPrev)
%FIXCLUSTERLABELS   Maintain the cluster labels of receptors between
%   iterations.
%   The call to the MATLAB function *cluster* changes the cluster 
%   labels for receptors even if there has not been any change in the 
%   makeup of the clusters from previous iteration. If there have not been 
%   any dissociation or aggregation events in the current iteration, keep 
%   the cluster labels for each receptors the same as previous iteration.
%   This results in maintaining the same cluster labels that were assigned
%   when the event occured, for the duration of that configuration.
%   Modification on 09/10/13 
%       1) Removed the receptor-by-receptor check via the
%          while loop since it is no longer necessary since dissociation 
%          and association events can not occur within the same iteration.
%       2) Using call to unique to sort the cluster labels for the
%          receptors in ascending order.
%
%   INPUT:  receptor2clusterNew:    The current receptor2cluster array
%                                   indicating the cluster assignment for 
%                                   each receptor
%           receptor2clusterPrev:   receptor2cluster from previous
%                                   iteration 
%
%   OUTPUT: receptor2cluster:       The modified (or not) receptor2cluster
%
%   Robel Yirdaw, 07/05/13
%       Modified, 09/10/13
%
        
    
    %If there is an aggregation, the max cluster label will be smaller than
    %what it was previously - not possible to have an aggregation event and
    %keep the same cluster labels for the receptors, as the simulation is
    %implemented right now.  Dissociation and association events within the
    %same iteration are also not allowed - this can keep the max cluster
    %label the same while cluster make-ups are different, i.e. 
    %r2cNew = [1;2;2;3;] r2cPrev = [1;2;3;3], where the last cluster
    %dissociated and the free receptor associates with recpetor #2, all
    %within the same iteration.
    
    if ((~isempty(receptor2clusterPrev)) && ...
            max(receptor2clusterNew) == max(receptor2clusterPrev) &&... 
            any(receptor2clusterNew - receptor2clusterPrev) )
        %Same number of clusters but assignment has changed. Revert.
        receptor2cluster = receptor2clusterPrev;
    elseif (isempty(receptor2clusterPrev) || ...
            ((~isempty(receptor2clusterPrev)) && ...
            max(receptor2clusterNew) ~= max(receptor2clusterPrev)) )
        %The first condition covers the very first state.  Otherwise,
        %there has been a change.  Order cluster labels - third output from
        %call to unique with 'stable' option returns labels in ascending
        %order.
        [~,~,receptor2cluster] = unique(receptor2clusterNew,'stable');
    else
        receptor2cluster = receptor2clusterNew;
    end
    
        
        

