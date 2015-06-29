function receptorInfoLabeled = genReceptorInfoLabeled(receptorInfoAll,...
    labelRatio,intensityQuantum)
%GENRECEPTORINFOLABELED sub-samples simulated receptor trajectories and outputs information for labeled subset, including compound tracks
%
%SYNPOSIS receptorInfoLabeled = genReceptorInfoLabeled(receptorInfoAll,...
%    labelRatio,intensityQuantum)
%
%INPUT  receptorInfoAll : Output of receptorAggregationSimple_new
%       labelRatio      : Column vector of labeling ratios.
%                         Optional. Default: 1 (i.e. all).
%       intensityQuantum: Row vector with 2 values, the mean and std of
%                         receptor labeling intensities.
%                         Optional. Default: [1 0].
%
%OUTPUT receptorInfoLabeled: Structure array with number of elements =
%                         number of different labeling ratios. For each
%                         labeling ratio, fields are similar to
%                         receptorInfoAll, but for the labaled
%                         receptors only. It has two additional fields:
%               .compTracks : The receptor trajectories and interactions over
%                             time, as converted by convReceptClust2CompTracks.
%                             Same format as output of trackCloseGapsKalman.
%               .labelRatio : The labeling ratio used to subsample.
%
%Code started by Robel Yirdaw in 2014, initially as a copy-paste from
%receptorAggregationSimple_new.
%
%Modified May 2015, Khuloud Jaqaman

%% Input/Output

receptorTraj = receptorInfoAll.receptorTraj;
recept2clustAssign = receptorInfoAll.recept2clustAssign;
clust2receptAssign = receptorInfoAll.clust2receptAssign;

[numReceptors,~,numIterSim] = size(receptorTraj);

%09/05/14 (ryirdaw) - modified to allow a vector labelRatio
numLabelRatio = length(labelRatio);
receptorInfoLabeled(numLabelRatio,1) = struct('receptorTraj',[],...
    'recept2clustAssign',[],...
    'clust2receptAssign',[],...
    'compTracks',[],...
    'labelRatio',[]);

%% Labeling and sub-sampling

for lRindx=1:numLabelRatio
    
    %if labeling ratio is less than one ...
    if labelRatio(lRindx) < 1
        
        %label some receptors based on the labeling ratio
        %1 = labeled, 0 = unlabeled
        %09/05/14 (ryirdaw) - modified to allow a vector labelRatio
        labelFlag = rand(numReceptors,1) <= labelRatio(lRindx);
        indxLabeled = find(labelFlag==1);
        indxNotLabeled = find(labelFlag==0);
        
        %extract the trajectories of the labeled receptors
        receptorTrajLabeled = receptorTraj(indxLabeled,:,:);
        
        %assign the intensities of the labeled receptors
        receptorIntensityLabeled = intensityQuantum(1) + ...
            randn(length(indxLabeled),numIterSim) * intensityQuantum(2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %12/06/13 (ryirdaw)
        %Reset intensity values that are < epsilon to epsilon
        receptorIntensityLabeled(receptorIntensityLabeled < eps) = eps;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %extract the cluster assignments of the labeled receptors
        recept2clustAssignLabeled = recept2clustAssign(indxLabeled,:);
        
        %modify the cluster-to-receptor assignments to include only labeled
        %receptors
        clust2receptAssignLabeled = clust2receptAssign;
        %convert receptors that are not labeled to zero
        for iNotLabeled = indxNotLabeled'
            clust2receptAssignLabeled(clust2receptAssignLabeled==iNotLabeled) = 0;
        end
        %update the indices of the labeled receptors
        for iLabeled = 1 : length(indxLabeled)
            clust2receptAssignLabeled(clust2receptAssignLabeled==indxLabeled(iLabeled)) = iLabeled;
        end
        %make sure that zeros come after receptor indices
        clust2receptAssignLabeled = sort(clust2receptAssignLabeled,2,'descend');
        
        %remove empty clusters (which include unlabeled receptors)
        %modify receptor-to-cluster assignments accordingly
        for iIter = 1 : numIterSim
            clustSize = sum(clust2receptAssignLabeled(:,:,iIter)~=0,2);
            indxFull = find(clustSize~=0);
            indxEmpty = find(clustSize==0);
            clust2receptAssignLabeled(:,:,iIter) = clust2receptAssignLabeled(...
                [indxFull;indxEmpty],:,iIter);
            for iFull = 1 : length(indxFull)
                recept2clustAssignLabeled(recept2clustAssignLabeled(:,iIter)...
                    ==indxFull(iFull),iIter) = iFull;
            end
        end
        
        %remove empty rows and columns from clust2receptAssign
        cluster2receptor = max(clust2receptAssignLabeled,[],3);
        columnSum = sum(cluster2receptor);
        clust2receptAssignLabeled = clust2receptAssignLabeled(:,columnSum~=0,:);
        rowSum = sum(cluster2receptor,2);
        clust2receptAssignLabeled = clust2receptAssignLabeled(rowSum~=0,:,:);
        
        convTracksStartTime = tic;
        
        %put labeled receptor trajectories and clusters into the format of the
        %output of trackCloseGapsKalman
        compTracksLabeled = convReceptClust2CompTracks(clust2receptAssignLabeled,...
            recept2clustAssignLabeled,receptorTrajLabeled,receptorIntensityLabeled);
        
        convTracksETime = toc(convTracksStartTime);
        fprintf('\n============\nTime for convReceptClust2CompTracks is %g seconds. \n============\n',convTracksETime);
        
        %put information in receptorInfoLabeled
        %09/05/14 (ryirdaw) - modified to allow a vector labelRatio
        %{
        %original
        receptorInfoLabeled = struct('receptorTraj',receptorTrajLabeled,...
            'recept2clustAssign',recept2clustAssignLabeled,...
            'clust2receptAssign',clust2receptAssignLabeled,...
            'compTracks',compTracksLabeled);
        %}
        receptorInfoLabeled(lRindx).receptorTraj = receptorTrajLabeled;
        receptorInfoLabeled(lRindx).recept2clustAssign = recept2clustAssignLabeled;
        receptorInfoLabeled(lRindx).clust2receptAssign = clust2receptAssignLabeled;
        receptorInfoLabeled(lRindx).compTracks = compTracksLabeled;
        receptorInfoLabeled(lRindx).labelRatio = labelRatio(lRindx);
      
    else %if all receptors are labeled
        
        %09/05/14 (ryirdaw) - modified to allow a vector labelRatio
        %{
        %original        
        %copy receptor information from all to labaled
        receptorInfoLabeled = receptorInfoAll;
        %}
        
        receptorInfoLabeled(lRindx).receptorTraj = receptorInfoAll.receptorTraj;
        receptorInfoLabeled(lRindx).recept2clustAssign = receptorInfoAll.recept2clustAssign;
        receptorInfoLabeled(lRindx).clust2receptAssign = receptorInfoAll.clust2receptAssign;
        receptorInfoLabeled(lRindx).labelRatio = labelRatio(lRindx);
        
        %assign receptor intensities
        receptorIntensity = intensityQuantum(1) + ...
            randn(numReceptors,numIterSim) * intensityQuantum(2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %12/06/13 (ryirdaw)
        %Reset intensity values that are < epsilon to epsilon
        receptorIntensity(receptorIntensity < eps) = eps;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        convTracksStartTime = tic;

        %put labeled receptor trajectories and clusters into the format of the
        %output of trackCloseGapsKalman
        compTracksLabeled = convReceptClust2CompTracks(clust2receptAssign,...
            recept2clustAssign,receptorTraj,receptorIntensity);

        convTracksETime = toc(convTracksStartTime);
        fprintf('\n============\nTime for convReceptClust2CompTracks is %g seconds. \n============\n',convTracksETime);

        %09/05/14 (ryirdaw) - modified to allow a vector labelRatio
        %{
        %original        
        %store the new compound tracks with the labeling intensity
        receptorInfoLabeled.compTracks = compTracksLabeled;
        %}
        receptorInfoLabeled(lRindx).compTracks = compTracksLabeled;

    end %(if labelRatio < 1 ... else ...)
    
end % for each labelRatio element

end