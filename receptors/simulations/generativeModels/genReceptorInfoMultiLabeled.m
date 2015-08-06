function receptorInfoLabeled = genReceptorInfoMultiLabeled(receptorInfoAll,...
    labelRatio,intensityQuantum)
%GENRECEPTORINFOLABELED sub-samples simulated receptor trajectories and outputs information for labeled subset, including compound tracks
%
%SYNPOSIS receptorInfoLabeled = genReceptorInfoLabeled(receptorInfoAll,...
%    labelRatio,intensityQuantum)
%
%INPUT  receptorInfoAll : Output of receptorAggregationSimple_new
%       labelRatio      : Matrix of labeling ratios.
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
%
%Modified June 2015, Paul Blazek, to include multi-colored labeling.

%% Input/Output

receptorTraj = receptorInfoAll.receptorTraj;
recept2clustAssign = receptorInfoAll.recept2clustAssign;
clust2receptAssign = receptorInfoAll.clust2receptAssign;

[numReceptors,~,numIterSim] = size(receptorTraj);

%09/05/14 (ryirdaw) - modified to allow a vector labelRatio
%06/04/15 (pblaze) -modified to allow matrix labelRatio (i.e. multi-color)
numLabelRatio = sum(sum(labelRatio>0));
receptorInfoLabeled(numLabelRatio,1) = struct('receptorTraj',[],...
    'recept2clustAssign',[],...
    'clust2receptAssign',[],...
    'compTracks',[],...
    'labelRatio',[],...
    'samplingNum',[]);

%% Labeling and sub-sampling
%06/04/15 (pblaze) the results from each sampling are placed in separate
%rows, so nextColor serves as a counter across multiple colors and
%samplings
nextColor = 1;
for lRindx=1:length(labelRatio)
    
    %06/04/14 (pblaze) - checks if total labeled not greater than 1 for each
    %sampling, otherwise it throws an error
    if sum(labelRatio(lRindx,:)) > 1
        for color=1:sum(labelRatio(lRindx,:)>0)
            receptorInfoLabeled(nextColor).labelRatio = labelRatio(lRindx,color);
            receptorInfoLabeled(nextColor).samplingNum = lRindx;
            nextColor = nextColor + 1;
        end
        fprintf('\n============\nLabeling greater than 1 in sampling %d, so it was skipped.\n============\n',lRindx);
        continue
    end
    
    
    %if labeling ratio is less than one ...
    if labelRatio(lRindx,1) < 1
        
        %number of colors in this sampling
        numColors = sum(labelRatio(lRindx,:)>0);
        
        %label some receptors based on the labeling ratio
        %n>0 = number of color, 0 = unlabeled
        %09/05/14 (ryirdaw) - modified to allow a vector labelRatio
        %06/04/15 (pblaze) - modified to allow matrix labelRatio (i.e.
        %multi-color)
        randLabel = rand(numReceptors,1);
        labelFlag = ones(numReceptors,1);
        for color = 1:numColors
            labelFlag = labelFlag + (randLabel >= sum(labelRatio(lRindx,1:color)));
        end
        labelFlag = mod(labelFlag,numColors+1);
        
        %06/04/14 (pblaze) iterate over all the colors in the sample
        for color = 1:numColors
            indxLabeled = find(labelFlag==color);
            indxNotLabeled = find(labelFlag~=color);
            
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
            fprintf('\n============\nTime for convReceptClust2CompTracks sample %d, color %d is %g seconds. \n============\n',...
                lRindx,color,convTracksETime);

            %put information in receptorInfoLabeled
            %09/05/14 (ryirdaw) - modified to allow a vector labelRatio
            %{
            %original
            receptorInfoLabeled = struct('receptorTraj',receptorTrajLabeled,...
                'recept2clustAssign',recept2clustAssignLabeled,...
                'clust2receptAssign',clust2receptAssignLabeled,...
                'compTracks',compTracksLabeled);
            %}
            receptorInfoLabeled(nextColor).receptorTraj = receptorTrajLabeled;
            receptorInfoLabeled(nextColor).recept2clustAssign = recept2clustAssignLabeled;
            receptorInfoLabeled(nextColor).clust2receptAssign = clust2receptAssignLabeled;
            receptorInfoLabeled(nextColor).compTracks = compTracksLabeled;
            receptorInfoLabeled(nextColor).labelRatio = labelRatio(lRindx,color);
            receptorInfoLabeled(nextColor).samplingNum = lRindx;
            
            nextColor = nextColor + 1;
        end
      
    else %if all receptors are labeled
        
        %09/05/14 (ryirdaw) - modified to allow a vector labelRatio
        %{
        %original        
        %copy receptor information from all to labaled
        receptorInfoLabeled = receptorInfoAll;
        %}
        
        receptorInfoLabeled(nextColor).receptorTraj = receptorInfoAll.receptorTraj;
        receptorInfoLabeled(nextColor).recept2clustAssign = receptorInfoAll.recept2clustAssign;
        receptorInfoLabeled(nextColor).clust2receptAssign = receptorInfoAll.clust2receptAssign;
        receptorInfoLabeled(nextColor).labelRatio = labelRatio(lRindx,1);
        receptorInfoLabeled(nextColor).samplingNum = lRindx;
        
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
        fprintf('\n============\nTime for convReceptClust2CompTracks sample %d, color 1 is %g seconds. \n============\n',...
                lRindx,convTracksETime);

        %09/05/14 (ryirdaw) - modified to allow a vector labelRatio
        %{
        %original        
        %store the new compound tracks with the labeling intensity
        receptorInfoLabeled.compTracks = compTracksLabeled;
        %}
        receptorInfoLabeled(nextColor).compTracks = compTracksLabeled;
        
        nextColor = nextColor + 1;

    end %(if labelRatio < 1 ... else ...)
    
end % for each labelRatio element

end