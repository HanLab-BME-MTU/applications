function plotClusters(recept2clust)
%PLOTCLUSTERS is a function to easily visualize clusters in a simulation.  
%   recept2clust is used to determine the clustering of each receptor at 
%   each iteration. Plot shows receptors along the Y axis with the 
%   iterations along the X. Clusters are indicated by parallel points for 
%   the respective receptors. Receptors that belong to the same cluster are
%   shown with the same color. For example, an iteration point with a 
%   plotted point for each receptor indicates all receptors are clustered
%   together (1 cluster and a single color on plot).  Blank areas for 
%   receptors indicate no clustering.
%   
%   INPUT:
%           recept2clust:   matrix of receptor-to-cluster relationships
%                           spanning all receptors and iterations or just
%                           smaller portions for easier visualization
%
%   OUPTUT:
%       none.
%
%   NOTE: 
%       1) Repetitive call to plot makes function slow.
%       2) Six different colors only - if more than six clusters at an
%       iteration point, a magenta diamond is used.
%
%   Robel Yirdaw, 07/05/13
%

    [numReceptors,iters] = size(recept2clust);
    %Use the following index array to avoid using find within loop
    recepIndx = (1:numReceptors)';
    %Plot color/point scheme set
    cpSet = {'ko';'r+';'c*';'gs';'b<';'yx'};

    figH = figure();    
    hold on
    for iterIndx=1:iters
        %Reset plot color/point counter
        cpSetCounter = 1;
        [unqR2C,~,~] = unique(recept2clust(:,iterIndx),'stable');
        if ( numel(unqR2C) < numReceptors)
            %Found clusters
            for clustIndx=1:numel(unqR2C)
                %Determine if cluster is shared - i.e. receptors in cluster
                %recepInClustIndx = find(recept2clust(:,iterIndx)==unqR2C(clustIndx));  
                recepInClust = recept2clust(:,iterIndx)==unqR2C(clustIndx);  
                %Count number of receptors in current receptor's cluster
                %numRecepInClust = numel(recepInClustIndx);                
                numRecepInClust = sum(recepInClust);                
                if (numRecepInClust > 1)
                    %plot(repmat(iterIndx,numRecepInClust,1),recepInClustIndx,'ro');
                    if (cpSetCounter > length(cpSet))
                        plot(repmat(iterIndx,numRecepInClust,1),recepIndx(recepInClust),'md');
                    else
                        plot(repmat(iterIndx,numRecepInClust,1),recepIndx(recepInClust),cpSet{cpSetCounter});
                    end
                    cpSetCounter = cpSetCounter + 1;
                end
            end
            
        end %if clusters
        
    end %for each iteration
    
    hold off
    %Get axes handle and make adjustments
    axesH = get(figH,'CurrentAxes');
    set(axesH,'XGrid','on','YGrid','on');
    set(get(axesH,'XLabel'),'String','Iteration');
    set(get(axesH,'YLabel'),'String','Receptor');
    set(axesH,'YTick',1:numReceptors);
    ylim(axesH,[0.8 numReceptors+0.2]);
    
end %function
