function [ idxCMap ] = gcaPlotLinksByCost( img,labelCandidates,labelMatSeedFilo,candFiloEPs,seedPtsx,seedPtsy,iSeg,weights,cMapLength,varargin )
%gcaPlotLinksByCost... 
% 

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('img'); 
ip.addRequired('labelCandidates'); 
ip.addRequired('labelMatSeedFilo'); 
ip.addRequired('candFiloEPs'); 
ip.addRequired('seedPtsx'); 
ip.addRequired('seedPtsy'); 
ip.addRequired('iSeg'); 
ip.addRequired('weights'); 
ip.addRequired('cMapLength'); 

ip.addParameter('cMapType','parula');


ip.parse(img,labelCandidates,labelMatSeedFilo,candFiloEPs,seedPtsx,seedPtsy,iSeg,weights,cMapLength,varargin{:});
cMapType =str2func(ip.Results.cMapType); 
% Start Plot
        imshow(-img,[]); 
        hold on 
        spy(labelMatSeedFilo>0,'k'); 
        spy(labelCandidates>0,'k',5); 
        scatter(seedPtsx(:),seedPtsy(:),'k','filled'); 
        allCandEPs = vertcat(candFiloEPs{:}); 
        scatter(allCandEPs(:,1),allCandEPs(:,2),'k','filled'); 
        
        % create distance mapper
        cMap=cMapType(cMapLength);
        mapper=linspace(min(weights),max(weights),cMapLength)';
        
        DMap=createDistanceMatrix(weights,mapper);
        [sD,idxCMap]=sort(abs(DMap),2);
        
        for k = 1:length(cMap);
            if sum(idxCMap(:,1)==k)~=0
                toPlot = iSeg(idxCMap(:,1) == k);
                
                cellfun(@(x) plot([x(1,1),x(end,1)],[x(1,2),x(end,2)],'color',cMap(k,:)),toPlot);
                clear toPlot
            else
            end
            % make colormap of costs.
            
            % show each segment in iSeg cell plotted by the costTotal color
            %
        end 
    

end

