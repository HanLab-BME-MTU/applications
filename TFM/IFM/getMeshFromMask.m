function [msh,borderE,borderSeg,exBndE,exBndSeg,numEdges,bndInd,ind] = getMeshFromMask(movieData,jj,curFlow, mask,minImgSize,numSubDoms)
        % get the image boundaries
        LeftUpperCorner(1:2) = [min(curFlow(:,2)), min(curFlow(:,1))];
        RightLowerCorner(1:2) = [max(curFlow(:,2)), max(curFlow(:,1))];
        
        mask(RightLowerCorner(2)+1:end,:,jj) = 0; % for bottom
        mask(:,RightLowerCorner(1)+1:end,jj) = 0; % for right
        mask(1:LeftUpperCorner(2)-1,:,jj) = 0; % for top
        mask(:,1:LeftUpperCorner(1)-1,jj) = 0; % for left

        B = bwboundaries(mask);
        % Build a Geometry Description Matrix - I'll use curve (polygon
        % solid)
        % find the straight lines first
        xmin = min(B{1}(:,2));
        xmax = max(B{1}(:,2));
        ymin = min(B{1}(:,1));
        ymax = max(B{1}(:,1));
        iFreeEdge = []; % the id number of free edge (1: left, 2:top, 3:right, 4:bottom, these can be multiple)
        
        % rotating from left, top, right to bottom, check the straightness
        % (if the line is inner boundary vs. free boundary (curved))
        % Previously we assumed that there is only one free boundary, but
        % it was not the case all the time. Thus, now we assume at least
        % one straight boundary should be there. 
        % First find the straight edge from four edges
        imgBoundaryRatio = 0.3; % ratio of straight edge compared to total height or width of image.
        % for left edge.
        xminIdx = B{1}(:,2)==xmin; % logical index
        if sum(xminIdx) > max(minImgSize, imgBoundaryRatio*(ymax-ymin)) && ...
            (sum(diff(sort(B{1}(xminIdx,1)))==1) == sum(xminIdx)-1  || ...
            sum(diff(sort(B{1}(xminIdx,1)))==1) == sum(xminIdx)-2) % see if they are consecutive
            curveL = [xmin max(B{1}(xminIdx,1)) xmin min(B{1}(xminIdx,1))]; % from bottom to top
            [~,curveLIdx1] = max(B{1}(xminIdx,1));
            [~,curveLIdx2] = min(B{1}(xminIdx,1));
            curveLIdxIdx = find(xminIdx);
            curveLIdx = curveLIdxIdx([curveLIdx1;curveLIdx2]);
        else
            curveL = [];
            xminIdx = false(size(xminIdx));
            iFreeEdge = [iFreeEdge 1];
        end
        % for top edge
        yminIdx = B{1}(:,1)==ymin; % logical index
        if sum(yminIdx) > max(minImgSize, imgBoundaryRatio*(xmax-xmin)) && ...
            (sum(diff(sort(B{1}(yminIdx,2)))==1) == sum(yminIdx)-1  || ...
            sum(diff(sort(B{1}(yminIdx,2)))==1) == sum(yminIdx)-2) % see if they are consecutive
            curveT = [min(B{1}(yminIdx,2)) ymin max(B{1}(yminIdx,2)) ymin]; % from left to right
            [~,curveTIdx1] = min(B{1}(yminIdx,2));
            [~,curveTIdx2] = max(B{1}(yminIdx,2));
            curveTIdxIdx = find(yminIdx);
            curveTIdx = curveTIdxIdx([curveTIdx1;curveTIdx2]);
        else
            % curve approximation with ...
            curveT = [];
            yminIdx = false(size(yminIdx));
            iFreeEdge = [iFreeEdge 2];
        end
        % for right edge
        xmaxIdx = B{1}(:,2)==xmax; % logical index
        if sum(xmaxIdx) > max(minImgSize, imgBoundaryRatio*(ymax-ymin)) && ...
            (sum(diff(sort(B{1}(xmaxIdx,1)))==1) == sum(xmaxIdx)-1  || ...
            sum(diff(sort(B{1}(xmaxIdx,1)))==1) == sum(xmaxIdx)-2) % see if they are consecutive
            curveR = [xmax min(B{1}(xmaxIdx,1)) xmax max(B{1}(xmaxIdx,1))]; % from top to bottom
            [~,curveRIdx1] = min(B{1}(xmaxIdx,1));
            [~,curveRIdx2] = max(B{1}(xmaxIdx,1));
            curveRIdxIdx = find(xmaxIdx);
            curveRIdx = curveRIdxIdx([curveRIdx1;curveRIdx2]);
        else
            curveR = [];
            xmaxIdx = false(size(xmaxIdx));
            iFreeEdge = [iFreeEdge 3];
        end
        % for bottom edge
        ymaxIdx = B{1}(:,1)==ymax; % logical index
        if sum(ymaxIdx) > max(minImgSize, imgBoundaryRatio*(xmax-xmin)) && ...
            (sum(diff(sort(B{1}(ymaxIdx,2)))==1) == sum(ymaxIdx)-1  || ...
            sum(diff(sort(B{1}(ymaxIdx,2)))==1) == sum(ymaxIdx)-2) % see if they are consecutive
            curveB = [max(B{1}(ymaxIdx,2)) ymax min(B{1}(ymaxIdx,2)) ymax]; % from bottom to top
            [~,curveBIdx1] = max(B{1}(ymaxIdx,2));
            [~,curveBIdx2] = min(B{1}(ymaxIdx,2));
            curveBIdxIdx = find(ymaxIdx);
            curveBIdx = curveBIdxIdx([curveBIdx1;curveBIdx2]);
        else
            curveB = [];
            ymaxIdx = false(size(ymaxIdx));
            iFreeEdge = [iFreeEdge 4];
        end
        
        % Once the strainght edge is identified, divide the free edges to
        % the rest of (4 - nFreeEdge) edges.
        numFreeEdges = length(iFreeEdge);
        iStrEdge = setdiff([1,2,3,4],iFreeEdge);
        if length(iStrEdge)==1 || ...
                (length(iStrEdge)>1 && sum(diff(iStrEdge)==1) == length(iStrEdge)-1) %if straight edge is only one or consecutively arranged,
            % get free curve indices
            if ismember(
        else
        end
%         % curve approximation with ...
%         freeIdx = ~(xminIdx | yminIdx | xmaxIdx | ymaxIdx);% index for free edge
%         freeIdx(max(find(freeIdx,1)-1,1)) = true;
%         freeIdx(min(find(freeIdx,1,'last')+1,end)) = true;
%         
%         % determine distance between nodes based on speckle density
%         numSpeckles=length(flow);
%         areaImg=prod(RightLowerCorner-LeftUpperCorner);
%         interSpecDist=ceil(sqrt(areaImg/numSpeckles));  % the number of skipping points for equi-spatial sampling of curves
% 
%         freeIdxIdx = find(freeIdx,1):interSpecDist:find(freeIdx,1,'last');
%         if freeIdxIdx(end) ~= find(freeIdx,1,'last')
%             freeIdxIdx(end) = find(freeIdx,1,'last');
%         end
%         
% %         freeIdxS = false(size(freeIdx));
% %         freeIdxS(freeIdxIdx) = true;
%         
%         switch nFreeEdge
%             case 1
%                 curveLIdx = freeIdxIdx;
%              case 2
%                 curveTIdx = freeIdxIdx;
%             case 3
%                 curveRIdx = freeIdxIdx;
%             case 4
%                 curveBIdx = freeIdxIdx;
%             case 0
%                 disp('Something is wrong. There is no free edge. Check your boundary condition.')
%         end
        
%         allEdgeIdx = xminIdx | yminIdx | xmaxIdx | ymaxIdx | freeIdxS;
        np = length(curveLIdx)+length(curveTIdx)+length(curveRIdx)+length(curveBIdx)-4; % the number of polygon segments
        gd = zeros(2+2*np,1); % initialization of geometry description matrix
        gd(1) = 2; % represents polygon solid
        gd(2) = np;
        gd(3:2+np) = B{1}([curveLIdx(1:end-1); curveTIdx(1:end-1)'; curveRIdx(1:end-1); curveBIdx(1:end-1)],2); % x-coordinates
        gd(3+np:2+2*np) = B{1}([curveLIdx(1:end-1); curveTIdx(1:end-1)'; curveRIdx(1:end-1); curveBIdx(1:end-1)],1); % y-coordinates
        dl = decsg(gd); % decompose constructive solid geometry into minimal regions
        %% mesh
        [msh.p,msh.e,msh.t]=initmesh(dl,'hmax',2*interSpecDist); 
        iActinChannel = 3;
        curActin = movieData.getChannel(iActinChannel).loadImage(jj);
        figure, imshow(curActin,[])
        hold on
        plot(curFlow(:,2),curFlow(:,1),'y.')
        quiver(curFlow(:,2),curFlow(:,1),3*(curFlow(:,4)-curFlow(:,2)),3*(curFlow(:,3)-curFlow(:,1)),0,'y')
        pdemesh(msh.p,msh.e,msh.t)
        %Separate the boundary segments into internal border and real boundaryies.
        borderE   = find(msh.e(6,:)~=0&msh.e(7,:)~=0);
        borderSeg = unique(msh.e(5,borderE));
        exBndE    = find(msh.e(6,:)==0|msh.e(7,:)==0);
        exBndSeg  = unique(msh.e(5,exBndE));
        numEdges  = length(exBndSeg);
        %Group subdomains and boundaries.
%         numSubDoms = 1;
        ind    = ones(1,numSubDoms);
        bndInd = {};
        for k = 1:numEdges
         bndInd{k} = exBndSeg(k);
        end

