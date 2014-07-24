function [path, cost] = findMinimalPathOnCostMatrix(costMatrix, srcNode, destNode, varargin)
%
% Example:
% 
%     costMatrix = zeros(100,100);
%     
%     costMatrix( 27:73, 27:73 ) = 100;  
%     srcNode = sub2ind(size(costMatrix), 25, 25);
%     destNode = sub2ind(size(costMatrix), 75, 75);
% 
%     [path, cost] = findMinimalPathOnCostMatrix(costMatrix, srcNode, destNode);
%     
%     pathNodeSubInd = cell(1, ndims(costMatrix));
%     [pathNodeSubInd{:}] = ind2sub( size(costMatrix), path )
%     
%     figure, imshow( costMatrix );
%     hold on;
%     plot( pathNodeSubInd{2}, pathNodeSubInd{1}, 'r-', 'LineWidth', 2.0 );
%     hold off;
%     
    
    p = inputParser;
    p.addRequired('costMatrix', @(x) (isnumeric(x) && ismember(ndims(x), [2,3])) );
    p.addRequired('srcNode', @(x) (isnumeric(x) && isscalar(x) && x >= 1 && x <= numel(costMatrix)) );
    p.addRequired('destNode', @(x) (isnumeric(x) && isscalar(x) && x >= 1 && x <= numel(costMatrix)) );
    p.addParamValue('costWeight', 0.5, @(x) (isnumeric(x) && isscalar(x) && x >= 0 && x <= 1) );
    p.addParamValue('roiMask', [], @(x) (ismember(ndims(x), [2,3])) );
    p.addParamValue('spacing', ones(1, ndims(costMatrix)), @(x) (isnumeric(x) && numel(x) == ndims(x)) );
    p.addParamValue('debugMode', false, @(x) (islogical(x) && isscalar(x)) );
    p.parse(costMatrix, srcNode, destNode, varargin{:} );
    
    spacing = p.Results.spacing;    
    flagDebugMode = p.Results.debugMode;
    roiMask = p.Results.roiMask;
    costWeight = p.Results.costWeight;
    
    imsize = size(costMatrix);
    imdims = ndims(costMatrix);

    % prepare offsets for pixel neighbors
    neighOffSetSubInd = cell(1, imdims);
    [neighOffSetSubInd{:}] = ind2sub( 3 * ones(1,imdims), (1:3^imdims)' );
    neighOffSetSubInd = cell2mat( neighOffSetSubInd ) - 2;
    neighOffSetSubInd((3^imdims + 1)/2, :) = [];
    
    % intialize distance and predecessor matrices
    dist = zeros( size(costMatrix) ) + Inf;
    dist(srcNode) = 0;
    
    parentIndex = zeros( size(costMatrix) );
    parentIndex( srcNode ) = srcNode;
    
    % find shortest path to destNode using fast-marching
    pq = IndexedPriorityQueue( numel(costMatrix), false );    
    pq.insert(srcNode, dist(srcNode));
    
    if flagDebugMode
        visitTime = zeros(size(costMatrix));
        visitTimer = 0;
        
        figure; 
        ptSrc = cell(1,ndims(costMatrix));        
        [ptSrc{:}] = ind2sub( size(costMatrix), srcNode );
        ptSrc = cell2mat(ptSrc);
        
        ptDest = cell(1,ndims(costMatrix));
        [ptDest{:}] = ind2sub( size(costMatrix), destNode );        
        ptDest = cell2mat(ptDest);
        
    end
    
    while ~pq.isEmpty()
        
        [curNodeIndex,~] = pq.pop();
        
        if flagDebugMode
            visitTime( curNodeIndex ) = visitTimer;
            visitTimer = visitTimer + 1;
            
            imshow(visitTime, []);
            hold on;
                plot(ptSrc(2), ptSrc(1), 'ro', 'LineWidth', 2.0);
                plot(ptDest(2), ptDest(1), 'ro', 'LineWidth', 2.0);
            hold off;
            
        end
        
        % terminate search if destination node has been reached
        if curNodeIndex == destNode
            break;
        end

        % prepare a list of its neighbors
        curNodeSubInd = cell(1, imdims);
        [curNodeSubInd{:}] = ind2sub(imsize, curNodeIndex);
        curNodeSubInd = cell2mat(curNodeSubInd);
        
        % update distances to all neighbors
        for i = 1:size( neighOffSetSubInd, 1 )
           
            % boundary neighbor check - if neighbor at offset is inside matrix
            curNeighSubInd = curNodeSubInd + neighOffSetSubInd(i,:);
            if any(curNeighSubInd < 1) || any(curNeighSubInd > imsize)
                continue;
            end
            curNeighSubInd = num2cell(curNeighSubInd);
            curNeighInd = sub2ind(imsize, curNeighSubInd{:});
            
            % ignore if it lies outside the roi mask
            if ~isempty(roiMask) && ~roiMask(curNeighInd)
                continue;
            end
            
            % update distance if a new minimum path is found            
            curNeighEucDist = norm(neighOffSetSubInd(i,:) .* spacing);
            curNeighEdgeWeight = costWeight * costMatrix(curNeighInd) + (1 - costWeight) * curNeighEucDist;
            
            if dist(curNeighInd) > dist(curNodeIndex) + curNeighEdgeWeight
                
                dist(curNeighInd) = dist(curNodeIndex) + curNeighEdgeWeight;                
                parentIndex(curNeighInd) = curNodeIndex;
                
                if pq.contains(curNeighInd)
                    pq.update(curNeighInd, dist(curNeighInd));
                else
                    pq.insert(curNeighInd, dist(curNeighInd));
                end
                    
            end
                
        end
        
    end
    
    % extract list of indices on the path using backtracking
    path = destNode;
    
    curNode = destNode;
    while curNode ~= srcNode        
        curNode = parentIndex( curNode );
        path = cat(1, path, curNode); 
    end
    
    path = flipud(path);
    cost = dist(destNode);
    
end
