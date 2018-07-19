function [nodeLocations, edges, edgeType, edgeTypeOffset] = buildGridGraphND( imsize, neighconn )

    imdims = numel(imsize);
    pixInd = (1:prod(imsize))';
    pixSubInd = cell(1, imdims);
    [pixSubInd{:}] = ind2sub(imsize, pixInd);
    
    switch neighconn
        
        case 0  % 4-connected (2D), 6-connected (3D), or 2N-connected (ND)
            
            edges = [];
            edgeType = [];
            edgeTypeOffset = eye(imdims);
            
            for i = 1:imdims
                
                % from node indices
                from_node_ind = pixInd( pixSubInd{i} + 1 <= imsize(i) );

                % to node indices
                neighPixSubInd = cell(1, imdims);
                for j = 1:imdims                    
                    neighPixSubInd{j} = pixSubInd{j}( from_node_ind );
                end
                neighPixSubInd{i} = neighPixSubInd{i} + 1;                
                
                to_node_ind = sub2ind(imsize, neighPixSubInd{:});                
                edges = cat(1, edges, [from_node_ind, to_node_ind]);
                
                edgeType = cat(1, edgeType, i + zeros(size(from_node_ind)));
                
            end  
            
            
        case 1  % 8-connected (2D), 26-connected (3D), or (3^N-1)-connected (ND)
        
            neighRasterVisitTimes = reshape( 1:3^imdims, 3 * ones(1,imdims) );
            neighOffSetInd = find( neighRasterVisitTimes > (3^imdims + 1)/2 );
            neighOffSetSubInd = cell(1, imdims);
            [neighOffSetSubInd{:}] = ind2sub( size(neighRasterVisitTimes), neighOffSetInd );
            neighOffSetSubInd = cell2mat( neighOffSetSubInd ) - 2;
            
            edges = [];
            edgeType = [];
            edgeTypeOffset = neighOffSetSubInd;
            
            for i = 1:numel( neighOffSetInd )
                
                % find pixels for which neighbor determined by offset lies
                % inside image
                flagNeighInside = true( size(pixInd) );
                for j = 1:imdims   
                    curNeighSubInd = pixSubInd{j} + neighOffSetSubInd(i,j);
                    flagNeighInside(curNeighSubInd < 1 | curNeighSubInd > imsize(j)) = false;
                end
                
                % from node indices
                from_node_ind = pixInd( flagNeighInside );
                
                % to node indices
                neighPixSubInd = cell(1, imdims);
                for j = 1:imdims   
                    neighPixSubInd{j} = pixSubInd{j}( flagNeighInside ) + neighOffSetSubInd(i,j);
                end
                
                to_node_ind = sub2ind(imsize, neighPixSubInd{:});                
                edges = cat(1, edges, [from_node_ind, to_node_ind]);
                
                edgeType = cat(1, edgeType, i + zeros(size(from_node_ind)));
                
            end
            
        otherwise
            
            error( 'ERROR: neighconn argument invalid - must be either 0 or 1' );
    end
    
    nodeLocations = cell2mat(pixSubInd) - 1;

end