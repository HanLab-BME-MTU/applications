
costMatrix = zeros(100,100);

%costMatrix( 27:73, 27:73 ) = 100;  
srcNode = sub2ind(size(costMatrix), 25, 25);
destNode = sub2ind(size(costMatrix), 75, 75);

profile on

    [path, cost] = findMinimalPathOnCostMatrix(costMatrix, srcNode, destNode);
    pathMask = zeros( size(costMatrix) );
    pathMask( path ) = 1;
    
profile off
profile viewer

imseriesmaskshow( costMatrix, pathMask );
