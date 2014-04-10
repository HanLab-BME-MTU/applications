function [ L, adhesionStats ] = segAdhesionsBlobColoring( imPalm, adhesionAreaThreshold )

    imInput = double( imPalm > 0 );
    
    % preprocessing
    imFiltered = medfilt2( imInput, [5,5] );    
    imFiltered = imclose( imFiltered, strel( 'disk', 3 ) );
    imFiltered = imfill( imFiltered );
    
    imseriesshow( imFiltered );
    set( gcf, 'Name', 'Preprocessing Result' );    
    
    % connected component analysis
    L = bwlabel( imFiltered > 0, 8 );
    stats = regionprops( L, 'Area' );
    AreaVec = [stats.Area];

    imArea = zeros(size(L));
    imArea(L > 0) = AreaVec(L(L>0));
    imseriesmaskshow( imArea, imArea > adhesionAreaThreshold )
    set( gcf, 'Name', 'AdhesionAreaMap' );

    smallObjInd = find( AreaVec < adhesionAreaThreshold );
    L( ismember(L, smallObjInd) ) = 0;
    imAdhesionMask = double( L > 0 );

    L = bwlabel( imAdhesionMask > 0, 8 );
    [imSegRGBMask, labelColorMap] = label2rgbND( L );
    adhesionStats = regionprops( L, { 'Area', 'ConvexHull', 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation' } );

    % get cell boundary
    [yind, xind] = find( imAdhesionMask > 0 );    
    [yindOrg, xindOrg] = find( medfilt2( imInput, [3,3] ) > 0 );
    ptBottomLeft = [ min(xindOrg), max(yindOrg) ];
    ptBottomright = [ max(xindOrg), max(yindOrg) ];
    ptCell = [ [xind, yind]; ptBottomLeft ; ptBottomright ];
    k = convhull( ptCell );
    ptBoundary = [ptCell(k,1), ptCell(k,2)];    
    
    % display result
    imseriesmaskshowrgb( imPalm, imSegRGBMask );
    hold on;
    for i = 1:max(L(:))
       ptCenter = adhesionStats(i).Centroid;
       a = 0.5 * adhesionStats(i).MajorAxisLength;
       b = 0.5 * adhesionStats(i).MinorAxisLength;
       theta = adhesionStats(i).Orientation;
       ptRect = ([ -a, -b; ...
                    a, -b; ...
                    a, b; ...
                   -a, b; ...
                   -a, -b ])';
       ptRect(end+1,:) = 1;
       transMat = [ cosd(theta), sind(theta), ptCenter(1); ...
                   -sind(theta) cosd(theta), ptCenter(2); ...
                    0, 0, 1 ];
       ptRect = round( transMat * ptRect );
       plot( ptRect(1,:), ptRect(2,:), '-', 'Color', labelColorMap(i,:), 'LineWidth', 2.0 ); 
       
       plot( ptCenter(1), ptCenter(2), 'o', 'Color', [0,0,0], 'MarkerFaceColor', labelColorMap(i,:), 'LineWidth', 2.0 );        
    end
    plot( ptBoundary(:,1), ptBoundary(:,2), 'r-', 'LineWidth', 2.0 );
    hold off;
    set( gcf, 'Name', 'Adhesion Segmentation Result' );

end