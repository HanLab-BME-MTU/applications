function [ imAdhesionSegLabel, adhesionStats ] = segAdhesionsClustering( imPalm, varargin )

    p = inputParser;
    p.CaseSensitive = false;
    p.addRequired( 'imPalm', @(x) (isnumeric(x) && ndims(x) == 2) );
    p.addParamValue( 'minAdhesionArea', 1, @isscalar );
    p.addParamValue( 'sparsityFactor', 2.0, @isscalar );
    p.addParamValue( 'minClusterDistance', 30, @(x) isscalar(x) );
    p.addParamValue( 'flagComputeCellBoundary', false, @(x) (islogical(x) && isscalar(x)) );
    p.addParamValue( 'flagShowBBox', false, @(x) (islogical(x) && isscalar(x)) );
    p.addParamValue( 'flagDebugMode', true, @(x) (islogical(x) && isscalar(x)) );
    p.parse( imPalm, varargin{:} );
    
    sparsityFactor = p.Results.sparsityFactor;
    minAdhesionArea = p.Results.minAdhesionArea;
    minClusterDistance = p.Results.minClusterDistance;
    flagComputeCellBoundary = p.Results.flagComputeCellBoundary;
    flagShowBBox = p.Results.flagShowBBox;
    flagDebugMode = p.Results.flagDebugMode;
    
    %imInput = double( imPalm > 0 );
    
    % preprocessing - get rid of sparse isolated points   
    imFiltered = double( imfilter( imPalm, fspecial('average', [5, 5]) ) > sparsityFactor & imPalm > 0 ); 
    
    imseriesshow( imFiltered );
    set( gcf, 'Name', 'Preprocessing Result' );    
    
    % run mean-shift    
    bandwidth = minClusterDistance / 2.0;
    ptInd = find( imFiltered );
    [yind, xind] = ind2sub( size( imFiltered ), ptInd );    
    [clusterInfo, pointToClusterMap] = MeanShiftClustering([xind, yind], bandwidth, ... 
                                                         'flagDebug', true, ...
                                                         'kernel', 'gaussian', ...
                                                         'method', 'optimized', ...
                                                         'minClusterDistance', minClusterDistance, ...
                                                         'flagUseKDTree', true );

    % run adaptive mean-shift    
%     maxBandwidth = minClusterDistance / 2.0;
%     ptInd = find( imFiltered );
%     [yind, xind] = ind2sub( size( imFiltered ), ptInd );    
%     [clusterInfo, pointToClusterMap] = AdaptiveMeanShiftClustering([xind, yind], 20, ...
%                                                                    'maxBandwidth', maxBandwidth, ... 
%                                                                    'flagDebug', true, ...
%                                                                    'kernel', 'gaussian', ...
%                                                                    'method', 'optimized', ...
%                                                                    'minClusterDistance', minClusterDistance, ...
%                                                                    'flagUseKDTree', true );

    % delete clusters smaller than a specified threshold
    numSignificantClusters = 0;
    flagSmallCluster = false( numel(clusterInfo), 1 );
    for i = 1:numel( clusterInfo )
        if clusterInfo(i).numPoints < minAdhesionArea              
            flagSmallCluster(i) = true;
            pointToClusterMap( clusterInfo(i).ptIdData ) = 0;
        else
            numSignificantClusters = numSignificantClusters + 1;
            pointToClusterMap( clusterInfo(i).ptIdData ) = numSignificantClusters;
        end
    end
    clusterInfo( flagSmallCluster ) = []; % delete small clusters    
                                                     
    imAdhesionSegLabel = zeros( size( imPalm ) );
    imAdhesionSegLabel(ptInd) = pointToClusterMap;
    adhesionStats = regionprops( imAdhesionSegLabel, {'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'ConvexHull', 'ConvexArea'} );
    
    % get cell boundary
    if flagComputeCellBoundary
        
        [yind, xind] = find( imAdhesionSegLabel > 0 );
    %     [yindOrg, xindOrg] = find( medfilt2( imInput, [3,3] ) > 0 );
    %     ptBottomLeft = [ min(xindOrg), max(yindOrg) ];
    %     ptBottomright = [ max(xindOrg), max(yindOrg) ];
    %     ptCell = [ [xind, yind]; ptBottomLeft ; ptBottomright ];
        ptCell = [ [xind, yind] ];
        k = convhull( ptCell );
        ptBoundary = [ptCell(k,1), ptCell(k,2)];

    end
    
    % display area map
    imArea = zeros( size(imPalm) );
    AreaVec = [adhesionStats.Area];
    imArea(imAdhesionSegLabel > 0) = AreaVec(imAdhesionSegLabel(imAdhesionSegLabel>0));
    imseriesshow( imArea );
    set( gcf, 'Name', 'AdhesionAreaMap' );    
    
    % display segmentation result
    if flagDebugMode
        
        [imSegRGBMask, labelColorMap] = label2rgbND( imAdhesionSegLabel );   
        
        imseriesmaskshowrgb( imPalm, imSegRGBMask );
        set( gcf, 'Name', 'Adhesion Segmentation Result' );
        
        imseriesmaskshowrgb( imPalm, imSegRGBMask );
        hold on;
        
        if flagShowBBox
            
            for i = 1:max(imAdhesionSegLabel(:))
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
            
        end
        
        if flagComputeCellBoundary
            plot( ptBoundary(:,1), ptBoundary(:,2), 'r-', 'LineWidth', 2.0 );
        end
        
        hold off;
        set( gcf, 'Name', 'Adhesion Segmentation Result with Boundinng Boxes' );
        
    end

end

