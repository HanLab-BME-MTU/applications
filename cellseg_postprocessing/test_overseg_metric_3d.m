function test_overseg_metric

    clc
    clear 
    close all;

    % load data
    defaultDataDir = 'Z:\intravital\data\Stefan_September2012 (Mouse 4)\imageanalysis\annotations\overseg';

    if ~exist( 'dataFilePath', 'var' )
        [fileName,pathName] = uigetfile( fullfile( defaultDataDir, '*.mat' ), 'Select the data file' );   
        dataFilePath = fullfile( pathName, fileName )
    end

    dataFileContents = load( dataFilePath );
    imInput = dataFileContents.imageData{1,1};

    % standardize the image
    imAdjusted = mat2gray( imInput ) * 4096;    

    % apply median filter to remove any spotty noise
    fprintf( '\n\n>> Appling median filter to remove noise ...\n\n' );

    imAdjusted = matitk( 'FMEDIAN', [1,1,1], imAdjusted );    

    % show seg result
    segMaskRGB = label2rgbND( dataFileContents.imLabelCellSeg, dataFileContents.CellSegColorMap );
    imseriesmaskshowrgb( imInput, segMaskRGB );

    % extract profiles for each cell
    adjMatrix = zeros( numel(dataFileContents.cellStats) * [1,1] );
    hProfileOverSeg = figure; title( 'Over segmented cells' );
    hProfileWellSeg = figure; title( 'Well segmented cells' );
    profileMatOverSeg = [];
    profileMatWellSeg = [];

    for i = 1:numel(dataFileContents.cellStats)

        curCellStats = dataFileContents.cellStats(i);

        fprintf( '\n\n>> Analyzing cell id - %d of type - %s\n\n', i, curCellStats.cellPatternType );

    %     % crop cell data
    %     ptCellPixel = cell(1, 3);
    %     [ptCellPixel{:}] = ind2sub( size(imInputimInput), curCellStats.PixelIdxList );
    %     ptCellPixel = cell2mat( ptCellPixel );
    %     
    %     ptTopLeft = max([ ones(1,3), min( ptCellPixel ) - 1]);
    %     ptBotRight = min([ size(imInput), max( ptCellPixel ) + 1]);
    %     
    %     indCellCrop = cell(1,3);
    %     for j = 1:3
    %         indCellCrop{j} = ptTopLeft(j) : ptBotRight(j);
    %     end
    %     
    %     imCurCellCropped = imAdjusted(indCellCrop{:});
    %     imCurCellMaskCropped = (dataFileContents.imLabelCellSeg(indCellCrop{:}) == i);

        % get cell mask
        imCurCellMask = (dataFileContents.imLabelCellSeg == i);

        % check if it has any other cells touching it
        imDilatedCellMask = imdilate(imCurCellMask, ones(5,5,5));
        neighCellIds = unique(dataFileContents.imLabelCellSeg( imDilatedCellMask > 0 ));
        neighCellIds = setdiff( neighCellIds, [i, 0] );

        if isempty( neighCellIds )
            fprintf( '\nisolated cell - no more analysis needed\n' );
        end

        % extract profiles for each neighboring cell
        indCurCellSeed = find( imCurCellMask & dataFileContents.imCellSeedPoints > 0 );
        ptCurCellSeed = cell(1,3);
        [ptCurCellSeed{:}] = ind2sub( size(imInput), indCurCellSeed(1) );
        ptCurCellSeed = cell2mat( ptCurCellSeed );

        for j = neighCellIds

            adjMatrix(i, j) = 1;

            if j < i 
                continue;
            end

            % get seedpoint location of neighboring cell
            indCurNeighCellSeed = find( dataFileContents.imLabelCellSeg == j & dataFileContents.imCellSeedPoints > 0 );
            ptCurNeighCellSeed = cell(1,3);
            [ptCurNeighCellSeed{:}] = ind2sub( size(imInput), indCurNeighCellSeed(1) );
            ptCurNeighCellSeed = cell2mat(ptCurNeighCellSeed);

            % extract profile
            vec = ptCurNeighCellSeed - ptCurCellSeed;
            t = linspace(0, 1.0, norm(vec) );
            ptLine = repmat(ptCurCellSeed, [numel(t), 1] ) + t(:) * vec;

            profileIm = interp3( imAdjusted, ptLine(:,2), ptLine(:,1), ptLine(:,3) );
            profileIm = mat2gray(profileIm);

            profileMask = interp3( dataFileContents.imLabelCellSeg, ptLine(:,2), ptLine(:,1), ptLine(:,3), 'nearest' );

            indChange1 = find( profileMask ~= i );
            indChange2 = find( profileMask == j );        
            indChangeMid = round( 0.5 * (indChange1(1) + indChange2(1)) );
            indLine = (1:numel(profileMask))' - indChangeMid;

            switch curCellStats.cellPatternType

                case 'Over_Segmentation'

                    figure(hProfileOverSeg);
                    hold all;
                    plot( indLine(indChangeMid:end), profileIm(indChangeMid:end) );
                    plot( -1 * indLine(indChangeMid:-1:1), profileIm(indChangeMid:-1:1) );
                    ylim( [0, 1.0] );
                    hold off;

                    profileMatOverSeg( end+1, indLine(indChangeMid:end) + 1 ) = profileIm(indChangeMid:end);
                    profileMatOverSeg( end+1, abs(indLine(indChangeMid:-1:1)) + 1) = profileIm(indChangeMid:-1:1);

                case 'Mono_G1'

                    figure(hProfileWellSeg);
                    hold all;
                    plot( indLine(indChangeMid:end), profileIm(indChangeMid:end) );
                    plot( -1 * indLine(indChangeMid:-1:1), profileIm(indChangeMid:-1:1) );
                    ylim( [0, 1.0] );
                    hold off;

                    profileMatWellSeg( end+1, indLine(indChangeMid:end) + 1 ) = profileIm(indChangeMid:end);
                    profileMatWellSeg( end+1, abs(indLine(indChangeMid:-1:1)) + 1) = profileIm(indChangeMid:-1:1);

            end

        end

    end
    
end





