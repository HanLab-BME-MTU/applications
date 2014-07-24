function [imPalm, imNMap, imZMap] = getPalmImage( dataFilePath, varargin )

    p = inputParser;
    p.addRequired( 'dataFilePath', @(x) (ischar(x)) );
    p.addParamValue( 'magnification', 4, @(x) (isscalar(x)) );
    p.addParamValue( 'pixelSize', 133.33, @(x) (isscalar(x)) );
    p.parse( dataFilePath, varargin{:} );
    
    mag = p.Results.magnification;
    pixelSize = p.Results.pixelSize;
    
    % read data file and extract X,Y, and Z positions of all molecules
    fprintf( 1, '\n\tReading iPalm data file ...\n' );
    
        % read headers and find column ids of Group X, Y, and Z positions
        strposheaders = { 'Group X Position', 'Group Y Position', 'Group Z Position' };
        fid = fopen(dataFilePath);
        hdrline = fgetl(fid);
        fclose(fid);
        dataFileContents.colheaders = textscan( strtrim(hdrline), '%s', 'delimiter', sprintf('\t') );
        dataFileContents.colheaders = strtrim( dataFileContents.colheaders{1} );

        for i = 1:numel( strposheaders )
            curcolind = find( strcmpi( dataFileContents.colheaders, strposheaders{i} ) );        
            if isempty( curcolind ) || numel( curcolind ) > 1
                error( 'ERROR: could not read %s from data file', strposheaders{i} );
            end               
            poscolind(i) = curcolind;
        end
        
        % read numeric data
        dataFileContents.data = dlmread( dataFilePath, '\t', 1, 0 );
    
        % extract required columns
        MoleculePositionList = dataFileContents.data(:,poscolind);    
        clear dataFileContents;
        
    % generate palm image
    fprintf( 1, '\n\tGenerating Palm Images ...\n' );
    minPosPhysp = min( MoleculePositionList(:,1:2) );
    maxPosPhysp = max( MoleculePositionList(:,1:2) );
    
    minPosImsp = floor( minPosPhysp * mag / pixelSize );
    maxPosImsp = ceil( maxPosPhysp * mag / pixelSize );
    
    imsize = fliplr( maxPosImsp - minPosImsp + 1 );
    imPalm = zeros( imsize );
    imNMap = zeros( imsize );
    imZMap = zeros( imsize );
    
    MoleculePositionList(:,1) = round(MoleculePositionList(:,1) * mag / pixelSize) - minPosImsp(1) + 1;
    MoleculePositionList(:,2) = round(MoleculePositionList(:,2) * mag / pixelSize) - minPosImsp(2) + 1;
    
    for i = 1:size( MoleculePositionList, 1 )
        xind = MoleculePositionList(i,1);
        yind = MoleculePositionList(i,2);
        imPalm( yind, xind ) = 1;
        imNMap( yind, xind ) = imNMap( yind, xind ) + 1;
        imZMap( yind, xind ) = imZMap( yind, xind ) + MoleculePositionList(i,3);
    end    
    imZMap( imNMap > 0 ) = imZMap( imNMap > 0 ) ./ imNMap( imNMap > 0 );
    
end