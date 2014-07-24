% copy_oif_data( 'Z:\intravital\data\Stefan_Mouse8\goodones',
% 'Z:\intravital\data\Stefan_Mouse8\original', '..' )
function copy_oif_data( dspRootDir, oifDataRootDir, outDirRelativePath )

    if nargin < 1
        
        dspRootDir = uigetdir( 'Z:\', 'Enter root dir for dsp files' );
        oifDataRootDir = uigetdir( fullfile(dspRootDir, '..'), 'Enter root dir for oif files' );
        
    end

    if nargin < 4
        outDirRelativePath = '..';
    end
    
    diary( fullfile( dspRootDir, [mfilename, '.log'] ) );
    diary on;
    
    dspRootDir
    oifDataRootDir
    outDirRelativePath

    % make a list of all oif files
    fprintf( '\n\nPreparing a list of all oif files present ... \n\n' );
    
    oifFileList = rdir( fullfile(oifDataRootDir, '**', '**', '*.oif') );
    
    oifFileDirHashMap = containers.Map('KeyType', 'char', 'ValueType', 'char');    
    for i = 1:numel(oifFileList)
        
        [pathstr, name, ext] = fileparts( oifFileList(i).name );        
        oifFileDirHashMap( name ) = pathstr;
        
    end
    
    fprintf( '\n\n\t%d oif files were found ... \n\n', numel(oifFileList) );
        
    % get a list of all dsp files
    dspFileList = rdir( fullfile(dspRootDir, '**', 'dsp*\*.tif') ); 
    
    fprintf( '\n\nCopying the oif data of each dsp file selected ... ' );
    
    numFilesNotFound = 0;
    numCopyErrors = 0;
    
    for i = 1:numel(dspFileList)
        
        [pathstr, name, ext] = fileparts( dspFileList(i).name );
        
        fprintf('\n\n\t%d/%d (%.2f%%) Processing File - %s ... ', ...
                 i, numel( dspFileList ), i * 100 / numel(dspFileList), dspFileList(i).name);
        
        % find oif file
        if ~isKey(oifFileDirHashMap, name) 
           fprintf( '\n\n\t\tCould not find the oif file - %s', [name, '.oif'] );
           numFilesNotFound = numFilesNotFound + 1;
        else            
            try
                
            tic            
            curOifDir = oifFileDirHashMap(name);
            curOutDir = fullfile(pathstr, outDirRelativePath);
            
            copyfile( fullfile( curOifDir, [name, '.oif.files'] ), fullfile( curOutDir,  [name, '.oif.files'] ) );
            copyfile( fullfile( curOifDir, [name, '.oif'] ), curOutDir );                           
            
            timeElapsed = toc;
            fprintf( 'took %.2f seconds', timeElapsed ); 
            
            catch err
                toc
                fprintf( '\n\n\t\tCOPY-ERROR - Could not copy the oif file - %s', [name, '.oif'] );
                fprintf( '\n\t\terror description - %s', err.message ); 
                numCopyErrors = numCopyErrors + 1;
            end            
        end
        
    end
    
    numFilesNotFound
    numCopyErrors
    
    diary off;

end