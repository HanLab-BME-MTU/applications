function reorganize_oif_data( dspdatadir, datasearchdirlist, outdir )

    diary( fullfile( outdir, [mfilename, '.log'] ) );
    diary on;
    
    fprintf( '\n\n' );    
    process_dir( '', dspdatadir, datasearchdirlist, outdir );    
    fprintf( '\n\n' );
    
    diary off;
    
end

function process_dir( curDir, dspdatadir, datasearchdirlist, outdir )
    
    fprintf( '\n\nProcessing Directory - %s ... ', curDir );

    tiffilelist = dir( fullfile( dspdatadir, curDir, '*.tif' ) );
    tiffilelist = { tiffilelist.name };    
    
    for i = 1:numel( tiffilelist )
        
        [pathstr, name, ext] = fileparts( tiffilelist{i} );
        
        fprintf( '\n\n\t%d/%d Processing File - %s ... ', i, numel( tiffilelist ), tiffilelist{i} );
        
        tic
        blnFileFound = false;
        for j = 1:numel(datasearchdirlist)
            
            if exist( fullfile( datasearchdirlist{j}, [ name, '.oif' ] ), 'file' )
                
               blnFileFound = true;
               curOutDir = fullfile( outdir, curDir );
               
               if ~isdir( fullfile( curOutDir, 'dsp' ) )
                    mkdir( fullfile( curOutDir, 'dsp' ) );
               end
               
               copyfile( fullfile( dspdatadir, curDir, tiffilelist{i} ), fullfile(curOutDir, 'dsp') );
               copyfile( fullfile( datasearchdirlist{j}, [name, '.oif.files'] ), fullfile( curOutDir,  [name, '.oif.files'] ) );
               copyfile( fullfile( datasearchdirlist{j}, [name, '.oif'] ), curOutDir );               
               
               break;
               
            end
        end
        
        if ~blnFileFound
           fprintf( '\n\n\t\tCould not find the oif file - %s', tiffilelist{i} );
        else
            timeElapsed = toc;
            fprintf( 'took %.2f seconds', timeElapsed ); 
        end
        
    end
    
    dirInfo = dir( fullfile( dspdatadir, curDir ) );
    subdir = { dirInfo( [dirInfo.isdir] ).name };
    subdir = setdiff( subdir, { '.', '..' } );
    
    for i = 1:numel( subdir )        
        process_dir( fullfile(curDir, subdir{i}), dspdatadir, datasearchdirlist, outdir );        
    end
end