function WriteFeatureMatrixToCSVFile( filepath, featureMatrix, featureNameList, featureClass, flagClassOnFront )

    if exist( 'featureClass', 'var' )       
        
        if ~exist( 'flagClassOnFront', 'var' )
            flagClassOnFront = true;
        end
        
        if flagClassOnFront
            featureNameList = cat(1, 'ClassLabel', featureNameList);
            featureMatrix = cat(2, featureClass, featureMatrix);        
        else
            featureNameList = cat(1, featureNameList, 'ClassLabel');
            featureMatrix = cat(2, featureMatrix, featureClass);        
        end
            
    end
    
    fid = fopen( filepath, 'w' );
    
    for i = 1:numel(featureNameList)       
        if i > 1
            fprintf( fid, ',' );
        end
        fprintf( fid, '%s', featureNameList{i} );
    end
    
    fprintf( '\nProgress: \n' );
    
    last_percent_done = 0;
    numPrint = 0;
    
    for i = 1:size(featureMatrix,1)       
        
        fprintf( fid, '\n' );
        for j = 1:size(featureMatrix,2)
           
            if j > 1
                fprintf( fid, ',' );
            end
                
            featureVal = featureMatrix{i,j};
            
            if isnumeric(featureVal) 
                
                if featureVal - floor(featureVal) > 0
                    fprintf( fid, '%f', featureVal );
                else
                    fprintf( fid, '%d', featureVal );
                end
                
            elseif islogical( featureVal )
                
                fprintf( fid, '%d', featureVal );
                
            elseif ischar( featureVal )
                
                fprintf( fid, '%s', featureVal );
                
            else
                fclose( fid );
                error( 'ERROR - invalid feature data type' );                
            end               
            
        end     
        
        percent_done = round(100*i/size(featureMatrix,1));       
        
        if percent_done > last_percent_done
            fprintf( '%.2d%%  ', percent_done );
            last_percent_done = percent_done;
            numPrint = numPrint + 1;
            if mod( numPrint, 10 ) == 0
               fprintf( '\n' ); 
            end
        end        
            
    end

    fprintf( '\n' );
    
    fclose( fid );

end