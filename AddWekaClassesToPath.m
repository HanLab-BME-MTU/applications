function AddWekaClassesToPath( wekaPackageList )
    
    if ~exist( 'wekaPackageList', 'var' )
        wekaPackageList = { 'weka' };
    else
        wekaPackageList = [ {'weka'}, wekaPackageList ];
    end
    
    % add classes
    for i = 1:numel(wekaPackageList)
        
        curJarFileName = [wekaPackageList{i} '.jar'];
        if ~isempty( cell2mat( strfind( javaclasspath('-dynamic'), curJarFileName) ) )
            continue;
        end
        
        packageMatlabPath = which( curJarFileName );
        if ~isempty( packageMatlabPath )
            javaaddpath( packageMatlabPath );
            continue;
        end       
            
        switch wekaPackageList{i}
           
            case 'weka'
            
                winPackageRootDir = 'C:\Program Files\Weka-3-7';
                
            otherwise
                
                winPackageRootDir = fullfile( 'C:\Program Files\Weka-3-7', 'packages', wekaPackageList{i}, 'lib' );
        end
        
        if isdir(winPackageRootDir)
            
            curJarFileLocList = rdir( fullfile(winPackageRootDir, '**', curJarFileName) );
        
            if ~isempty(curJarFileLocList)
                javaaddpath( curJarFileLocList(1).name );
                continue;
            end      
        
        end
        
        error( 'Warning: unable to find the package %s. Place it in the current directory or in matlab path', curJarFileName );        
        
    end
end