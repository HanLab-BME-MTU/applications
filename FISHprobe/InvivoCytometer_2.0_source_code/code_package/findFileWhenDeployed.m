function fpath = findFileWhenDeployed(fname, flagGenerateErrorIfNotFound)

    if ~exist('flagGenerateErrorIfNotFound', 'var')
        flagGenerateErrorIfNotFound = false;
    end
    
    if exist(fname, 'file')
        fpath = fname;
        return;
    end
        
    fpath = which(fname);
    if ~isempty(fpath)
        return;
    end
    
    if flagGenerateErrorIfNotFound
        error('Unable to find file - %s', fname);
    end
    
end