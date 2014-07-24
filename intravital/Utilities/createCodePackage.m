function createCodePackage( mainFuncName, varargin )

    p = inputParser;
    p.addParamValue('outDir', [], @ischar);
    p.addParamValue('additionalMatlabFiles', {}, @iscell);
    p.addParamValue('additionalNonMatlabFiles', {}, @iscell);
    p.parse( varargin{:} );
    PARAMETERS = p.Results;

    outDir = PARAMETERS.outDir;
    if isempty( outDir )
        outDir = uigetdir;
    end
    
    fid = fopen( fullfile(outDir, 'readme.txt'), 'w' );
    
    mainFuncNameWithoutExt = strrep(mainFuncName, '.m', '');
    fprintf(fid, 'This folder contains the source code for the software %s\n',mainFuncNameWithoutExt );
    
    % analyze program dependencies dependencies
    PrettyPrintStepDescription( 'Analyzing Program Dependencies' );
    p = fdep( mainFuncName );

    for i = 1:numel(PARAMETERS.additionalMatlabFiles)
        curDep = fdep(PARAMETERS.additionalMatlabFiles{i});
        p.fun = union(p.fun, curDep.fun);
        p.toolbox = union(p.toolbox, curDep.toolbox);
    end
    
    % copy matlab files to output folder    
    PrettyPrintStepDescription( 'Copying matlab files to output folder' );    
    h = waitbar(0, 'Copying matlab files to output folder');
    
    for i = 1:numel(p.fun)        
        
        copyfile( p.fun{i}, outDir );
        
        % look for corresponding fig file
        figFileName = strrep(p.fun{i}, '.m', '.fig');
        if exist( figFileName, 'file' )
            copyfile( figFileName, outDir );
        end

        % look for mex files
        mexFileList = rdir( strrep(p.fun{i}, '.m', '.mex*') );
        if ~isempty(mexFileList)
            for j = 1:numel(mexFileList)
                copyfile( mexFileList(j).name, outDir );
            end
        end
        
        waitbar(i / numel(p.fun), h);
        
    end
    closeStatusDialog(h);

    % copy non-matlab files to output folder    
    if ~isempty( PARAMETERS.additionalNonMatlabFiles )
        
        PrettyPrintStepDescription( 'Copying non-matlab files to output folder' );    
        h = waitbar(0, 'Copying non-matlab files to output folder');

        for i = 1:numel(PARAMETERS.additionalNonMatlabFiles)        

            curFilePath = which( PARAMETERS.additionalNonMatlabFiles{i} );

            if numel(curFilePath) == 0
                fprintf('\nCould not find specified dependency file - %s\n', ...
                        PARAMETERS.additionalNonMatlabFiles{i} );
            else
                copyfile(curFilePath, outDir );
            end

            waitbar(i / numel(PARAMETERS.additionalNonMatlabFiles), h);
        end
        
    end    
    closeStatusDialog(h);
    
    % generate a readme file
    fprintf(fid, '\nThis software was developed using MATLAB %s\n', version );
    
    fprintf(fid, '\nThe following matlab toolboxes are required to run this software:\n' );
    for i = 1:numel(p.toolbox)
        fprintf(fid, '\n%s', p.toolbox{i});
    end
    fprintf(fid, '\n\nPlease make sure that the aforementioned toolboxes are available\n' );

    fprintf(fid, '\nUse the following steps to run the software: \n' );
    fprintf(fid, '\n1) Start matlab' );
    fprintf(fid, '\n2) Set the folder containing the code as the current working directory' );
    fprintf(fid, '\n3) Type %s and press the ENTER key', mainFuncNameWithoutExt );
    
    fclose(fid);   
    
end