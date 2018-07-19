function [ dicomRegister elastixTransforms itersInfo ] = elastixRegister(dicomFix, dicomMov, elastixParameters, varargin)
% ELASTIXREGISTER Perform a registration between a fixed volume and a moving volume, using Elastix.
% Author: Paul Balança
%
% [ dicomRegister elastixTransforms ] = ELASTIXREGISTER(dicomFix, dicomMov, elastixParameters, param1, val1, param2, val2, ...)
%
%     Input:
%        dicomFix            Dicomdata data containing fixed volume.
%        dicomMov            Dicomdata data containing moving volume.
%        elastixParameters   Parameters used by Elastix (can be load from a text file with readElastixParameters).
%
%     Parameters:
%        'logFile'           Filename where to save the log file.
%        'elastixExe'        Filename of the Elastix executable.
%        'outputDir'         Directory where to save Elastix outputs.
%                            The default value is a temporary directory which will be deleted at the end.
%        'initialTransform'  Initial transform to perform before registration.
%        'fixedMask'         Fixed mask use to compute the metric.
%        'movingMask'        Moving mask use to compute the metric.
%
%     Output: 
%        dicomRegister       Moving volume registered on the fixed one.
%        elastixTransforms   Transforms parameters produced by Elastix.
%        itersInfo           Convergence information at each iteration.
%
%  See Elastix documentation : http://elastix.isi.uu.nl
%

    % Parameters
    p = inputParser;
    p.addParamValue('logFile', '');
    p.addParamValue('elastixExe', '');
    p.addParamValue('outputDir', '');
    p.addParamValue('initialTransform', []);
    p.addParamValue('fixedMask', []);
    p.addParamValue('movingMask', []);
    p.parse(varargin{:});
    parameters = p.Results;

    % Temp directory, if no output directory
    if isempty(parameters.outputDir)
        if ~isempty(getenv('TEMP'))
            tmpOutDir = getenv('TEMP');
        elseif ~isempty( strfind( computer('arch'), 'glnxa' ) ) && isdir( '/tmp' )
            tmpOutDir = '/tmp';
        else
            tmpOutDir = pwd;            
        end
        c = clock;
        uniqDirSuffix = sprintf( '.%d_%d_%d_%f_%d', c(3), c(4), c(5), c(6), randi(100) );
        pathdir = fullfile(tmpOutDir, ['elastix.tmp', uniqDirSuffix] );
    else
        pathdir = parameters.outputDir;
    end
    if ~exist(pathdir, 'dir')
        mkdir(pathdir);
    end

    % Write dicom volumes
    fprintf('Write dicom volumes in MHD/RAW format ...\n');
    fnameFix = fullfile(pathdir, 'fixed' );
    mhdRawWrite(int16(dicomFix.im), dicomFix.Spacing, fnameFix);
    fnameMov = fullfile(pathdir, 'moving');
    mhdRawWrite(int16(dicomMov.im), dicomMov.Spacing, fnameMov);

    % Force write result image
    elastixParameters{numel(elastixParameters)}.WriteResultImage{1} = 'true';

    % Write parameters
    cmdParameters = '';
    fprintf('Write Elastix parameters files ...\n');
    for i = 1:numel(elastixParameters)
        fname = fullfile(pathdir, sprintf('parameters.%d.txt', i));
        writeElastixParameters(elastixParameters{i}, fname);
        cmdParameters = [cmdParameters '  -p "' fname '"']; %#ok<AGROW>
    end

    % Use initial transform
    if ~isempty(parameters.initialTransform)
        fname = fullfile(pathdir, 'InitialTransform.txt');
        writeElastixTransform(parameters.initialTransform, fname);
        cmdParameters = [cmdParameters '  -t0 "' fname '"']; %#ok<AGROW>
    end

    % Fixed mask
    if ~isempty(parameters.fixedMask)
        fname = fullfile(pathdir, 'fixedMask');
        mhdRawWrite(logical(parameters.fixedMask), dicomFix.Spacing, fname);
        cmdParameters = [cmdParameters '  -fMask "' fname '.mhd"']; %#ok<AGROW>
    end

    % Moving mask
    if ~isempty(parameters.movingMask)
        fname = fullfile(pathdir, 'movingMask');
        mhdRawWrite(logical(parameters.movingMask), dicomMov.Spacing, fname);
        cmdParameters = [cmdParameters '  -mMask "' fname '.mhd"']; %#ok<AGROW>
    end

    % Call Elastix ...
    if ~isempty(parameters.elastixExe)
        cmdElastix = parameters.elastixExe;
    else
        % M-file directory
        M = inmem('-completenames');
        for i = 1:numel(M)
            if ~isempty(strfind(M{i}, 'elastixRegister.m'))
                folderExe = fileparts(M{i});
            end
        end
        
       switch computer('arch')
           
           case { 'win32', 'win64'}
               
                cmdElastix = fullfile(folderExe, 'elastix.exe');
                
           case 'glnxa64' 
               
                cmdElastix = fullfile(folderExe, 'elastixLinux');
                
           case 'maci64'     
               
               cmdElastix = fullfile(folderExe, 'elastixMac');
          
           otherwise
           
               error( 'Incompatible computer architecture' );
       end
        
    end
    
    % make sure we have permission to execute elastix
    if ismember(computer('arch'), {'glnxa64', 'maci64'})

        try

            [status, message, messageid] = fileattrib([cmdElastix, '*'], '+x', 'u');

            if ~status
                error( 'Tried to assign executable permission but failed. Probably you dont have rights to do that. Contact your system administrator.' );
            end

        catch ME

            if strcmpi(computer('arch'), 'glnxa64')
                emsg = sprintf( 'Please make sure to give executable permission for files named elastixLinux and elastixLinuxExec in folder %s. Please contant the system administrator if you dont understand what this means.', folderExe );
            else
                emsg = sprintf( 'Please make sure to give executable permission for files named elastixMac and elastixMacExec in folder %s. Please contant the system administrator if you dont understand what this means.', folderExe );    
            end

            errordlg( emsg );
            error( emsg );        

        end

    end

    cmdElastix = [cmdElastix ' -f "' fnameFix '.mhd"'];
    cmdElastix = [cmdElastix ' -m "' fnameMov '.mhd"'];
    cmdElastix = [cmdElastix ' -out "' pathdir '"'];
    cmdElastix = [cmdElastix cmdParameters];

    cmdElastix    
    
    status = system(cmdElastix);

    % Error while registering
    if status ~= 0
       error('>>> Error (%d) while registering dicom volumes. <<<', status); 
    end

    % Load results
    fprintf('Load Elastix registered volume ...\n');
    fname = fullfile(pathdir, sprintf('result.%d.mhd', numel(elastixParameters)-1));
    dicomRegister = mhdRawRead(fname);

    % Load parameters
    fprintf('Load Elastix transforms ...\n');
    fname = fullfile(pathdir, sprintf('TransformParameters.%d.txt', numel(elastixParameters)-1));
    elastixTransforms = readElastixTransform(fname);

    % Load log
    if ~isempty(parameters.logFile)
        fprintf('Save Elastix log ...\n');
        copyfile(fullfile(pathdir, 'elastix.log'), parameters.logFile);
    end

    % Load iterations information
    fprintf('Load iterations information ...\n');
    for i = 1:numel(elastixParameters)
        for j = 1:elastixParameters{i}.NumberOfResolutions
            fname = sprintf('IterationInfo.%d.R%d.txt', i-1, j-1);
            fname = fullfile(pathdir, fname);
            itersInfo{i,j} = readIterationInfo(fname); %#ok<AGROW>
        end
    end

    % Remove temp files
    if isempty(parameters.outputDir)
        rmdir(pathdir, 's');
    end

end
