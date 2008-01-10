function [success,dataStructOut] = makiSaveDataFile(serverType,dataStruct,secureName)
%MAKISAVEDATAFILE saves maki-data to disk
%
% SYNOPSIS: [success,dataStructOut] = makiSaveDataFile(serverType,dataStruct,secureName)
%
% INPUT serverType: 'TEST', 'HERCULES', 'DANUSER', 'MERALDI', 'SWEDLOW' or
%                   'MCAINSH'
%       dataStruct: maki-data structure (see makiMakeDataStruct)
%       secureName: fieldname of data item that should be saved via
%                   secureSave, such as 'initCoord' if the initial
%                   coordinate detection has just been performed.
%                   SecureSave ensures that old data will not
%                   be overwritten
%                   Pass multiple fieldNames as a cell array
%
% OUTPUT success: 1 if successful save of data
%        dataStructOut: updated dataStruct -- which contains the secure-saved
%                    filenames
%
% REMARKS If the dataFilePath does not exist, makiSaveDataFile attempts to
%         create it
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn
% DATE: 28-Jun-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check for dataStruct
if nargin < 2 || isempty(dataStruct) || ~isstruct(dataStruct)
    error('no dataStruct supplied')
end
if nargin < 3 || isempty(secureName)
    secureName = 'noName';
end
if ischar(secureName)
    secureName = {secureName};
end

if nargout == 2
    % make a copy of the incoming dataStruct onto dataStructOut; 
    % this is necessary because the incoming dataStruct is emptied of
    % actual data so that only minimal information is saved to the data file on disk;
    % however, functions, which retrieve from this function a dataStruct with the
    % updated secure-saved filenames, expect the full dataStruct returned.
    returnUpdatedStruct = 1;
    dataStructOut = dataStruct;
else
    returnUpdatedStruct = 0;
end

% check for dataFilePath and create if necessary
if ~isdir(dataStruct.dataFilePath)
    warning('MAKISAVEDATASTRUCT:PATHNOTFOUND',...
        'path %s not found. Folder will be created.',...
        dataStruct.dataFilePath)
    try
        mkdir(dataStruct.dataFilePath)
    catch
        warning('MAKISAVEDATASTRUCT:PATHNOTFOUND',...
            'unable to create directory %s ',...
            dataStruct.dataFilePath)
        if nargout > 0
            success = false;
        end
        return
    end
end


try
    % save individual files, such as dataProperties
    fn = fieldnames(dataStruct);

    for i=1:length(fn)
        % only save the xxx/xxxName pairs, and only if they're not empty
        fileName = [fn{i},'Name'];
        if any(strmatch(fileName,fn)) && ~isempty(dataStruct.(fn{i}))
            eval([fn{i},'=dataStruct.',fn{i},';']);

            switch fn{i}
                case secureName
                    % do not overwrite files we just changed
                    savedName = secureSave(...
                        fullfile(dataStruct.dataFilePath,dataStruct.(fileName)),fn{i});
                    % in case there was an older version, update the name.
                    dataStruct.(fileName) = savedName;
                    if returnUpdatedStruct
                        dataStructOut.(fileName) = savedName;
                    end
                otherwise
                    % overwrite, b/c these didn't/shouldn't change
                    save(...
                        fullfile(dataStruct.dataFilePath,dataStruct.(fileName)),fn{i});
            end
            % empty the xxx-fields here to save space & time
            dataStruct.(fn{i}) = [];
        end
    end

    % check rawMoviePath
    % if rawMoviePath is the same as dataFilePath, make it empty
    % if rawMoviePath contains some path stored in an environment variable
    % or similar, replace that part with an identifier
    if strmatch(dataStruct.rawMoviePath,dataStruct.dataFilePath)
        dataStruct.rawMoviePath = [];
    else
        try
            % make platform independant, if possible
            dataStruct.rawMoviePath = makiPathDef(dataStruct.rawMoviePath,serverType);
        catch
        end
    end

    % indicate that dataStruct has been saved in status (and last save)
    dataStruct.status(2) = 1;
    dataStruct.statusHelp{2,2} = date;

    % save dataFile
    save(fullfile(dataStruct.dataFilePath,dataStruct.dataFileName),'dataStruct');

catch
    err = lasterror;
    warning('MAKISAVEDATAFILE:FAILEDSAVE',err.message);
    if nargout > 0
        success = false;
    end
    return
end

if nargout > 0
    success = true;
end
