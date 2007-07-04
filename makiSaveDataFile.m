function success = makiSaveDataFile(dataStruct,secureName)
%MAKISAVEDATAFILE saves maki-data to disk
%
% SYNOPSIS: success = makiSaveDataFile(dataStruct)
%
% INPUT dataStruct: maki-data structure (see makiMakeDataStruct)
%       secureName: fieldname of data item that should be saved via
%                   secureSave, such as 'initCoord' if the initial
%                   coordinate detection has just been performed.
%                   SecureSave ensures that old data will not  
%                   be overwritten
%
% OUTPUT success: 1 if successful save of data
%
% REMARKS
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn
% DATE: 28-Jun-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check for dataStruct
if nargin == 0 || isempty(dataStruct) || ~isstruct(dataStruct)
    error('no dataStruct supplied')
end
if nargin < 2 || isempty(secureName)
    secureName = 'noName';
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
                    % in case there was an older version, remember the name.
                    dataStruct.(fileName) = savedName;
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
    elseif false
    end

    % indicate that dataStruct has been saved in status (and last save)
    dataStruct.status(2) = 1;
    dataStruct.statusHelp{2,2} = date;

    % save dataFile
    save(fullfile(dataStruct.dataFilePath,dataStruct.dataFileName),'dataStruct');

catch
    if nargout > 0
        success = false;
    end
    return
end

if nargout > 0
    success = true;
end