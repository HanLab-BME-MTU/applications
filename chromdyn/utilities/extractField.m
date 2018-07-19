function data=extractField(field,fileStringInclude,fileStringExclude,dim)
%EXTRACTFIELD extracts the contents of a specific field from a set of stored structures 
%
% SYNOPSIS data=extractField(fieldname,fileStringInclude,fileStringExclude,dim)
%
% INPUT    field             : name of the field to be displayed/returned
%          fileStringInclude : strings included in the filename of the
%                              saved structures
%          fileStringExclude : (opt) strings not to be included in the
%                              filename of the saved structures
%
% OUTPUT   data, structure with entries for every found value
%
%c: jonas 05/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=============
% TEST INPUT
%=============

% not finished yet

if nargin < 2
    error('not enough input arguments')
end

if nargin < 3 | isempty(fileStringExclude)
    fileStringExclude = '';
end

%=============

%=======================
% COLLECT LIST OF FILES
%=======================

listOfFiles = searchFiles(fileStringInclude,fileStringExclude,'ask',1,'all');

%=======================

%==================
% COLLECT FIELDS
%==================

nFilesMax = size(listOfFiles,1);
data(1:nFilesMax,1) = struct(field,[],'fileName',[]);
fileCt = 0;

for iFile = 1:nFilesMax
    
    % load stuff
    currentDir = [listOfFiles{iFile,2}, filesep];
    currentDataFile = load([currentDir,listOfFiles{iFile,1}]);
    
    % find the field we're looking for
    filesInData = fieldnames(currentDataFile);
    nVars = length(filesInData);
    iVar = 1;
    while iVar < nVars
        isRightVar = eval(['isfield(currentDataFile.' filesInData{iVar} ', ''' field ''');']);
        if isRightVar
            
            % update counters
            
            fileCt = fileCt + 1;
            
            % store data
            eval(['data(fileCt).' field '= currentDataFile.' filesInData{iVar} '.' field ';']);
            data(fileCt).fileName = [currentDir,listOfFiles{iFile,1},'   ' filesInData{iVar}];
            
            % disp
            
            % update the other counter
            iVar = nVars+1;
        else
            iVar = iVar + 1;
        end
            
    end
end

% cleanup
data(fileCt+1:end) = [];