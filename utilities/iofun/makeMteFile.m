function makeMteFile

% buld identifier cell
identifierCell = cell(3,3);

identifierCell{3,1} = 'HOME';
h = (getenv('HOME'));
if ~isempty(h) && strcmp(h(end),filesep) %should actually never happen
    h = h(1:end-1);
end
identifierCell{3,2} = h;
identifierCell{3,3} = length(identifierCell{3,2});

identifierCell{2,1} = 'BIODATA';
b = (getenv('BIODATA'));
if ~isempty(b) && strcmp(b(end),filesep)
    b = b(1:end-1);
end
identifierCell{2,2} = b;
identifierCell{2,3} = length(identifierCell{2,2});

identifierCell{1,1} = 'SIMDATA';
s = (getenv('SIMDATA'));
if ~isempty(s) && strcmp(s(end),filesep)
    s = s(1:end-1);
end
identifierCell{1,2} = s;
identifierCell{1,3} = length(identifierCell{1,2});

clear h b s

tagList = {'spb1';'cen1'};
writeShortFileList = 1;

% get directory
directory = uigetdir;

fileNameList = searchFiles('-data','log',directory,1,'new');

selectIdx = listSelectGUI(fileNameList(:,1),[],'move');
if isempty(selectIdx)
    return
else
    fileNameList = fileNameList(selectIdx,:);
end

[saveTxtName,saveTxtPath] = uiputfile({'*.mte','experimental MT data'},'save results as text file');


%if user cancelled, nothing will be saved
if saveTxtName == 0
    return
end


%create the file
    fidTxt = fopen([saveTxtPath,saveTxtName],'w');
    
    %if we selected to write the short version of the fileList, we do not
    %use a line break between identifier/options and the rest of the
    %filename
    if writeShortFileList
        separationString = '   ';
    else
        separationString = '\n';
    end
    
    %write introduction to file
    fprintf(fidTxt,'%s\n%s\n%s\n%s\n','%  MICROTUBULE DYNAMICS ANALYSIS - list of filenames',...
        '%    this file contains a list of filenames that can be used for MT dynamic analysis',...
        '%    the list is made up as: Identifier#{''tag1'',''tag2''}# \n  restOfPathIncludingFileName, where Identifier is the environment variable for the particular path',...
        '%    Please do not uncomment this header or delete the ''***'' that mark the end of the filenames. Fileseps can be windows or linux type.');
    
    
    %now loop through the fileNameList, find the identifier and write the
    %file
    for iFile = 1:size(fileNameList,1)
        
        %init variables
        identifier = '';
        restOfFileName = '';
       
        
        %read fileName
        fileName = [fileNameList{iFile,2},filesep,fileNameList{iFile,1}];
        lengthFileName = length(fileName);
        
        %now sieve the fileName until we find the identifier, or know for
        %sure that there isn't any
        %use if/elseif, because we want biodata/simdata to override home
        if ~isempty(identifierCell{1,2}) && strcmpi(identifierCell{1,2},fileName(1:min(identifierCell{1,3},lengthFileName))) %check for SIMDATA
            
            %read identifier, restOfFileName. There is no filesep at the
            %end of the identifier path, so the restOfFileName should start
            %with one
            identifier = identifierCell{1,1};
            restOfFileName = fileName(identifierCell{1,3}+1:end);
            
        elseif ~isempty(identifierCell{2,2}) && strcmpi(identifierCell{2,2},fileName(1:min(identifierCell{2,3},lengthFileName))) %check for BIODATA
            
            %read identifier, restOfFileName. There is no filesep at the
            %end of the identifier path, so the restOfFileName should start
            %with one
            identifier = identifierCell{2,1};
            restOfFileName = fileName(identifierCell{2,3}+1:end);
            
        elseif ~isempty(identifierCell{3,2}) && strcmpi(identifierCell{3,2},fileName(1:min(identifierCell{3,3},lengthFileName))) %check for HOME
            
            %read identifier, restOfFileName. There is no filesep at the
            %end of the identifier path, so the restOfFileName should start
            %with one
            identifier = identifierCell{3,1};
            restOfFileName = fileName(identifierCell{3,3}+1:end);
            
        elseif exist(fileName,'file') %check for NONE
            
            %assign none
            identifier = 'NONE';
            restOfFileName = fileName;
            
        else %we have no valid filename at all
            
            identifier = 'NOFILE';
            restOfFileName = fileName;
            
        end %check for identifier and restOfFilename
        
        %now write everything to file
        fprintf(fidTxt,['%s#%s#%s',separationString,'%s\n'],...
            identifier, tagList{1}, tagList{2}, restOfFileName);
        
        
    end %for nFile = 1:length(fileNameList)
    
    %close off the file by writing '***'
    fprintf(fidTxt,'\n***\n\n');
    fclose(fidTxt);


