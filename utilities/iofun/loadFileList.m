function fileListOut = loadFileList(inputChoicesCell,inputTypeList,fileListFile,helpTitle)
%LOADFILELIST allows a user to generate a fileList for subsequent automatic loading of data
%
%SYNOPSIS fileListOut = loadFileList(inputChoicesCell,inputTypeList,fileListFile)
%
%INPUT    if the first two input arguments are empty, the third must not be
%         inputChoicesCell : Cell containing the input choices being passed
%                            to the matlab function uigetfile
%         inputTypeList    : list with the same number of rows as
%                            inputChoicesCell and 1-2 columns. The first
%                            column specifies if the file is a data file
%                            (1) or if it is a file containing a fileList
%                            (2). The second column (optional) is an identification of
%                            the type of file that will be passed to the
%                            output.
%         fileListFile     : [pathName,fileName] of a text file containing
%                            a fileList that is to be loaded. The format of
%                            this file is as follows: The first line for every
%                            file is an identifier specifying a first
%                            part of the path that has been stored in an
%                            environment variable. If there is no such
%                            path, the string on the identifier line is
%                            'NONE'. If data has come from e.g a
%                            simulation, the identifier is 'NOFILE'
%                            Optional additional parameters are
%                            stored in the same line, separated by #. The
%                            expressions inbetween will be read as string:
%                            fileListOut.opt{i} = string.
%                            The options have to end with #!
%                            The second line specifies the rest of
%                            the path including the filename.
%                            The list terminates with three asteriks ***.
%                            Lines commented in matlab style will not be
%                            read.
%
%                            Example:
%
%                            % Explanation of file
%
%                            BIODATA#spb1#cen1#
%                            \Wildtype30C\1secMovies\WT_1sec_30C_12\WT_1sec_30C_12-data.mat
%                            ...
%                            %alternative forms (spaces are important!)
%                            BIODATA   \Wildtype30C\1secMovies\WT_1sec_30C_12\WT_1sec_30C_12-data.mat
%                            BIODATA#spb1#cen1#   \Wildtype30C\1secMovies\WT_1sec_30C_12\WT_1sec_30C_12-data.mat
%                            ...
%                            ***
%
%        helpTitle         : (opt) title of the help dialog box. Default: ''
%
%OUTPUT  fileListOut       : struct with fields .file, .type of length n,
%                            where n is the number of files to be loaded.
%                   .file  : full path for the file including filename
%                   .type  : (numeric) filetype as specified in the input.
%                   .opt   : cell array with additional information, as
%                            stored on the identifier-line. opt{1} is the
%                            identifier
%        if everything cancelled by user, fileListOut will be empty (actually, it will be a 1x0 structure)
%
%WARNING: because the function exist(filename,'file') is very slow, the
%         program does not check wheter the returned files actually exist.
%         You have to use load within a try-catch clause to make your
%         program work smoothly!
%
%c: 01/04 jonas dorn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function could be made faster by bypassing the input-check in textread

%----------TEST INPUT----------
allowedInputTypeListFirstCol = [1,2];

%test inputChoicesCell and inputTypeList first.
if nargin < 2
    error('at least two input arguments needed')
end

if isempty(inputChoicesCell)
    verbose = 0;
else
    verbose = 1;
    %we need both inputChoicesCell and inputTypeList to be ok
    if ~iscell(inputChoicesCell) | size(inputChoicesCell,2) ~= 2
        error('wrong list of choices - see help for builtin function uigetfile')
    end
    if  isempty(inputTypeList) | size(inputTypeList,1) ~= size(inputChoicesCell,1) | size(inputTypeList,2) > 2 | ~all(ismember(inputTypeList(:,1),allowedInputTypeListFirstCol))
        error('missing or badly specified inputTypeList')
    end
    if size(inputTypeList,2) == 1
        inputTypeList = [inputTypeList,zeros(size(inputTypeList))];
    end
end

if nargin < 3 | isempty(fileListFile)
    if verbose == 0
        error('if the first two input arguments are empty, fileListFile has to be specified')
    end
else
    if ~exist(fileListFile)
        error('can not find fileListFile')
    end
end

if nargin < 4 | isempty(helpTitle) | ~isstr(helpTitle)
    helpTitle = '';
end


%------END TEST INPUT----------

%--------ASK FOR FILES---------

%define variables
%output
store(1:100) = struct('file','','type',[],'toutDoux',[]);
nextStoreIdx = 1;
INITSTORELENGTH = 100;
storeLength = 100;

%while-loop var
done = 0;
%lookup list for what to do - add a 0th entry for cancel
toutDouxList = [-1;inputTypeList(:,1)];
choiceSequence = [1:size(inputTypeList,1)];
inputChoicesNow = inputChoicesCell;

%tell the rules to the user
h = helpdlg(...
    ['please select the filetypes and the files you want to load.',...
    ' Once you are done, press ''cancel'' in the selection dialogue'], helpTitle);
uiwait(h)

% remember where we were
oldPath = pwd;

%start while loop: loop until user does not want to load further data
while ~done

    %get file and entry number of typeList filterIndex returns the
    %choice or 0 if cancel
    [fileName,pathName,filterIdx] = uigetfile(inputChoicesNow,'select file or cancel!');

    %lookup the choice
    whatToutDoux = toutDouxList(filterIdx+1); %filterIdx could be 0 if cancel

    %switch according to the choice
    switch whatToutDoux
        %possible cases
        %-1: done = 1
        % 1: store file
        % 2: read files from list

        case -1
            %cancelled by user
            done = 1;

        case {1,2}
            %store the file. We will load files from the list later
            store(nextStoreIdx).file = [pathName, fileName];
            store(nextStoreIdx).type = inputTypeList(filterIdx,2);
            store(nextStoreIdx).toutDoux = whatToutDoux;


            %update nextOutIdx
            nextStoreIdx = nextStoreIdx + 1;

            %make sure we do not go over lenght of output
            if nextStoreIdx > storeLength
                storeLength = storeLength + INITSTORELENGTH;
                tmp = store;
                store(1:storeLength) = struct('file','','type',[],'toutDoux',[]);
                store(1:storeLength-INITSTORELENGTH) = tmp;
            end

            % change to the new directory, and go up if project dir (allow for nonBiodata)
            cd(pathName);
            cdBiodata(3);

        otherwise
            error('bad entry in chooseFileType was not correctly detected during input check')
    end

    % rearrange the input choices, so that the last selection is topmost
    if ~done
        % user chose choice#filteridx. copy #1 to that position and the
        % choice to position 1
        newTopChoice = choiceSequence(filterIdx);
        choiceSequence(filterIdx) = choiceSequence(1);
        choiceSequence(1) = newTopChoice;

        % rearrange lists
        inputChoicesNow = inputChoicesCell(choiceSequence,:);
        toutDouxList = [-1;inputTypeList(choiceSequence,1)];
    end

end %while ~done

% go back to old path
cd(oldPath);

store(nextStoreIdx:end) = [];

%----END ASK FOR FILES---------

%--------BUILD OUTPUT FILELIST----------

%define fileListOutTmp
fileListOutTmp(1:(nextStoreIdx + INITSTORELENGTH)) = struct('file','','type',[],'opt',[]);
nextFileIdx = 1;
fileListLength = length(fileListOutTmp);

%loop through store and build fileListOutTmp
for s = 1:nextStoreIdx - 1

    %switch according to field toutDoux
    switch store(s).toutDoux
        case 1
            %just take the file
            fileListOutTmp(nextFileIdx).file = store(s).file;
            fileListOutTmp(nextFileIdx).type = store(s).type;

            nextFileIdx = nextFileIdx + 1;

            %check that we have allocated enough space to fileListOutTmp.
            %otherwise, increase variable size
            if nextFileIdx > fileListLength
                fileListLength = fileListLength + INITSTORELENGTH;
                tmp = fileListOutTmp;
                fileListOutTmp(1:fileListLength) = struct('file','','type',[],'opt',{});
                fileListOutTmp(1:fileListLength-INITSTORELENGTH) = tmp;
            end

        case 2
            %load from fileList

            %read text with textread. Read the whole file (better than to read with fgetl)
            list = textread(store(s).file,'%s',-1,'commentstyle','matlab');

            %the list is now a cell array of strings, containing (maybe) empty cells at the beginning, then an indicator
            %(HOME,BIODATA,SIMDATA,NONE) followed by the relative path starting with the filesep
            %now loop through the list and recover the files
            done = 0;

            while ~done

                %if empty: drop entry. If three stars: done. Else, do
                %something

                if isempty(list{1})
                    %empty first entry
                    list(1) = [];

                else
                    %read the data
                    firstLine = list{1};
                    list(1) = [];

                    %firstLine is of the form
                    %IDENTIFIER#expr1#expr2# etc. we now look for all
                    %#-signs to read identifiers and the optional data,
                    %that will be read into opt

                    opt = {};
                    separatorIdx = findstr(firstLine,'#');
                    if separatorIdx(end) == length(firstLine);
                        separatorIdx = separatorIdx(1:end-1);
                    end

                    if isempty(separatorIdx)
                        %easy. only one entry for opt
                        identifier = firstLine;
                        opt{1} = identifier;
                    else
                        %assign identifier
                        identifier = firstLine(1:separatorIdx(1)-1);
                        opt{1} = identifier;

                        %loop through all separators and read the
                        %respective options
                        for sn = 1:length(separatorIdx)
                            opt{sn+1} = firstLine(separatorIdx(sn)+1:separatorIdx(sn+1)-1);
                        end

                    end % if isempty(separatorIdx)


                    firstPath = getenv(identifier);
                    isNone = strcmpi(identifier,'NONE');
                    isNofile = strcmpi(identifier,'NOFILE');

                    %check whether firstPath exists
                    if isempty(firstPath) & ~isNone & ~isNofile
                        noValidID = 1;
                        disp([identifier,' not found - ',list{1},' not loaded']);
                    else
                        noValidID = 0;
                    end



                    %read second part of path
                    secondPath = list{1};
                    list(1) = [];

                    %if there is nofile, we do not load: tell user
                    if isNofile
                        disp(['file ',secondPath,' not loaded: data was not saved']);
                    end

                    %if there was a valid identifier, we start with
                    %an filesep, else there is no point in changing
                    %the fileseps anyway: the drivenames are
                    %different between windows and linux

                    if ~isNone & ~isNofile & ~noValidID


                        %check first filesep: if it is not the correct one,
                        %replace all fileseps
                        if secondPath(1)~=filesep
                            %find all 'fileseps'
                            filesepList = strfind(secondPath,secondPath(1));
                            secondPath(filesepList) = filesep;
                        end
                    end

                    %store the file anyway (saves time)
                    file = [firstPath,secondPath];
                    if ~isNofile & ~noValidID
                        % automatically remove '_corr' from the fileName.
                        % If the file would really be '_corr', then it's
                        % high time the file was updated!
                        fileList(iFile).file = regexprep(fileList(iFile).file,'_corr','');
                        % store the file
                        fileListOutTmp(nextFileIdx).file = file;
                        fileListOutTmp(nextFileIdx).type = store(s).type;
                        fileListOutTmp(nextFileIdx).opt  = opt;

                        nextFileIdx = nextFileIdx + 1;

                        %check that we have allocated enough space to fileListOutTmp.
                        %otherwise, increase variable size
                        if nextFileIdx > fileListLength
                            fileListLength = fileListLength + INITSTORELENGTH;
                            tmp = fileListOutTmp;
                            fileListOutTmp(1:fileListLength) = struct('file','','type',[],'opt',{});
                            fileListOutTmp(1:fileListLength-INITSTORELENGTH) = tmp;
                        end
                        %                         else
                        %                             disp(['file ',file,' not found!']);
                        %                         end
                    end

                end %if isempty(list{1})

                %stop if there is no list left
                if isempty(list)|strmatch(list{1},'***')
                    done = 1;
                end
            end %while ~done
    end %switch store(s).toutDoux


end %for i = 1:nextStoreIdx - 1

%chop of surplus lenght of fileList
fileListOut = fileListOutTmp(1:nextFileIdx-1);


%----END BUILD OUTPUT FILELIST----------

