function [xlistList, sortNumbers] = readSynthXlist(listID, path, sortExpression, sortOrder)
% READSYNTHXLIST  reads synthetic slists or idlists into a structure.
%
% the output structure is sorted according to sortExpression and sortOrder,
% which are both optional
%
% listID is either "s" or "slist" or "idl" or "idlist" or "idt" or 
% "idlisttrack" It is assumed that the filenames are of the form 
% Xlist_?#_?# etc, where ? is a letter and # is a double
%
% sortexpression for slists with amplitude, SNR, and iteration:
% '_A([\d.]+)_S([\d.]+)_i([\d]+)' A, S, and i are the identifiers after
% which the doubles are written
%
% sortOrder gives the order of columns after which the columns are sorted.
% A negative column-index means that the rows are sorted in descending
% order
%


% TEST INPUt
if nargin < 2 || isempty(path)
    path = pwd;
end
if nargin < 3 || isempty(sortExpression)
    sortExpression = 0;
end
if nargin < 4 || isempty(sortOrder)
    sortOrder = 0;
end

% READ FILES
switch listID(1)
    case 's'
        listID = 'slist_';
    case 'i'
        if strcmp(listID,'idt') || strcmp(listID,'idlisttrack')
            listID = 'idlisttrack_';
        else
            listID = 'idlist_';
        end

    otherwise
        error('listID not recognized')
end

% search files. fileList(:,1) is a cell array with file names
% if path is a list of paths, search all of them
if iscell(path)
    % there are not going to be too many paths, so reassignment within the
    % loop is not totally bad
    fileList = {};
    for i=1:length(path)
        fileList = cat(1,fileList,searchFiles(listID,[],path{i},0));
    end
else
    fileList = searchFiles(listID,[],path,0);
end

% if required, files are sorted
if sortExpression
    
    % read numbers with sortExpression
    if strmatch('6.5',version) % ensure backwards compatibility
        [dummy,dummy,idxList] = regexp(fileList(:,1),sortExpression);
        if isempty(idxList)
            xlistList = [];
            disp('Warning: no files found')
            return
        end
        for i=length(idxList):-1:1
             tokenList = idxList{i};
             tokenList = tokenList{1};
             currentFile = fileList{i,1};
             for j=size(tokenList,1):-1:1
                 sortNumbers(i,j) = ...
                     str2double(currentFile(tokenList(j,1):tokenList(j,2)));
             end
         end
           % sort
    if ~sortOrder
        sortOrder = 1:size(sortNumbers,2);
    end
    % sortrows doesn't accept negative numbers in the sortOrder in 6.51
    sig = repmat(sign(sortOrder),size(sortNumbers,1),1);
    sortNumbers = sortNumbers .* sig;
    [sortNumbers,rowIdx] = sortrows(sortNumbers,abs(sortOrder));
    sortNumbers = sortNumbers.* sig;
      
                 
    else
        sortNumbers = regexp(fileList(:,1),sortExpression,'tokens');
        
        % transform strings to doubles. This requires that the same number of
        % elements are found for every filename
        % quit if empty
        if isempty(sortNumbers)
            xlistList = [];
            disp('Warning: no files found')
            return
        end
        sortNumbers = cat(1,sortNumbers{:});
        sortNumbers = cat(1,sortNumbers{:});
        sortNumbers = str2double(sortNumbers);
        % sort
    if ~sortOrder
        sortOrder = 1:size(sortNumbers,2);
    end
    [sortNumbers,rowIdx] = sortrows(sortNumbers,sortOrder);

    end

    
    % reorder fileList
    fileList = fileList(rowIdx,:);

end

% load all the files
for i=size(fileList,1):-1:1
    xlistList(i,1) = load([fileList{i,2},filesep,fileList{i,1}]);
end