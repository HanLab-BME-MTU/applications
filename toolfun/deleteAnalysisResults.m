function deleteAnalysisResults(data, dirList)
if ~iscell(dirList);
    dirList = {dirList};
end

for i = 1:numel(data)
    for d = 1:numel(dirList)
        delpath = [data(i).source dirList{d}];
        if exist(delpath, 'dir')==7
            str = [];
            while ~any(strcmpi(str, {'y','n'}))
                str = input(['Delete ' getShortPath(data(i)) dirList{d} ' ? [y/n]: '],'s');
            end
            if strcmpi(str, 'y')
                rmdir(delpath, 's');
                fprintf('Deleted %s\n', delpath);
            end
        end
    end
end