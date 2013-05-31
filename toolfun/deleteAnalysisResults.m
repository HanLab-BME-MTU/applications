function deleteAnalysisResults(data)

for i = 1:numel(data)
    delpath = [data(i).source 'Detection'];
    if exist(delpath, 'dir')==7
        str = [];
        while ~any(strcmpi(str, {'y','n'}))
            str = input(['Delete results for ' getShortPath(data(i)) ' ? [y/n]: '],'s');
        end
        if strcmpi(str, 'y')
            rmdir(delpath, 's');
            fprintf('Deleted %s\n', delpath);
        end
    end
end