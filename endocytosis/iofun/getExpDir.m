%[fpath] = getExpDir(data) returns the parent directory of the movies in 'data'

function [fpath] = getExpDir(data)

nd = numel(data);
fpath = cell(1,nd);
for k = 1:nd
    [~,fpath{k}] = getCellDir(data(k));
end
fpath = unique(fpath);
if numel(fpath)>1
    fprintf(2, 'Data sets from different conditions were combined.\n');
end
fpath = fpath{1};