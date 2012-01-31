function cellDir = getCellDir(data)

chParents = cellfun(@(c) getParentDir(c), data.channels, 'UniformOutput', false);
sParent = getParentDir(data.source);

v = strcmp(sParent, chParents);

if numel(data.channels)==1 || ~all(v)
    cellDir = data.source;
else % all channels at same level, below 'source'
    cellDir = sParent;
end
