function cellDir = getCellDir(data)

chParents = cellfun(@(c) getParentDir(c), data.channels, 'UniformOutput', false);
sParent = getParentDir(data.source);

v = strcmp(sParent, chParents);

if all(v) % all channels at same level
    cellDir = sParent;
elseif sum(v)==1
    cellDir = data.source;
else
    error('Incompatible data structure.');
end