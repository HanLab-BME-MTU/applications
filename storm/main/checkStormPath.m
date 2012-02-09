function pathOut = checkStormPath(pathIn)

[localPath,remotePath] = getStormPath();
pathOut = strrep(pathIn,remotePath,localPath);

% Convert windows file paths to linux file paths
if isunix
    pathOut = strrep(pathOut,'\','/');
else
    pathOut = strrep(pathOut,'/','\');
end

end


