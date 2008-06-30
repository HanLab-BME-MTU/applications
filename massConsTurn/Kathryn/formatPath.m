function [newPath]=formatPath(oldPath)
% FORMATPATH converts between linux and window paths and vice versa

% assumes pwd contains some overlap with the oldPath
% for example: if pwd = 'S:\scripps\analysis\thompson' (Windows) and
% oldPath =
% '/mnt/fsm/scripps/analysis/thompson/2008-06-26/050224D04Af/images' (Linux),
% formatPath will replace oldPath up through /thompson with pwd, in
% addition to changing the direction of the slash

idxLinux=strfind(oldPath,'/');
idxWindows=strfind(oldPath,'\');


if ispc && ~isempty(idxLinux) % Linux --> Windows
    oldPath=strrep(oldPath, '/', '\');
    
    chunk=pwd;
    pwdSlashIdx=strfind(pwd,'\');
    chunk(1:pwdSlashIdx(end))=[];

    pathIdx=strfind(oldPath,chunk);
    newPath=[pwd oldPath(pathIdx+length(chunk):end)];
    
elseif ~ispc && ~isempty(idxWindows) % Windows --> Linux
    oldPath=strrep(oldPath, '\', '/');
    
    chunk=pwd;
    pwdSlashIdx=strfind(pwd,'/');
    chunk(1:pwdSlashIdx(end))=[];

    pathIdx=strfind(oldPath,chunk);
    newPath=[pwd oldPath(pathIdx+length(chunk):end)];
    
else % same OS, no change needed
    newPath=oldPath;
end

