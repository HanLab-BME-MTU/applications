function [path1,path2] = getStormPath()

try
    [localPath,orchestraPath] = stormPath();
catch
    disp('Main: The file stormPath.m does not exist!');
    disp('Main: Create it and adapt its content.');
    disp('Main: An example can be found in <storm project directory>/main/stormPathExample.m');
end

if exist(localPath,'file')
    path1 = localPath;
    path2 = orchestraPath;
else
    path1 = orchestraPath;
    path2 = localPath;
end

end