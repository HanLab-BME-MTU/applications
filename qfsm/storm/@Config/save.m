function save(obj,fullPath)
% Save the object in the .mat format
if exist(fullPath,'file')
    disp('Config: WARNING! Configuration file overwritten!');
    msgbox('Configuration file overwritten!','Warning','warn');
    copyfile(fullPath,[fullPath '.bak']);
end
save(fullPath,'obj');
disp('Config: Config object written to file!');
end