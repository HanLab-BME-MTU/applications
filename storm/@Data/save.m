function save(obj,fullPath)
% Save the object in the .mat format
save(fullPath,'obj');
disp('Data: Data object written to file!');
end
