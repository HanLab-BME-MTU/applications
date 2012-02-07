function save(obj,fullPath)
% Save the object in the .mat format
save(fullPath,'obj');
disp('Timing: Timing object written to file!');
end
