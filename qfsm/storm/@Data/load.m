function data = load(fullPath)
% Load the object
tmp = load(fullPath,'-mat');
data = tmp.obj;
disp('Data: Data object read from file!');
end

