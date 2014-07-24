function cfg = load(fullPath)
% Load the object
tmp = load(fullPath,'-mat');
cfg = tmp.obj;
disp('Config: Config object read from file!');
end

