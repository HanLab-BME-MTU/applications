function I = loadDepolyMap(filename)

load(filename);
if ~exist('depolyMap', 'var')
    error([filename, ' does not contain depolyMap.']);
end
I = -depolyMap;
clear depolyMap;

