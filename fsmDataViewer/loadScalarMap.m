function I = loadScalarMap(filename)

data = load(filename);
if ~isstruct(data)
    error([filename, ' does not contain any data.']);
end
fieldNames = fieldnames(data);
for iField = 1:numel(fieldNames)
    I = data.(fieldNames{iField});
    
    if (isa(I, 'double') || isa(I, 'single')) && numel(size(I)) == 2 && numel(I) ~= 1
        return;
    end
end

error([filename, ' does not contain any data.']);