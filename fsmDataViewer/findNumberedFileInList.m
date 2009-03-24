function [file found] = findNumberedFileInList(fileList, no)

found = 0;
file = '';

% Check as a first guess fileList{no}
if no > 0 && no <= numel(fileList)
    [dummy, body, no2] = getFilenameBody(fileList{no});
    
    if str2double(no2) == no
        file = fileList{no};
        found = 1;
        return;
    end
end

% Search into the whole list
for iFile = 1:numel(fileList)
    [dummy, body, no2] = getFilenameBody(fileList{iFile});
    
    if str2double(no2) == no
        file = fileList{iFile};
        found = 1;
        return;
    end
end
end