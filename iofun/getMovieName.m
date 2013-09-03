function str = getMovieName(data)
if isempty(data.date)
    str = [' ' getCellDir(data)];
else
    str = [' ' num2str(data.date) filesep getCellDir(data)];
end