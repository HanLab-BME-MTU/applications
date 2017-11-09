function saveWithHeaders(folder, fileName, headers, variable)

% FileName has to be a propper path
% headers should be a cell will text entries {'header1', header2',...}
% variable is the matrix to save
% there in no consistency check in terms of numeber of columns

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


tabbedHearders=[];
for it=1: size(headers, 2)
    thisHeader=[headers{1, it} ' \t '];
    tabbedHearders=[tabbedHearders thisHeader];
end
tabbedHearders=[tabbedHearders '\r' ' \n'];

fileName=[datestr(now, 'yy.mm.dd.HH.MM') fileName];
fileHandle=fopen([folder fileName], 'w');
fprintf(fileHandle, tabbedHearders);
fclose(fileHandle);
save([folder fileName], 'variable', '-ASCII', '-APPEND', '-DOUBLE');