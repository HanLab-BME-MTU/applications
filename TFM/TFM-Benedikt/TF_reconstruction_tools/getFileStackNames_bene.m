% crappy function by benedikt for course
function [listo]=getFileStackNames_bene(firstfilename)

[fpath,fname,ext,shit]=fileparts(firstfilename);
fl = dir(fullfile(fpath,'*.tif'));


if ~isempty(fl)
    iEntry = 1;
    fileList = {};
    for( i = 1:length(fl))
       if(~fl(i).isdir)
          fileList(iEntry) = lower({fl(i).name});
          iEntry = iEntry + 1;
       end;
    end;
    fileList = sort(fileList);
    for i = 1:size(fileList,2)
        listo(i) = {strcat(fpath,filesep,fileList{i})};
    end
else
    listo = [];
    disp('Error!!! No files found');
end


