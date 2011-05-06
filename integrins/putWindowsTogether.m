function windowsAll = putWindowsTogether(firstWindowFile)

%ask user for window files
if nargin < 1 || isempty(firstWindowFile)
    [fName,dirName] = uigetfile('*.mat','specify first window file in stack');
else
    if iscell(firstWindowFile)
        [fpath,fname,fno,fext] = getFilenameBody(firstWindowFile{1});
        dirName = [fpath,filesep];
        fName = [fname,fno,fext];
    elseif ischar(firstImageFile)
        [fpath,fname,fno,fext] = getFilenameBody(firstWindowFile);
        dirName = [fpath,filesep];
        fName = [fname,fno,fext];
    end
end

%if input is valid ...
if(isa(fName,'char') && isa(dirName,'char'))
    
    %get all file names in stack
    outFileList = getFileStackNames([dirName,fName]);
    numFiles = length(outFileList);
        
else %else, exit
    
    disp('--putWindowsTogether: Bad file selection');
    return
    
end

%put the windows together
for iFile = 1 : numFiles
    
    load(outFileList{iFile});
    windowsAll(iFile,:) = windows;
    
end