oldFwdMapBFFile = [reslDir filesep 'AbfId.mat'];

if ~exist(oldFwdMapBFFile,'file')
   fprintf(1,'Old fwdMapBFFile does not exist.\n');
   return;
end

s = load(oldFwdMapBFFile);
A = s.A{1};

fwdMapBFFile = [fwdMapBFDir filesep 'A' sprintf(DTIndexForm,jj) '.mat'];
save(fwdMapBFFile,'A');
