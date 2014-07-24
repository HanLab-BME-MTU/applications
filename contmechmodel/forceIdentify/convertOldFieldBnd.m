%oldFieldBndFile = [mechDir filesep 'fieldGeom.mat'];

%if exist(oldFieldBndFile,'file') == 2
%   load(oldFieldBndFile);
%   fieldBnd.x = fieldPGx;
%   fieldBnd.y = fieldPGy;
%else
%   fprintf(1,'There is no old field boundary defined.');
%   return;
%end

oldRecBndFile = [mechDir filesep 'recGeom.mat'];
if exist(oldRecBndFile,'file') == 2
   load(oldRecBndFile);
   fieldBnd.x  = recPGx;
   fieldBnd.y  = recPGy;
   fieldBnd.VI = recPGVI;
else
   fprintf(1,'There is no old recover boundary defined.');
   return;
end

%Create the directory for storing field boundary if it does not exist.
fieldBndDir = [mechDir filesep 'fieldBnd'];
if ~isdir(fieldBndDir)
   success = mkdir(mechDir,'fieldBnd');
   if ~success
      error('Trouble making directory for field boundary.');
   end
end

fieldBndFile = [fieldBndDir filesep 'fieldBnd.mat'];
save(fieldBndFile,'fieldBnd');
