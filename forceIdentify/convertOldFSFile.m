oldFSFile = [reslDir filesep 'coefBFId.mat'];

if ~exist(oldFSFile,'file')
   fprintf(1,'Old ''fs'' file does not exist.\n');
   return;
end

s = load(oldFSFile);
fs.fem = s.fs;

%Find the indices of the DOFs whose shape functions have support disconnected 
% from the boundary.
N = assemble(fs.fem,'out','N');
indDomDOF = find(full(sum(N,1))==0);

%The Degree of Freedom vector for the basis or shape function in the finite
% element space.
%Dimention of the function space where the basis can be nonzero on the 
% boundary.
dimFS  = flngdof(fs.fem);

%Dimention of the function space where the basis function is kept zero on the
% boundary.
dimBF = length(indDomDOF);

fs.dimFS     = dimFS;
fs.dimBF     = dimBF;
fs.indDomDOF = indDomDOF;

if strcmp(isFieldBndFixed,'yes')
   fsFile       = [femSolBasisBFDir filesep 'fs.mat'];
else
   DTDir = ['DT' sprintf(DTIndexForm,jj)]; 
   if ~isdir([femSolBasisBFDir filesep DTDir])
      [success,msg,msgId] = mkdir(femSolBasisBFDir,DTDir);
      if ~success
         error('Trouble making directory for ''femSolBasisBF''.');
      end
   end
   fsFile = [femSolBasisBFDir filesep DTDir filesep 'fs.mat'];
end

save(fsFile,'fs');
