%Identify the body force with the constraint on the adhesion location). We 
% solve the nonlinear optimization problem through iterations.

localStartTime = cputime;
fprintf(1,'Solve the constrained force identification problem:\n');
%First get the displacements on the meshing points of 'fs' for the basis
% functions of the body force.
if strcmp(fsAComputed,'yes') == 1 & strcmp(fwdOpComputed,'none') == 0
   fprintf(1,'   Load the matrix fsA.\n');
   load([resultPath 'fsABFId']);
else
   load([modelPath 'solBFId'],'sol');

   fsA = cell(numTimeSteps,1);
   fsMeshPtsIn = fs.mesh.p(:,indDomDOF);
   for jj = 1:numTimeSteps
      fprintf(1,'  Constructing the matrix fsA at time step %d :\n', jj);
      fsA{jj} = zeros(dimBF,2,dimBF,2);

      col = 0;
      ll  = 0;
      for k = 1:dimBF
         ll = ll+1;

         %'bspU1' and 'bspU2' is used to temporarily store the solution at the
         % data points.
         fem.sol = sol{col+1};
         [bspU1 bspU2] = postinterp(fem,'u1','u2',fsMeshPtsIn);
         fsA{jj}(:,1,ll,1) = bspU1.';
         fsA{jj}(:,2,ll,1) = bspU2.';

         fem.sol = sol{col+2};
         [bspU1 bspU2] = postinterp(fem,'u1','u2',fsMeshPtsIn);
         fsA{jj}(:,1,ll,2)  = bspU1.';
         fsA{jj}(:,2,ll,2)  = bspU2.';

         col = col+2;
         if rem(col,100) == 0
            fprintf(1,'    Columns %d out of %d finished.  %f sec.\n', ...
               col,2*dimBF,cputime-localStartTime);
         end
      end

      %Reshape 'fsA';
      fsA{jj} = reshape(fsA{jj},2*dimBF,2*dimBF);
   end
   save([resultPath 'fsAbfId'],'fsA','fsMeshPtsIn');
end

coefBFx  = zeros(dimFS,numTimeSteps);
coefBFy  = zeros(dimFS,numTimeSteps);
coefMCFx = zeros(dimFS,numTimeSteps);
coefMCFy = zeros(dimFS,numTimeSteps);
coefAdh  = zeros(dimFS,numTimeSteps);
for jj = 1:numTimeSteps
   B = zeros(2*numDP(jj),2*dimMCF+dimAdh);
   %First solve the linear system (without constraint) to get the initial
   % displacement on the mesh of 'fs'.
   n       = size(A{jj},2);
   coefA   = (A{jj}.'*A{jj}+sigma*eye(n))\(A{jj}.'*rightU{jj});
   coefBFx(indDomDOF,jj) = coefA(1:dimBF);
   coefBFy(indDomDOF,jj) = coefA(dimBF+1:2*dimBF);
   fsMeshU = fsA{jj}*coefA;

   errDispU   = norm(fsMeshU(inAdhDOF));
   tolDispU   = errDispU*relTolDispU;
   fprintf(1,'   Iterations       Residue    \n');
   k = 1;
   while errDispU > tolDispU & k <= maxIterations
      %Collect the columns for the contraction force. They are not changed.
      B(:,1:2*dimMCF) = A{jj}(:,[inMcfDOF inMcfDOF+dimBF]);

      %Calculate the column for each adhesion basis function.
      B(:,2*dimMCF+1:end) = -A{jj}(:,inAdhDOF)*diag(fsMeshU(inAdhDOF)) - ...
         A{jj}(:,inAdhDOF+dimBF)*diag(fsMeshU(inAdhDOF+dimBF));

      n     = size(B,2);
      coefB = (diag([sigma*ones(1,2*dimMCF) sigmaAdh*ones(1,dimAdh)])+ ...
         B.'*B)\(B.'*rightU{jj});

      %Update the identified displacement.
      coefA([inMcfDOF inMcfDOF+dimBF]) = coefB(1:2*dimMCF);
      coefA(inAdhDOF)       = coefB(2*dimMCF+1:end).*fsMeshU(inAdhDOF);
      coefA(inAdhDOF+dimBF) = coefB(2*dimMCF+1:end).*fsMeshU(inAdhDOF+dimBF);

      newFSMeshU = fsA{jj}*coefA;
      errDispU   = norm(newFSMeshU(inAdhDOF)-fsMeshU(inAdhDOF));
      fsMeshU    = newFSMeshU;

      fprintf(1,'      %d           %f     \n',k,errDispU);
      k = k+1;
   end

   %Coefficient for the contraction force.
   coefMCFx(mcfDOF,jj) = coefB(1:dimMCF);
   coefMCFy(mcfDOF,jj) = coefB(dimMCF+1:2*dimMCF);

   %Coefficient for the adhesion.
   coefAdh(adhDOF,jj) = coefB(2*dimMCF+1:end);
end

%Save the coefficient.
save([resultPath 'coefBFId'],'fs','coefBFx','coefBFy');
save([resultPath 'coefAdhLConstrnBF'],'coefMCFx','coefMCFy','coefAdh');
