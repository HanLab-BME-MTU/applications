%This script file solves the linear system assembled by 'calFwdOpBF' and
% 'calRHVecBF' for the identification of body force.
fprintf(1,'Solving the linear system : ');

localStartTime = cputime;

coef = zeros(2*dimBF,numTimeSteps);
for jj = 1:numTimeSteps
   [m,n] = size(A{jj});

   coef(:,jj) = (A{jj}.'*A{jj}+sigma*eye(n))\(A{jj}.'*rightU{jj});
end

%Save the coefficient.
coefBFx = zeros(dimFS,numTimeSteps);
coefBFy = zeros(dimFS,numTimeSteps);
for jj = 1:numTimeSteps
   coefBFx(indDomDOF,jj) = coef(1:dimBF,jj);
   coefBFy(indDomDOF,jj) = coef(dimBF+1:2*dimBF,jj);
end
save([resultPath 'coefBFId'],'fs','coefBFx','coefBFy');

fprintf(1,'%f sec.\n',cputime-localStartTime);

