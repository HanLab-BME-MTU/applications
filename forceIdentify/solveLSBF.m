%This script file solves the linear system assembled by 'calFwdOpBF' and
% 'calRHVecBF' for the identification of body force.
fprintf(1,'Solving the linear system : ');

localStartTime = cputime;

%Regularization parameter for 'solveLS'.
sigma = 1e2; %1e6;

coef = zeros(2*dimBF,numTimeSteps);
for jj = 1:numTimeSteps
   [m,n] = size(A{jj});

   coef(:,jj) = (A{jj}.'*A{jj}+sigma*eye(n))\(A{jj}.'*rightU{jj});
end

fprintf(1,'%f sec.\n',cputime-localStartTime);

