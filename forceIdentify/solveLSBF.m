%This script file solves the linear system assembled by 'calFwdOpBF' and
% 'calRHVecBF' for the identification of body force.
fprintf(1,'Solving the linear system : ');

localStartTime = cputime;

%Reshape 'A';
A = reshape(A,2*numDP,2*dimBF);
[m,n] = size(A);

%Regularization parameter for 'solveLS'.
sigma = 1e6;

coef  = (A.'*A+sigma*eye(n))\(A.'*rightU);

fprintf(1,'%f sec.\n',cputime-localStartTime);

