%This script file solves the linear system assembled by 'calFwdOpTF' and
% 'calRHVecTF' for the identification of the boundary traction force.
fprintf(1,'Solving the linear system: ');

localStartTime = cputime;

%Regularization parameter for 'solveLS'.
sigma = 1e2;

coef = cell(numEdges,1);
for k = 1:numEdges
   %Reshape 'A';
   A{k} = reshape(A{k},2*numDP,2*dimTF{k});
   [m,n] = size(A{k});

   coef{k} = (A{k}.'*A{k}+sigma*eye(n))\(A{k}.'*rightU{k});
end

fprintf(1,'%f sec.\n',cputime-localStartTime);

