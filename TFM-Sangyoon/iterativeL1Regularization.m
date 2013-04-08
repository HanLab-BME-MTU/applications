function mreg=iterativeL1Regularization(G,GTG,d,L,alpha,maxiter,tolx,tolr)
% adopted from Aster et al.
% Default for tolr=1.0e-6
if (nargin < 7)
  tolr=1.0e-6;
end

% Default for tolx=1.0e-4;
if (nargin < 6)
  tolx=1.0e-4;
end

% Default for maxiter=100
if (nargin < 5)
  maxiter=100;
end

% unchanging constants in the system that is solved repeatedly
% GTG=G'*G;
GTd=G'*d;
% Start with an initial unweighted solution.
iter=1;
display(['L1 Norm regularization (iteration:' num2str(iter) ')'])
tic
m=(2*GTG+alpha*(L'*L))\(2*G'*d);
toc
% iterate until maxiter or we converge and return
while (iter < maxiter)
  iter=iter+1;

  % get get the magnitude of Lm, but don't let any element be less than tolr
  absLm=abs(L*m);
  absLm(absLm<tolr)=tolr;

  % build the diagonal weighting matrix for this iteration
  R=diag(1./absLm);

  mold=m;
  
%   % use LSQR to get m that minimizes argmin||[            M              ]f - [d]||2
%   %                                           sqrt(alpha/2)*sqrt(absLm)*L      0
%   display('LSQR...')
%   tic
%   A_bottom = sqrt(alpha/2)*sqrt(R)*L;
%   A = [G;A_bottom];
%   b = [d;zeros(size(A_bottom,1),1)];
%   m = lsqr(A,b,tolx,maxiter);
%   toc
%   % get the new iterate and check for convergance
%   display('QR...')
%   tic
%   [Q,Rr] = qr(2*GTG+alpha*L'*R*L);
%   m=Rr\(Q'*(2*GTd));
%   toc
%   display(norm(m))
  display(['(iteration:' num2str(iter) ')'])
  tic
  m=(2*GTG+alpha*L'*R*L)\(2*GTd);
  toc
  display(['norm(m-mold)=' num2str(norm(m-mold)) ', 1+norm(mold)=' num2str(1+norm(mold)) ', norm(m-mold)/(1+norm(mold))=' ...
      num2str(norm(m-mold)/(1+norm(mold)))])

  if (norm(m-mold)/(1+norm(mold)) < tolx)
    mreg=m;
    return
  end
end

% Give a warning, if desired, but return best solution.
display('L1 norm regularization: maximum iterations exceeded.');
mreg=m;