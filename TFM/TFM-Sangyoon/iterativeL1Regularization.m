function mreg=iterativeL1Regularization(G,GTG,d,L,alpha,maxiter,tolx,tolr,m_diff)
% mreg=iterativeL1Regularization(G,GTG,d,L,alpha,maxiter,tolx,tolr,m_diff)
% solves L1 regularization problem with forward matrix G, displacement
% vector d, regularization parameter alpha, semi-norm matrix L.
% tolx and m_diff are used for convergence criteria

% Default for tolr=1.0e-6
if (nargin < 9)
  m_diff=400; % unit: Pa. 
end

if (nargin < 8)
  tolr=1.0e-6;
end

% Default for tolx=1.0e-4;
if (nargin < 7)
  tolx=1.0e-3;
end

% Default for maxiter=100
if (nargin < 6)
  maxiter=100;
end

% unchanging constants in the system that is solved repeatedly
% GTG=G'*G;
GTd=G'*d;
% Start with an initial unweighted solution.
iter=1;
display(['L1 Norm regularization (iteration:' num2str(iter) ')....'])
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

  if (norm(m-mold)/(1+norm(mold)) < tolx) || norm(m-mold)<m_diff
    mreg=m;
    display(['norm(m-mold)=' num2str(norm(m-mold)) ', 1+norm(mold)=' num2str(1+norm(mold)) ', norm(m-mold)/(1+norm(mold))=' ...
      num2str(norm(m-mold)/(1+norm(mold)))])
    return
  end
end

% Give a warning, if desired, but return best solution.
display('L1 norm regularization: maximum iterations exceeded.');
mreg=m;