function []=plotLcurve(M,u_M,forceMesh,fac)
display('This is correct only for meshes with a single basis class!');
display('In the general case one has to substitute: f_rec=sol_coef * vol');

if nargin < 4 || isempty(fac)
    fac=10;
end

lowerLim=-12;
upperLim=-2;
numPoints=fac*(upperLim-lowerLim)+1;
L=logspace(lowerLim,upperLim,numPoints);

residuals_u=zeros(numPoints,1);
norm_f=zeros(numPoints,1);

Msq  = M'*M;
Mu   = M'*u_M;
Imat = eye(2*forceMesh.numBasis);

for j=1:length(L)
    %sol_coef=(L(j)*eye(2*forceMesh.numBasis)+M'*M)\(M'*u_M);
    sol_coef=(L(j)*Imat+Msq)\Mu;
    u_rec=M*sol_coef;
    f_rec=sol_coef;

    residuals_u(j)=sum((u_M-u_rec).^2);
    norm_f(j)=sum((f_rec).^2);
end
figure(100)
%loglog(residuals_u,norm_f)
plot(log(residuals_u),log(norm_f))
hold on
plot(log(residuals_u(1:fac:end)),log(norm_f(1:fac:end)),'.k')
for k=1:fac:length(L)
    text(log(residuals_u(k)),log(norm_f(k)),num2str(L(k)))
end
xlabel('log|u_{data}-u_{model}|^2')
ylabel('log|f|^2')
hold off
title('L-curve')