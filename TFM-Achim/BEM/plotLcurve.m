function []=plotLcurve(M,sol_mats,u_M,forceMesh,fac)
% INPUT:
% M         : M is either the fwd Map, or a force field contain all
%             information needed for the remaining inputs except fac.
%             The L-curve will be plot ONLY FOR THE FIRST FRAME. Execute
%             for e.g. plotLcurve(forceField(5)) to plot it for the 5th
%             frame.
% sol_mats  : structure that contains the matrices that are need to
%             reconstruct the solutions. sol_mats has to be in the form
%             that is understood by the function calcSolFromSolMatsFastBEM.
% u_M       : displacement information for which M has been build up.
% forceMesh : contains the force mesh for which M has been calculated.
% fac       : the refinement factor. The higher fac the finer the plot will
%             be.

% display('This is correct only for meshes with a single basis class!');
% display('In the general case one has to substitute: f_rec=sol_coef * vol');

if isstruct(M)
    forceField=M(1);
    M         = forceField.par.M;
    u_M       = forceField.par.u;
    sol_mats  = forceField.par.sol_mats;
    forceMesh = forceField.par.forceMesh;
end

if nargin < 4 || isempty(fac)
    fac=10;
end

lowerLim=-15;
upperLim=0;
numPoints=fac*(upperLim-lowerLim)+1;
L=logspace(lowerLim,upperLim,numPoints);
L=[0 L Inf];

residuals_u=zeros(numPoints,1);
norm_f=zeros(numPoints,1);

% Msq  = M'*M;
% Mu   = M'*u_M;
% Imat = eye(2*forceMesh.numBasis);

for j=1:length(L)
    textMess=['Calculate: ',num2str(numel(L)),' values for regularization parameter'];
    progressText(j/numel(L),textMess);
    
    regParam=L(j);
    %sol_coef=(L(j)*eye(2*forceMesh.numBasis)+M'*M)\(M'*u_M);
    %sol_coef=(L(j)*Imat+Msq)\Mu;
    [~,~,sol_coef]=calcSolFromSolMatsFastBEM(M,sol_mats,u_M,forceMesh,regParam,[],[]);
    u_rec=M*sol_coef;
    
    % basis functions with a larger support have to be weighted more!
    [normWeights]=getNormWeights(forceMesh);
    eyeWeights =diag(normWeights);    
    f_rec=eyeWeights*sol_coef;

    residuals_u(j)=sum((u_M-u_rec).^2);
    norm_f(j)=sum((f_rec).^2);
end
figure(100)
%loglog(residuals_u,norm_f)
plot(log(residuals_u),log(norm_f),'r')
hold on
plot(log(residuals_u(2:fac:end-1)),log(norm_f(2:fac:end-1)),'.k')
for k=2:fac:length(L)-1
    text(log(residuals_u(k)),log(norm_f(k)),num2str(L(k)))
end
% plot extremal values:
plot(log(residuals_u(1)),log(norm_f(1)),'xr')
text(log(residuals_u(1)),log(norm_f(1)),num2str(L(1)))
plot(log([residuals_u(end) residuals_u(end)]),[min(log(norm_f(1:end-1))) max(log(norm_f(1:end-1)))],'--r')
text(log(residuals_u(end)),mean(log(norm_f(1:end-1))),num2str(L(end)))
xlabel('log|u_{data}-u_{model}|^2')
ylabel('log|f|^2')
hold off
title('L-curve')