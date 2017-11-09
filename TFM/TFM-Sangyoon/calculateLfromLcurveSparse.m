function [sol_coef,reg_corner] = calculateLfromLcurveSparse(L,M,MpM,u,eyeWeights,maxIter,tolx,tolr,LcurveDataPath,LcurveFigPath,LcurveFactor)
%examine a logarithmically spaced range of regularization parameters
alphas=10.^(log10(L)-2.5:1.25/LcurveFactor:log10(L)+2);
rho=zeros(length(alphas),1);
eta=zeros(length(alphas),1);
% eta0=zeros(length(alphas),1);
msparse=zeros(size(M,2),length(alphas));
% if matlabpool('size')==0
%     matlabpool open
% end
tolFactor = 20; % make the L-curve calculation faster with generous tolx with this factor
for i=1:length(alphas);
    disp(['testing L = ' num2str(alphas(i)) '... '])
    msparse(:,i)=iterativeL1Regularization(M,MpM,u,eyeWeights,alphas(i),maxIter,tolx*tolFactor,tolr);
    rho(i)=norm(M*msparse(:,i)-u);
    eta(i)=norm(msparse(:,i),1);
%     eta0(i)=sum(abs(msparse(:,i))>1);
end

% Find the L-corner
% [reg_corner,ireg_corner,~]=l_curve_corner(rho,eta,alphas);
save(LcurveDataPath,'rho','eta','alphas','L','msparse','-v7.3'); % saving before selection.
[reg_corner,ireg_corner,~]=regParamSelecetionLcurve(rho,eta,alphas,L);

% Also, I can use L0 norm information to choose regularization parameter

% Plot the sparse deconvolution L-curve.
hLcurve = figure;
set(hLcurve, 'Position', [100 100 500 500])

loglog(rho,eta,'k-');
ylim([min(eta) max(eta)])
xlabel('Residual Norm ||Gm-d||_{2}');
ylabel('Solution Norm ||m||_{1}');
hold on
% mark and label the corner
if mod(ireg_corner,1)>0 % if ireg_corner is interpolated
    rho_corner = rho(floor(ireg_corner))+mod(ireg_corner,1)*(rho(floor(ireg_corner)+1)-rho(floor(ireg_corner)));
    eta_corner = eta(floor(ireg_corner))+mod(ireg_corner,1)*(eta(floor(ireg_corner)+1)-eta(floor(ireg_corner)));
else
    rho_corner = rho(ireg_corner);
    eta_corner = eta(ireg_corner);
end    
H=loglog(rho_corner,eta_corner,'ro');
set(H,'markersize',6)
H=text(rho_corner,1.1*eta_corner,...
    ['    ',num2str(reg_corner,'%5.1e')]);
set(H,'Fontsize',7);
% axis([1e-2 100 0.001 1e8])
disp('Displaying the 1-norm L-curve')
% print -deps2 nameSave
% print(hLcurve,strcat(nameSave,'.eps'),'-depsc')
saveas(hLcurve,LcurveFigPath);
save(LcurveDataPath,'rho','eta','reg_corner','ireg_corner','alphas','rho_corner','eta_corner','msparse','-v7.3');

if mod(ireg_corner,1)>0 % if ireg_corner is interpolated
    disp(['L-corner regularization parmater L = ' num2str(reg_corner) '... final solution calculation ...'])
    sol_coef=iterativeL1Regularization(M,MpM,u,eyeWeights,reg_corner,maxIter,tolx,tolr);
else
    sol_coef = msparse(:,ireg_corner);
end