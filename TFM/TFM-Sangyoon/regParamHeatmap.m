function []=regParamHeatmap(alphas,msparse,forceMesh,outputPath,fmin,fmax,quiverTrue)
%regParamHeatmap(rho,eta,alpha,msparse,outputPath) creates heatmaps for
%each alpha and store them in outputPath
% rho: residual norm
% eta: semi norm
% alphas: regularization parameters
if nargin < 7
    quiverTrue=false;
end
x_out=forceMesh.p(:,1);
y_out=forceMesh.p(:,2);
formatSpec = '%10.2e\n';
for ii=1:length(alphas);
    disp(['Drawing L = ' num2str(alphas(ii)) '... '])
    sol_coef = msparse(:,ii);
    [fx,fy,x_out,y_out]=calcForcesFromCoef(forceMesh,sol_coef,x_out,y_out,'new');
    h = generateHeatmapFromGridData(x_out,y_out,fx,fy,[outputPath filesep 'Map' num2str(alphas(ii),formatSpec)],0,fmin,fmax,quiverTrue);
    close(h)
end
