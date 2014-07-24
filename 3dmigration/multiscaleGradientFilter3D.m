function [maxResp,dX,dY,dZ,maxRespScale] = multiscaleGradientFilter3D(imageIn)

sMin = 1;
sMax = 3;
nSig = 10;

sigmas = linspace(sMin,sMax,nSig);

maxResp = zeros(size(imageIn));
maxRespScale = zeros(size(imageIn));
dX = zeros(size(imageIn));
dY = zeros(size(imageIn));
dZ = zeros(size(imageIn));

for j = 1:nSig
            
    [dXtmp,dYtmp,dZtmp] = gradientFilterGauss3D(imageIn,sigmas(j));            
    
    sMag = sqrt(dXtmp.^2 + dYtmp .^2 + dZtmp.^2);
    
    isBetter = sMag > maxResp;    
    maxRespScale(isBetter) = sigmas(j); 
    maxResp(isBetter) = sMag(isBetter);
    dX(isBetter) = dXtmp(isBetter);
    dY(isBetter) = dYtmp(isBetter);
    dZ(isBetter) = dZtmp(isBetter);
                
end