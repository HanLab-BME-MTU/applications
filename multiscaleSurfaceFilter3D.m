function [maxResp,d2X,d2Y,d2Z,maxRespScale] = multiscaleSurfaceFilter3D(imageIn)

sMin = 1;
sMax = 5;%The high scales are very rarely helpful - decrease this!!
%nSig = 10;

%sigmas = logspace(log10(sMin),log10(sMax),nSig);
sigmasXY = [1.5 2 4 ];%TTTEEEEMMMPPPP!!!
sigmasZ = sigmasXY*4/3;

nSig = numel(sigmasXY);

maxResp = zeros(size(imageIn));
maxRespScale = zeros(size(imageIn));
d2X = zeros(size(imageIn));
d2Y = zeros(size(imageIn));
d2Z = zeros(size(imageIn));

for j = 1:nSig
            
    [d2Xtmp,d2Ytmp,d2Ztmp] = surfaceFilterGauss3D(imageIn,[sigmasXY(j) sigmasXY(j) sigmasZ(j)]);            
    
    d2Xtmp(d2Xtmp<0) = 0;
    d2Ytmp(d2Ytmp<0) = 0;
    d2Ztmp(d2Ztmp<0) = 0;        
        
    %Get magnitude and normalize response based on sigma to give comparable
    %responses at different scales
    %sMag = sqrt(d2Xtmp.^2 + d2Ytmp .^2 + (d2Ztmp * sigmasZ(j)/sigmasXY(j)).^2) * sigmasXY(j);%TEMP - IS THIS RIGHT????? 
    sMag = sqrt((d2Xtmp * sigmasXY(j)) .^2 + ...
                (d2Ytmp * sigmasXY(j)) .^2 + ...
                (d2Ztmp * sigmasZ(j)*4)  .^2);
    
    isBetter = sMag > maxResp;    
    maxRespScale(isBetter) = sigmasXY(j); 
    maxResp(isBetter) = sMag(isBetter);
    d2X(isBetter) = d2Xtmp(isBetter);
    d2Y(isBetter) = d2Ytmp(isBetter);
    d2Z(isBetter) = d2Ztmp(isBetter);
                
end