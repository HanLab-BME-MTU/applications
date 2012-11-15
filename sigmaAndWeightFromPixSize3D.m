function [sigXY,sigZ,wtZ] = sigmaAndWeightFromPixSize3D(pixXY,pixZ,NA,n,lambda)
%SIGMAANDWEIGHTFROMPIXSIZE3D heuristic selection of surface filter sigmas and weighting from pixel size
%
% [sigXY,sigZ,wtZ] = sigmaAndWeightFromPixSize3D(pixXY,pixZ)
%
% Gives XY and Z multi-scale surface filter sigmas based on the input
% imaging parameters. Minimum sigmas are determined by approximate PSF,
% weighting and larger sigmas are determined by stupid trial and error
% method. Maybe someday I'll have time to base this on simulated images....
%
% Hunter Elliott
% 10/2012


%Hard-coded parameters
minSig = 1;%Minimum sigma to use. Numerical form of kernel gets shitty below this size.
nSig = 3;%Number of sigmas to return
sigRat = 1;%Ratio of filter sigmas to PSF sigmas;
zFact = 1;%Arbitrary increase of Z sigmas because it works...

%Get PSF approximation sigma in XY
[sPSFxy,sPSFz] = calcFilterParms(lambda,NA,n,[],[],[pixXY, pixZ]);

%Calc sigmas to use. Don't use sigma smaller than minimum size.
sigXY1 = max(sPSFxy*sigRat,minSig);
sigXY = sigXY1 .* 2.^(0:1:nSig-1);

sigZ1 = max(sPSFz*sigRat,minSig);
sigZ = sigZ1 .* 2.^(0:1:nSig-1) * zFact;

%wtZ = (sigZ(1)/sigXY(1))^2;
wtZ = (sPSFz(1)/sPSFxy(1))^2;%This should give approximately normalized response in the different dimensions.
%wtZ = 1;
%wtZ = (sPSFxy(1)/sPSFz(1));



