function [cands2] = candsToSubpixelN(img,cands,sigma)
% convert integer speckle position values in cands to subpixel accuracy
% using mixture-model fitting
% 
% INPUT:
% img: normalized image
% cands: cands determined by qfsm
% sigma: Gaussian sigma, determined by PSF
% bitdepth: necessary for image scaling, e.g. 12, 16
% 
% OUTPUT:
% cands2: new cands with double-precision pixel positions
%
% DEPENDENCiES   candsToSubpixelN uses fitMixModel
%               candsToSubpixelN is used by fsmPrepMainSecondarySpeckles
%
% Dinah Loerke, Mar 1, 2005

maxbitdepth=2^16;
%pointer to significant speckles and fitting inputs
signi=find([cands.status]==1);
xycand_all=cat(1,cands.Lmax);
xycand_sig=xycand_all(signi,1:2);
ampcand_all=cat(1,cands.deltaI);
ampcand_sig=ampcand_all(signi);
bgcand_all=cat(1,cands.IBkg);
bgcand_sig=bgcand_all(signi);

%fitting input matrix
xycand_fit=[xycand_sig];
b_fit=mean(bgcand_sig)*maxbitdepth;
amp_fit=ampcand_sig*maxbitdepth;

% fit w/ Jacobian with previous information about all parameters
% (xy, bg and amplitudes) taken from 
% instead of normalizing the image, rescaled the intensity (amp_fit), this has
% better effects in the fitting
img=img*maxbitdepth;

%Mixture model fitting
[estimates] = fitMixModel(img,xycand_fit,sigma, amp_fit, b_fit); 
amps_new=estimates(:,3);
xy_new=estimates(:,1:2);

%normalize resulting amplitudes back to cand format with bitdethp
% and rewrite cands
cands2=cands;
for i=1:length(signi)
    cands2(signi(i)).Lmax=xy_new(i,1:2);
    deltaInew=amps_new/maxbitdepth;
    cands2(signi(i)).deltaI=deltaInew(i);
end

end %of funciton
        