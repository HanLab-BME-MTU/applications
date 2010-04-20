function cands=fsmPrepMainSecondarySpeckles(I,strg,counter,noiseParam,Speckles,fsmParam)

% fsmPrepMainSecondarySpeckles is the main function of the fsmPrepSecondarySpeckles sub-module
%
% Plots PSFs on the positions of the primary speckles (found by locmax operator) and substructs
% them from the filtered data. Appled again on the resulting image, the locmax-operator finds 
% new (secondary) speckles (intensity and distance significance tests applied)
%
%
% SYNOPSIS   cands=fmsPrepMainSecondarySpeckles(I,strg,counter,noiseParam,Speckles)
%
% INPUT      I          :  filtered image
%            strg       :  format string for the correct file numbering
%            counter    :  image number
%            noiseParam :  noise parameters for statistical speckle selection
%            Speckles   :  (1) contains information for the hierarchical level 
%                          (2) minimal increase in (%) of new speckles
%                          (before stopping)
%            fsmParam   :  (optional) fsmParam structure for SpeckTackle
%
% OUTPUT     cands      :  augmented cands structure (see fsmPrepTestLocalMaxima.m)
%
%
%
% DEPENDENCES   fsmPrepMainSecondarySpeckles uses { fsmPrepConfirmSpeckles, fsmPrepSubstructMaxima, 
%                                         fsmPrepCheckDistance, fsmPrepUpdateImax, fsmPrepCheckInfo}
%               fsmPrepMainSecondarySpeckles is used by { fsmPrepMain }
%
% Aaron Ponti, October 4th, 2002

IG=I;

if nargin==5
    fsmParam=[];
    userROIbw=[];
end

if ~isempty(fsmParam)
    if fsmParam.prep.drawROI~=0 % Either 1 (drawn) or 2 (loaded)
        
        % Load user-defined ROI from disk
        ROIname=[fsmParam.main.path,filesep,'userROI.mat'];
        if exist(ROIname, 'file')==2 % File found
            tmp=load(ROIname);
            userROIbw=tmp.userROIbw;
            clear tmp;
        else
            userROIbw=[];
        end
    else
        userROIbw=[];    
    end
end

SIG=1.88; % for the twice convolved image (or 1.77)

% SB: We need to add virtual points from the cell edge to Imin with special
% value (i.e. -1000). These points fall ouside the cell footprint, i.e.
% I(p) == 0. If no cell mask is provided, there is no additional point.

% local minima
Imin=locmin2d(IG,[3,3]);

% intial (filtered) image
[cands,triMin,pMin] = fsmPrepConfirmSpeckles(IG,Imin,noiseParam,userROIbw);

aux=length(cands);
for i=1:aux
    cands(i).speckleType=1;
end

Inew=IG;
candsS=cands;
HierLevel=2;

while HierLevel<=Speckles(1) && length(candsS)>(Speckles(2)*length(cands)) && any([candsS.status])
    
    Inew=fsmPrepSubstructMaxima(Inew,SIG,candsS); % prednite Cands
    candsS = fsmPrepConfirmLoopSpeckles(Inew,noiseParam,triMin,pMin,IG,userROIbw);
    
    aux=length(candsS);
    for i=1:aux
        candsS(i).speckleType=HierLevel; % type flag
    end
    
    candsS=fsmPrepCheckDistance(candsS,cands);
    
    HierLevel=HierLevel+1;
    
    if ~isempty(candsS)
        cands=cat(2,cands,candsS); % concatenating secondary and primary cands structures
    end
end

% remove repetitions because of secondary speckles apearing on the same
% positions as primary (because of floating background)
cands=fsmPrepCheckInfo(cands);

% obtain updated IM from candstTot
locMax=zeros(size(IG));

validCands = [cands(:).status] == 1;
if any(validCands)
    validLmax = vertcat(cands(validCands).Lmax);
    validILmax = [cands(validCands).ILmax];
    validIdx = sub2ind(size(IG), validLmax(:,1), validLmax(:,2));
    locMax(validIdx) = validILmax; %#ok<NASGU>
end

% Save speckle information (cands and locMax) to disk for future use
indxStr=sprintf(strg,counter);
save(['cands',filesep,'cands',indxStr,'.mat'], 'cands');
save(['locMax',filesep,'locMax',indxStr,'.mat'], 'locMax');

%-----------------------------------------
% Estimate subpixel positions of speckles
%-----------------------------------------
% The size of the GaussKernel used for filtering the image (in
% fsmPrepareImage) has, as of March 2005, been set to the real value of 
% psfsigma as determined by the physical parameters,
% psfsigma=0.21*(lambda/NA)/pixelsize
% Thus, the the mixture model fitting may be performed without loss of
% information on the filtered image (called I), rather than the original
% image (called oI)
%
% However, it is important to note that performing the Gauss fit in the
% mixture model on the filtered image (which is broadened), rather than on
% the original one, requires to modify the sigma of the mixture-model fit!

if isstruct(fsmParam) && fsmParam.prep.subpixel==1 
    psfsigma     = fsmParam.prep.psfSigma;       % true physical sigma of the image point-spread function, caluclated by sigma=0.21*(lambda/NA)/pixelsize
    filtersigma  = fsmParam.prep.filterSigma;    % sigma used for the low-pass filtering; except where specifically
                                                % stated differently by the user, filtersigma should have the same value as psfsigma; 
                                                % for filtersigma>psfsigma, image information is lost during filtering!!                                            % same value as 
    %mixture model Gauss sigma (mmsigma) is calculated from psfsigma and
    %filtersigma; in the usual case where psfsigma=filtersigma, then
    %mmsigma=sqrt(2)*psfsigma
    mmsigma=sqrt(psfsigma^2+filtersigma^2);
    bitDepth = log2(fsmParam.main.normMax + 1);
    disp(['psfsigma=',num2str(psfsigma),'   filtersigma=',num2str(filtersigma),'   mixmodsigma=',num2str(mmsigma)]);
    disp('calculating sub-pixel locations...');
    candsSP = candsToSubpixelN(I, cands, mmsigma, bitDepth); %#ok<NASGU>
    eval( (strcat('save cands',filesep,'cands',indxStr,'_spa.mat candsSP;')) );
end
