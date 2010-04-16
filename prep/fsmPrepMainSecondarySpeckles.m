function [IMfinal,candsTot]=fsmPrepMainSecondarySpeckles(I,strg,counter,noiseParam,Speckles,fsmParam)

% fsmPrepMainSecondarySpeckles is the main function of the fsmPrepSecondarySpeckles sub-module
%
% Plots PSFs on the positions of the primary speckles (found by locmax operator) and substructs
% them from the filtered data. Appled again on the resulting image, the locmax-operator finds 
% new (secondary) speckles (intensity and distance significance tests applied)
%
%
% SYNOPSIS   [IMfinal,candsTot]=fmsPrepMainSecondarySpeckles(I,strg,counter,noiseParam,Speckles)
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
% OUTPUT     candsTot   :  augmented cands structure (see fsmPrepTestLocalMaxima.m)
%            IMfinal    :  local maxima map
%
%
%
% DEPENDENCES   fsmPrepMainSecondarySpeckles uses { fsmPrepConfirmSpeckles, fsmPrepSubstructMaxima, 
%                                         fsmPrepCheckDistance, fsmPrepUpdateImax, fsmPrepCheckInfo}
%               fsmPrepMainSecondarySpeckles is used by { fsmPrepMain }
%
% Aaron Ponti, October 4th, 2002

IG=I;
SAVEINFO=1;
if strg == 0 % if you provide all the fields but dont wanna write to disc
    SAVEINFO=0;
end

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

% local minima
Imin=locmin2d(IG,[3,3]);

% intial (filtered) image
[yi,xi,y,x,Imax,candsP,triMin,pMin]=fsmPrepConfirmSpeckles(IG,Imin,noiseParam,userROIbw); % TO DO: update cands

aux=length(candsP);
for i=1:aux
    candsP(i).speckleType=1;
end

candsTot=candsP;

Inew=IG;
candsS=candsP;
HierLevel=2;

while HierLevel<=Speckles(1) && length(candsS)>(Speckles(2)*length(candsTot)) && any([candsS.status])
    
    Inew=fsmPrepSubstructMaxima(Inew,SIG,candsS); % prednite Cands
    [yni,xni,yn,xn,Imax,candsS]=fsmPrepConfirmLoopSpeckles(Inew,noiseParam,triMin,pMin,IG,userROIbw);
    
    aux=length(candsS);
    for i=1:aux
        candsS(i).speckleType=HierLevel; % type flag
    end
    
    candsS=fsmPrepCheckDistance(candsS,candsTot);
    
    HierLevel=HierLevel+1;
    
    if ~isempty(candsS)
        candsTot=cat(2,candsTot,candsS); % concatenating secondary and primary cands structures
    end
end

% remove repetitions because of secondary speckles apearing on the same
% positions as primary (because of floating background)
candsTot=fsmPrepCheckInfo(candsTot);

% obtain updated IM from candstTot
IMfinal=zeros(size(IG));

% Replace loop
validCands = [candsTot(:).status] == 1;
if any(validCands)
    validLmax = vertcat(candsTot(validCands).Lmax);
    validILmax = [candsTot(validCands).ILmax];
    validIdx = sub2ind(size(IG), validLmax(:,1), validLmax(:,2));
    IMfinal(validIdx) = validILmax;
end
% for i=1:length(candsTot)
%     if candsTot(i).status==1
%         IMfinal(candsTot(i).Lmax(1),candsTot(i).Lmax(2))=candsTot(i).ILmax;
%     end
% end

% Save speckle information (cands and locMax) to disk for future use
if SAVEINFO==1
    locMax=IMfinal; %#ok<NASGU>
    cands=candsTot; %#ok<NASGU>
    indxStr=sprintf(strg,counter);
    save(['cands',filesep,'cands',indxStr,'.mat'], 'cands');
    save(['locMax',filesep,'locMax',indxStr,'.mat'], 'locMax');    
end

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
    cands=candsTot; 
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
