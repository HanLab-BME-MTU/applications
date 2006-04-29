function tnsl=trackSpots(mov,sl,region)
%TRACKFSPOTS main function for single/multiple spot tracker
%
% SYNOPSIS  nsl=trackFrame(mov,sl)
%
% INPUT mov : raw microscopy data (2 frames)
%            slC:    list of known spots in current frame
%            slN:   start values in second frame
%
% OUTPUT nsl : tracked frame spots info
%   
% c: 24/1/02	dT

%init params/ model

imgTemp=mov(:,:,:,1,2);
imgMatch=mov(:,:,:,1,1);

 %resize data if necessary (odd size)
idx=~rem(size(imgTemp),2);
imgTemp=imgTemp(1:end-idx(1),1:end-idx(2),1:end-idx(3));
imgMatch=imgMatch(1:end-idx(1),1:end-idx(2),1:end-idx(3));


%mask regions        
matchCoords=cat(1,sl(1).sp.cord);
tempCoords=cat(1,sl(2).sp.cord);

%TODO  make same number of newCoords as origCoords 
trans=tempCoords-matchCoords;
amp=sl(2).sp.amp/sl(1).sp.amp;

%---------------extract masked template coords --------------------
%
[spotsidx mask] = discernspots(tempCoords,size(imgTemp));
idxListTemp=find(mask);
mskTempData=imgTemp(idxList);

% find minimal box for tracking
[ig jg kg]=ind2sub(size(imgTemp),idxList);

%dermine minimal surrounding box
minVTemp=[min(ig) min(jg) min(kg)];
maxVTemp=[max(ig) max(jg) max(kg)];
mskTempBoxSize=maxVTemp-minVTemp+1;

%shift coord
% shiftCTemp=([minVTemp(1)-1, minVTemp(2)-1, minVTemp(3)-1]+ceil(mskTempBoxSize/2));
% tempCoords=tempCoords-ones(size(tempCoords,1),1)*shiftCTemp;

%---------------extract masked match coords --------------------
%

[spotsidx mask] = discernspots(matchCoords,size(imgTemp));
idxListMatch=find(mask);
mskMatchData=imgMatch(idxList);
 
% find minimal box for tracking
[ig jg kg]=ind2sub(size(imgMatch),idxList);

%dermine minimal surrounding box
minVMatch=[min(ig) min(jg) min(kg)];
maxVMatch=[max(ig) max(jg) max(kg)];
mskMatchBoxSize=maxVMatch-minVMatch+1;

%shift coord
% shiftCMatch=([minVMatch(1)-1, minVMatch(2)-1, minVMatch(3)-1]+ceil(mskMatchBoxSize/2));
% matchCoords=matchCoords-ones(size(matchCoords,1),1)*shiftCMatch;



% % if type== mult
% % separate spots: rebuild model with params
% % take % of spot at each pixel
% I_model= multiGaussFit(mskDataSize,parms);
% for i=1:length(sl(1))
%     I_cspot= multiGaussFit(mskDataSize,parms);
%     I_ratio=I_cspot./I_model;
%     spdat(:,:,:,i)=maskmov(:,:,:,:,1).*I_ratio;
% end;

%--------------------main tracking loop-----------------------
%

% init params from fit
parms=[trans amp];

% Set intial parameter change
dpar=[0 ];
%Counter for iteration
iterCt=0;
while(isempty(dpar) | (any(abs(dpar(:,iterCt).*osz_mask)>pm.cmd.prec) & iterCt<pm.cmd.maxIter) )
   iterCt=iterCt+1;
   %Compute new coords according to model
   [ncord, imgModel,modelGrad]=transSpot(mskMatchData,idxListMatch,mskMatchBoxSize,parms);
   % error with current params: e = mA*params-mB
   mB=imgModel-imgTemp;
	mA=[];
   for l=1:size(dIp,1)
      mA=[mA ; dIp(l,:)*modelGrad(:,:,l)];
   end;
	dpar(:,iterCt)= mA\mB(:);

end;

%map spots
