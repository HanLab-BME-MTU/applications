function [coveredDist,coveredArea,toDoList,badFlag]=calcCoveredArea(sheetMask,dPix,toDoList,doCheck)
if nargin < 3 || isempty(toDoList);
    toDoList=1:length(sheetMask);
end
badFlag=0;

dPixMax=max(dPix);
% First check if there is a really bad data set.
if dPixMax>100
    firstBadFrame=find(dPix==dPixMax,1,'first');
    display(['ToDoList cut off because of bad frame (dPix>100): ',num2str(firstBadFrame)]);
    toDoList=1:firstBadFrame-1;
    badFlag=1;
    dPixMax=max(dPix(toDoList));
end

[rowsOrg,colsOrg]=size(sheetMask(1).mat);
for frame=toDoList(1:end)
    % first crop the masks according to the largest dPix
    sheetMaskCrop(frame).mat =sheetMask(frame).mat((1+dPixMax):(rowsOrg-dPixMax),(1+dPixMax):(colsOrg-dPixMax));   
    diffMask   =sheetMaskCrop(frame).mat-sheetMaskCrop(1).mat;
       
    coveredArea(frame)=sum(diffMask(:));    
    coveredDist(frame)=coveredArea(frame)/rowsOrg;
end

if nargin<4 || doCheck==1
    % statistical variance cutoff numSigma
    numSigma=10;
    % minimal allowed offset of variance (to not penalize slowly moving sheets)
    minVar=10;

    dA=coveredDist(2:end)-coveredDist(1:end-1); 

    minCutOff=min(median(dA)-minVar,median(dA)-numSigma*mad(dA));
    maxCutOff=max(median(dA)+minVar,median(dA)+numSigma*mad(dA));

    checkVec= ((dA>maxCutOff) | (dA<minCutOff));
    cutOff=find(checkVec,1);

    if ~isempty(cutOff)
        display(['ToDoList cut off because sheet propagation is suspicious at frame: ',num2str(cutOff)]);
        toDoList=1:cutOff;
        badFlag=1;
    end
    coveredDist=coveredDist(toDoList);
    coveredArea=coveredArea(toDoList);
end
