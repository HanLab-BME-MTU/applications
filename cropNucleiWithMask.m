function [nucleiCrop,movieInfoCrop]=cropNucleiWithMask(imageFileList,r,nuclei,sheetMask,toDoList)
if nargin < 5 || isempty(toDoList);
    toDoList=1:length(nuclei);
end

for frame=toDoList
    currI=double(imread(imageFileList{frame}));    
    [rowsOrg, colsOrg]=size(currI);
    
    nucleiCrop(frame).img=nuclei(frame).img & sheetMask(frame).mat;
    
    ptsID=find(nucleiCrop(frame).img==1);
    [rpos,cpos]=ind2sub(size(nucleiCrop(frame).img),ptsID);
    nucleiCrop(frame).pos=[cpos rpos];
    
	movieInfoCrop(frame).xCoord(:,1)=nucleiCrop(frame).pos(:,1);
    movieInfoCrop(frame).xCoord(:,2)=r/2;
    movieInfoCrop(frame).yCoord(:,1)=nucleiCrop(frame).pos(:,2);
    movieInfoCrop(frame).yCoord(:,2)=r/2;
    movieInfoCrop(frame).amp(:,1)=currI(ptsID);
    movieInfoCrop(frame).amp(:,2)=0;
end