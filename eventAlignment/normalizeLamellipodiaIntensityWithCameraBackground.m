function y1Norm = normalizeLamellipodiaIntensityWithCameraBackground(y1,samplesAvg,cameraBackground)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[n m p]=size(samplesAvg);

nWin=n;
nLayer=m;
nTime=p;


for j=1:nTime

    offset=mean(cameraBackground);
    
    %for l=1:nWin
    %    minVal(l)=min(samplesAvg(l,3:4,j));
    %end
    %laMean=nanmean(minVal);
    
    laMean=nanmean(samplesAvg(:,4,j));
    %laMean=nanmean(nanmean(samplesAvg(:,6:8,j)));
    
    %laMean=nanmean(nanmean(samplesAvg(:,4:6,j)));
    %laMean=getLamellarBackground(img,maskCell);
    
    for k=1:nWin
        y1Norm(k,j)=(y1(k,j)-offset)/(laMean-offset);
    end
end



