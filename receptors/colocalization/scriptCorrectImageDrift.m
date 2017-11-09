function [tracksFinalT] = scriptCorrectImageDrift(tracksFinal,MD)
[trackedFeatureInfo,~,~,numSegments,~] = convStruct2MatIgnoreMS(tracksFinal,[]);
trackedFeatureInfoC = trackedFeatureInfo;
trackLife = getTrackSEL(trackedFeatureInfo);
startPt=MD.processes_{1,1}.funParams_.firstImageNum;
endPt = MD.processes_{1,1}.funParams_.lastImageNum;
fullLength = endPt-startPt+1;
indx = find(trackLife(:,3)==fullLength); %Change to reflect full time span always
[gapInfo,gapInfoSpace] = findTrackGaps(trackedFeatureInfo);
test = ismember(indx,gapInfo(:,1));
indx(test) = [];
fullTrack = trackedFeatureInfo(indx,:);


ks = 8*startPt+1:8:8*endPt-7;
xDiff = cell(size(ks));
yDiff = cell(size(ks));

j =1;
for k = ks
%     exist = ~isnan(trackedFeatureInfo(:,k));
%     xDiff = fullTrack(:,k)-fullTrack(:,k-8);
%     yDiff = fullTrack(:,k+1)-fullTrack(:,k-7);
%     scatter(j*ones(length(yDiff),1),yDiff);
%     xDiff{j} = fullTrack(:,k)-fullTrack(:,k-8);
%     yDiff{j} = fullTrack(:,k+1)-fullTrack(:,k-7);
    
    xDiff{j} = fullTrack(:,k)-fullTrack(:,ks(1)-8);
    yDiff{j} = fullTrack(:,k+1)-fullTrack(:,ks(1)+1-8);
   
    
    xAvg = nanmedian(xDiff{j});
    yAvg = nanmedian(yDiff{j});
    trackedFeatureInfoC(:,k) = trackedFeatureInfo(:,k)-(xAvg*ones(size(trackedFeatureInfo,1),1));
    trackedFeatureInfoC(:,k+1) = trackedFeatureInfo(:,k+1)-(yAvg*ones(size(trackedFeatureInfo,1),1));
    j = j+1;
    
end

tracksFinalT =tracksFinal;

k =1;
j =1;
while k <= length(tracksFinal)%check appropriate length, might have to be trackedFeatureInfo..
    parts = numSegments(k);
    nPart = j:j+parts-1;
    if length(nPart) >1
       tracksFinalT(k).tracksCoordAmpCG = nan(length(nPart),8*length(min(trackLife(nPart,1)):max(trackLife(nPart,2))));
    end
%   tracksFinalT(k).tracksCoordAmpCG = trackedFeatureInfoC(nPart,(trackLife(nPart,1)*8)-7:(trackLife(nPart,2)*8));
 tracksFinalT(k).tracksCoordAmpCG = trackedFeatureInfoC(nPart,(8*min(trackLife(nPart,1))-7:8*max(trackLife(nPart,2))));
  j = j +parts;
  k = k+1;
end

end