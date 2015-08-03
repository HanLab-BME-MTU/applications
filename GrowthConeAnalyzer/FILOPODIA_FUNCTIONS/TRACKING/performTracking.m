function [ ] = performTracking(analInfo,costMatParams,imDir,saveDir)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

 
% convert filoInfo such that if nan fit just use 

[movieInfo,analInfoFilt] = convertAnalInfo2MovieInfo(analInfo); 

[trackedFeatureIndx,trackedFeatureInfolink] = linkFeaturesMODFORFILO(movieInfo,costMatParams,analInfoFilt,[],2,[],[],[],saveDir);
save([saveDir filesep 'tracking.mat'],'trackedFeatureIndx','trackedFeatureInfolink','costMatParams','analInfoFilt'); 
makeTrackMovie(trackedFeatureInfolink,trackedFeatureIndx, analInfoFilt, imDir,saveDir,5,0.216); 
% trackingOutput = plotHists(trackedFeatureInfo,saveDir); 
end 

