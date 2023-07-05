function h2 = showSingleAdhesionTrackSummaryFromOnlyTrack(curTrack,MD)
%  SHOWSINGLEADHESIONTRACKSUMMARYFROMONLYTRACK This function uses only curTrack to show the sammary of the adhesion track 
% using function showSingleAdhesionTrackSummary.
% 
% Detailed explanation of this function.
% curTrack      for example, tracksNA(1). You can get this by running
%               tracksNA = assembleTracksFromFAPackage(MD); And run
%               tracksNA.owner = MD;

if nargin<2
    % It assumes curTrack to contain MD as owner.
    if isfield(curTrack,'owner')
        MD = curTrack.owner;
    else
        disp('MD is not in a field of curTrack. Please open it from the browser.')
        [filename, pathname] = uigetfile( ...
           {'*.mat','MAT-files (*.mat)'}, ...
            'Pick a MD file', 'movieData.mat');
        MDstruct = load([pathname filesep filename]);
        MD = MDstruct.MD;
    end
end
% Now Have to read imgMap and/or tMap
%% persistent set up for large memory-requiring variable
persistent imgStack tMap imgStack2 labelFAs curChanPath

if isempty(imgStack) || ~strcmp(curChanPath, MD.channels_(1).channelPath_)
    [imgStack, tMap, imgStack2, labelFAs] = getAnyStacks(MD);
    curChanPath = MD.channels_(1).channelPath_;
end

% This is for testing and improvement. I want to re-track the curTrack and
% see how the bkgAmp is estimated.

curTrack2 = reestimateBkgFromTracks(curTrack, imgStack, labelFAs);
figure('Position',[50 100 600 700]); 
subplot(3,2,1), plot(curTrack.bkgAmp),title('bkgAmp')
subplot(3,2,3), plot(curTrack.amp), title('Amp')
subplot(3,2,5), plot(curTrack.ampTotal), title('AmpTotal')

subplot(3,2,2), plot(curTrack2.bkgAmp),title('bkgAmp')
subplot(3,2,4), plot(curTrack2.amp), title('Amp')
subplot(3,2,6), plot(curTrack2.ampTotal), title('AmpTotal')

h2 = showSingleAdhesionTrackSummaryRateConstFitting(MD,curTrack2,imgStack,tMap,imgStack2);
% function [h2, timeLagMasterAgainstForce,timeLagMasterAgainstMainSlave] = ...
%     showSingleAdhesionTrackSummary(MD,curTrack,imgMap,tMap,imgMap2,IDtoInspect, ...
%     gPath,additionalName,imgStackBS,imgStackBS2)
end