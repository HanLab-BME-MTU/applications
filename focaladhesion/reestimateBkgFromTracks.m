function tracksNA = reestimateBkgFromTracks(tracksNA, imgStack, labelAdhesion)
%function curTrack2 = reestimateBkgFromTracks(curTrack, imgStack, maskFAs)
%reestimates bkgAmp from the adhesion segmentation.
%   input:
%       tracksNA            tracksNA track
%       imgStack            image stack where the signal is read
%       maskFAs             mask of focal adhesions
%   output:
%       tracksNA            same as curTrack but it has updated values for
%                           bkgAmp and amp.
%  Sangyoon Han, July 4, 2023

% go with curTrack
nTracks = numel(tracksNA);
dilateDist = 2;
dilateDistBG = 2;
adhRad = 3; %diffraction limited NA radius (generous estimation)
for ii= 1:nTracks
    % now go with frames
    curTrack = tracksNA(ii);
    for jj = curTrack.startingFrameExtraExtra:curTrack.endingFrameExtra
        % get the associated FAs
        curImg = imgStack(:,:,jj);
        if curTrack.refineFAID(jj)>0
            curMaskAdh = labelAdhesion(:,:,jj)==curTrack.refineFAID(jj); %it matched! It need to be re-segemented though because it's very small
        else %When there is no associated FA segmentation, we will use NA posotion
            masksize = size(labelAdhesion(:,:,jj));
            [x, y] = meshgrid(1:masksize(2), 1:masksize(1));  %create pixel coordinates
            curMaskAdh = hypot(x - curTrack.xCoord(jj), y - curTrack.yCoord(jj)) <= adhRad;  %require R2016b or later
%             figure, imshow(mask)
        end
            
        % dilate and make a band mask
        adhDilated = bwmorph(curMaskAdh,'dilate',dilateDist);
        %overlayMaskBoundaries(adhDilated,'r');
        bg= bwmorph(adhDilated,'dilate',dilateDistBG);
        % remove other FG from the band
        bkgMask = bg & ~(labelAdhesion(:,:,jj)>0) & ~adhDilated; 

%         figure, imshow(curImg,[])
%         hold on
%         plot(curTrack.xCoord(jj),curTrack.yCoord(jj),'co')
%         overlayMaskBoundaries(bkgMask,'k');

        % bkg estimation
        bkgAvg = mean(curImg(bkgMask));
%         curTrack.bkgAmp(jj)
%         bkgAvg
        curTrack.bkgAmp(jj) = bkgAvg;
        curTrack.ampTotal(jj) = mean(curImg(curMaskAdh));
        curTrack.amp(jj) = curTrack.ampTotal(jj) - bkgAvg;
    end
    tracksNA(ii) = curTrack;
end

end

