function tracksNA = readIntensityFromTracks(tracksNA,imgStack, attribute)
% tracksNA = readIntensityFromTracks(tracksNA,imgStack, attribute) records
% pixel intensity from imgStack for entire track, even after ANA state in the last x,y
% position, and store it to attribute in tracksNA (1 means ampTotal, 2
% means forceMag
% assumes that tracks reflects coordinates from stage-drift-corrected
% images. imgStack is also from SDC output.

% get stack size
N = size(imgStack,3);
for ii=1:N
    curImg = imgStack(:,:,ii);
    for k=1:numel(tracksNA)
        curStartingFrame = tracksNA(k).startingFrame;
        curEndingFrame = tracksNA(k).endingFrame;
        if ii<curStartingFrame
            x = tracksNA(k).xCoord(curStartingFrame);
            y = tracksNA(k).yCoord(curStartingFrame);
        elseif ii>curEndingFrame
            x = tracksNA(k).xCoord(curEndingFrame);
            y = tracksNA(k).yCoord(curEndingFrame);
        else
            x = tracksNA(k).xCoord(ii);
            y = tracksNA(k).yCoord(ii);
        end
        if attribute==1 %intensity
            tracksNA(k).ampTotal(ii) = curImg(round(y),round(x));
        elseif attribute==2 %forceMag
            tracksNA(k).forceMag(ii) = curImg(round(y),round(x));
        else
            disp('Please choose 1 or 2 for attribute.')
        end
    end
end
