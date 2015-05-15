function tracksNA = readIntensityFromTracks(tracksNA,imgStack, attribute, varargin)
% tracksNA = readIntensityFromTracks(tracksNA,imgStack, attribute) records
% pixel intensity from imgStack for entire track, even after ANA state in the last x,y
% position, and store it to attribute in tracksNA (1 means ampTotal, 2
% means forceMag
% assumes that tracks reflects coordinates from stage-drift-corrected
% images. imgStack is also from SDC output.
ip =inputParser;
ip.addParamValue('extraLength',0,@isscalar); % selcted track ids
ip.parse(varargin{:});
extraLength=ip.Results.extraLength;
% get stack size
numFrames = size(imgStack,3);
w4 = 8;
for k=1:numel(tracksNA)
%     startFrame = max(1, min(arrayfun(@(x) x.startingFrame,tracksNA))-extraLength);
%     endFrame = min(numFrames, max(arrayfun(@(x) x.endingFrame,tracksNA))+extraLength);
    startFrame = max(1, tracksNA(k).startingFrame-extraLength);
    endFrame = min(numFrames,tracksNA(k).endingFrame+extraLength);
    curStartingFrame = tracksNA(k).startingFrame;
    curEndingFrame = tracksNA(k).endingFrame;
    tracksNA(k).startingFrameExtra = startFrame;
    tracksNA(k).endingFrameExtra = endFrame;
    curRange = curStartingFrame:curEndingFrame;
    if extraLength==0
        if attribute==1
            tracksNA(k).ampTotal(curRange) = arrayfun(@(x) imgStack(round(tracksNA(k).yCoord(x)),round(tracksNA(k).xCoord(x)),x),curRange);
        elseif attribute==2
            tracksNA(k).forceMag(curRange) = arrayfun(@(x) imgStack(round(tracksNA(k).yCoord(x)),round(tracksNA(k).xCoord(x)),x),curRange);
        end
    else
        for ii=startFrame:endFrame
            curImg = imgStack(:,:,ii);
            if ii<curStartingFrame
                x = tracksNA(k).xCoord(curStartingFrame);
                y = tracksNA(k).yCoord(curStartingFrame);
                extrapolState=1;
            elseif ii>curEndingFrame
                x = tracksNA(k).xCoord(curEndingFrame);
                y = tracksNA(k).yCoord(curEndingFrame);
                extrapolState=1;
            else
                x = tracksNA(k).xCoord(ii);
                y = tracksNA(k).yCoord(ii);
                extrapolState=0;
            end
            if attribute==1 && extrapolState %intensity
                xi = floor(x);
                xres = x-xi;
                yi = floor(y);
                yres = y-yi;
                curAmpTotal = curImg(yi,xi);
                window = curImg(yi-w4:yi+w4, xi-w4:xi+w4);
                [prmVect]=fitGaussianMixture2D(window,[x-xi,y-yi,curAmpTotal-min(window(:)),1.35,min(window(:))],'xyas');
                if (abs(prmVect(1)-xres)<2 && abs(prmVect(2)-yres)<2 && prmVect(3)>0 && prmVect(3)<2*(curAmpTotal-min(window(:))))
                    tracksNA(k).xCoord(ii) = xi + prmVect(1);
                    tracksNA(k).yCoord(ii) = yi + prmVect(2);
                    tracksNA(k).amp(ii) = prmVect(3);
                    tracksNA(k).bkgAmp(ii) = prmVect(5);
                    tracksNA(k).ampTotal(ii) =  prmVect(3)+prmVect(5);
                else
                tracksNA(k).ampTotal(ii) = curAmpTotal;
            end
        elseif attribute==1 && ~extrapolState %intensity
            xi = round(x);
            yi = round(y);
            curAmpTotal = curImg(yi,xi);
            tracksNA(k).ampTotal(ii) = curAmpTotal;
        elseif attribute==2 %forceMag
            tracksNA(k).forceMag(ii) = curImg(round(y),round(x));
        elseif extrapolState
            disp('Please choose 1 or 2 for attribute.')
        end

        end
    end
end
end
