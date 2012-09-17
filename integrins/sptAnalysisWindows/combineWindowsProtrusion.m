function [protrusionCombinedWindows,posInfoCombinedWindows] = ...
    combineWindowsProtrusion(sliceActivityGroup,winPositions,...
    firstMaskFile,protSamples,funcParam)


%REMARKS There are various hard-coded numbers in the function. Should
%        modify for better programming practice.
%
%Khuloud Jaqaman, July 2012

%% Input

%read cell mask images
[fpath,fname,fno,fext]=getFilenameBody(firstMaskFile);
dirName=[fpath,filesep];
fName=[fname,fno,fext];
outFileList = getFileStackNames([dirName,fName]);
numFiles = length(outFileList);
maskImage = imread(outFileList{1});
maskImage = repmat(maskImage,[1 1 numFiles]);
for iFile = 2 : numFiles
    maskImage(:,:,iFile) = imread(outFileList{iFile});
end
maskImage = double(maskImage);

%get protrusion normals
avgNormal = protSamples.avgNormal;

%determine image size and number of frames
[imageSizeX,imageSizeY,numFrames] = size(maskImage);

%generate pixel coordinates
pixelCoordX = repmat((1:imageSizeX)',imageSizeY,1);
pixelCoordY = repmat((1:imageSizeY),imageSizeX,1);
pixelCoordY = pixelCoordY(:);

%hard-coded numbers
numTypeProt = funcParam.numTypeProt; %9; %number of protrusion types, classified based on the events before protrusion
numPixInBand = funcParam.numPixInBand; %2; %number of pixels in each band, as used throughout the analysis
numSmallBandsInBigBand = funcParam.numSmallBandsInBigBand; %5; %number of bands making ~1 um, based on the fact that 1 pixel = 111 nm and each band is ~2 pixels deep
numBigBands = funcParam.numBigBands; %10; %number of bands of bands to look inside the cell; 10 "bands of bands" looks into ~10 um from cell edge
maxNegInc = funcParam.maxNegInc; %3; %number of time points to look before protrusion onset
maxPosInc = funcParam.maxPosInc; %numFrames - 1; %number of time points to look after protrusion onset (i.e. during protrusion)

%% Combine windows

%initialize output
protrusionCombinedWindows = repmat(struct('eventInfo',[]),numTypeProt,1);
posInfoCombinedWindows = repmat(struct('posInfo',[]),numTypeProt,1);

%go over all protrusion types
for iType = 1 : numTypeProt
    
    %get array storing protrusion events
    protEvents = sliceActivityGroup(iType).edgeClassInfo;
    numEvents = size(protEvents,1);
    
    %initialize structure array for output
    eventInfo = repmat(struct('onset',[],'befStatic',[],'befDynamic',[],...
        'aftStatic',[],'aftDynamic',[]),numEvents,1);
    posInfo = repmat(struct('onset',[],'befStatic',[],'befDynamic',[],...
        'aftStatic',[],'aftDynamic',[]),numEvents,1);
    
    %go over each protrusion event and assemble its window time series
    %before and after protrusion onset
    for iEvent = 1 : numEvents
        
        %get event information
        onsetFrame = protEvents(iEvent,1);
        sliceIndx = protEvents(iEvent,2);
        eventDur = protEvents(iEvent,3);
        prevEventDur = protEvents(iEvent,4);
        prevEventDur(isnan(prevEventDur)) = 0;
        
        %put together windows at protrusion onset in big bands
        indxWindowsOnset = cell(1,numBigBands);
        posWindowsOnset = NaN(1,numBigBands,2);
        nBands = length(winPositions{onsetFrame,sliceIndx});
        for iBigBand = 1 : numBigBands
            bandStartCurrent = (iBigBand-1)*numSmallBandsInBigBand + 1;
            bandEndCurrent = iBigBand*numSmallBandsInBigBand;
            if bandStartCurrent <= nBands
                bandEndCurrent = min(nBands,bandEndCurrent);
                indxWindowsOnset{iBigBand} = [(bandStartCurrent:bandEndCurrent)' ...
                    repmat([sliceIndx onsetFrame onsetFrame],bandEndCurrent-bandStartCurrent+1,1)];
                posWindowsOnset(1,iBigBand,:) = [(bandStartCurrent-1) bandEndCurrent]*numPixInBand;
            end
        end
        
        %put together windows before protrusion onset
        indxWindowsBefDynamic = cell(maxNegInc,numBigBands);
        indxWindowsBefStatic = cell(maxNegInc,numBigBands);
        posWindowsBefDynamic = NaN(maxNegInc,numBigBands,2);
        posWindowsBefStatic = NaN(maxNegInc,numBigBands,2);
        for iInc = 1 : min(prevEventDur,maxNegInc)
            
            %specify current frame
            frameCurrent = onsetFrame - iInc;
            
            %put together windows in big bands
            nBands = length(winPositions{frameCurrent,sliceIndx});
            for iBigBand = 1 : numBigBands
                
                %dynamic windows that move with the cell edge
                bandStartCurrent = (iBigBand-1)*numSmallBandsInBigBand + 1;
                bandEndCurrent = iBigBand*numSmallBandsInBigBand;
                if bandStartCurrent <= nBands
                    bandEndCurrent = min(nBands,bandEndCurrent);
                    indxWindowsBefDynamic{iInc,iBigBand} = [(bandStartCurrent:bandEndCurrent)' ...
                        repmat([sliceIndx frameCurrent frameCurrent],bandEndCurrent-bandStartCurrent+1,1)];
                    posWindowsBefDynamic(iInc,iBigBand,:) = [(bandStartCurrent-1) bandEndCurrent]*numPixInBand;
                end
                
                %static windows fixed at location of protrusion onset
                tmp = indxWindowsOnset{iBigBand};
                if ~isempty(tmp)
                    tmp(:,4) = frameCurrent;
                    indxWindowsBefStatic{iInc,iBigBand} = tmp;
                    posWindowsBefStatic(iInc,iBigBand,:) = posWindowsOnset(1,iBigBand,:) - sum(avgNormal(sliceIndx,frameCurrent:onsetFrame-1));
                end
                
            end
            
        end
        
        %put together windows after protrusion onset
        indxWindowsAftStatic = cell(maxPosInc,numBigBands);
        indxWindowsAftDynamic = cell(maxPosInc);
        posWindowsAftStatic = NaN(maxPosInc,numBigBands,2);
        posWindowsAftDynamic = NaN(maxPosInc,maxPosInc,2);
        for iInc = 1 : eventDur - 1
            
            %specify current frame
            frameCurrent = onsetFrame + iInc;

            %fixed location of protrusion onset
            for iBigBand = 1 : numBigBands
                tmp = indxWindowsOnset{iBigBand};
                if ~isempty(tmp)
                    tmp(:,4) = frameCurrent;
                    indxWindowsAftStatic{iInc,iBigBand} = tmp;
                    posWindowsAftStatic(iInc,iBigBand,:) = posWindowsOnset(1,iBigBand,:) + sum(avgNormal(sliceIndx,onsetFrame:frameCurrent-1));
                end
            end
            
            %previous "new" cell-ECM contact areas
            for iNew = 1 : iInc - 1
                tmp = indxWindowsAftDynamic{iNew,iNew};
                if ~isempty(tmp)
                    tmp(:,4) = frameCurrent;
                    indxWindowsAftDynamic{iInc,iNew} = tmp;
                    posWindowsAftDynamic(iInc,iNew,:) = posWindowsAftDynamic(iNew,iNew,:) + sum(avgNormal(sliceIndx,onsetFrame+iNew:frameCurrent-1));
                end
            end
            
            %new cell-ECM contact area
            nBands = length(winPositions{frameCurrent,sliceIndx});
            iBand = 0;
            newArea = 1;
            while newArea && (iBand < nBands)
                
                %update band number
                iBand = iBand + 1;
                
                %get position of window in current band
                windowPoly = [winPositions{frameCurrent,sliceIndx}{iBand}{:}];
                winX = windowPoly(1,:);
                winY = windowPoly(2,:);
                
                %find pixels inside window
                pixelsWindow = find(inpolygon(pixelCoordX,pixelCoordY,winY,winX));
                numPixWindow = length(pixelsWindow);
                xCoordWindow = pixelCoordX(pixelsWindow);
                yCoordWindow = pixelCoordY(pixelsWindow);
                timeWindow = (frameCurrent - 1)*ones(numPixWindow,1);
                
                %determine which of these pixels represent new cell-ECM
                %contact area
                linearIndx = sub2ind([imageSizeX imageSizeY numFrames],...
                    xCoordWindow,yCoordWindow,timeWindow);
                numPixelsNew = length(find(maskImage(linearIndx)==0));
                
                %if new area >= 0.5 window size, set newArea to 1 in order
                %to look into next band
                newArea = (numPixelsNew >= 0.5*numPixWindow);
                
            end
            
            %store window information
            maxBand = iBand - 1;
            if maxBand >= 1
                indxWindowsAftDynamic{iInc,iInc} = [(1:maxBand)' ...
                    repmat([sliceIndx frameCurrent frameCurrent],maxBand,1)];
                posWindowsAftDynamic(iInc,iInc,:) = [0 avgNormal(sliceIndx,frameCurrent-1)];
            end
            
        end %(for iInc = 1 : eventDur - 1)
        
        %combine all new cell area events into one series with the events
        %aligned by their start frames
        indxWindowsAftDynamicComb = cell(maxPosInc,1);
        posWindowsAftDynamicComb = NaN(maxPosInc,1,2);
        for iInc = 1 : eventDur - 1
            colIndx = (1 : eventDur-iInc)';
            rowIndx = colIndx + iInc - 1;
            linIndx = sub2ind([maxPosInc maxPosInc],rowIndx,colIndx);
            indxWindowsAftDynamicComb{iInc} = vertcat(indxWindowsAftDynamic{linIndx});
            tmp = posWindowsAftDynamic(:,:,1);
            posWindowsAftDynamicComb(iInc,1,1) = mean(tmp(linIndx));
            tmp = posWindowsAftDynamic(:,:,2);
            posWindowsAftDynamicComb(iInc,1,2) = mean(tmp(linIndx));
        end
        tmp = mean([posWindowsAftDynamicComb(2:eventDur-1,1,1) posWindowsAftDynamicComb(1:eventDur-2,1,2)],2);
        posWindowsAftDynamicComb(2:eventDur-1,1,1) = tmp;
        posWindowsAftDynamicComb(1:eventDur-2,1,2) = tmp;
        
        %store information for output
        eventInfo(iEvent).onset = indxWindowsOnset;
        eventInfo(iEvent).befStatic = indxWindowsBefStatic;
        eventInfo(iEvent).befDynamic = indxWindowsBefDynamic;
        eventInfo(iEvent).aftStatic = indxWindowsAftStatic;
        eventInfo(iEvent).aftDynamic = indxWindowsAftDynamic;
        eventInfo(iEvent).aftDynamicComb = indxWindowsAftDynamicComb;
        posInfo(iEvent).onset = posWindowsOnset;
        posInfo(iEvent).befStatic = posWindowsBefStatic;
        posInfo(iEvent).befDynamic = posWindowsBefDynamic;
        posInfo(iEvent).aftStatic = posWindowsAftStatic;
        posInfo(iEvent).aftDynamic = posWindowsAftDynamic;
        posInfo(iEvent).aftDynamicComb = posWindowsAftDynamicComb;
        
    end %(for iEvent = 1 : numEvents)
    
    %final output
    protrusionCombinedWindows(iType).eventInfo = eventInfo;
    posInfoCombinedWindows(iType).posInfo = posInfo;
    
end %(for iType = 1 : numTypeProt)

%% ~~~ the end ~~~

