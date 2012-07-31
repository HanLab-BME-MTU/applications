function protrusionCombinedWindows = combineWindowsProtrusion(sliceActivityGroup,...
    winPositions,firstMaskFile)


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

%determine image size and number of frames
[imageSizeX,imageSizeY,numFrames] = size(maskImage);

%generate pixel coordinates
pixelCoordX = repmat((1:imageSizeX)',imageSizeY,1);
pixelCoordY = repmat((1:imageSizeY),imageSizeX,1);
pixelCoordY = pixelCoordY(:);

%hard-coded numbers
numTypeProt = 9;
numBands1um = 5;
maxNegInc = 3;
maxPosInc = numFrames - 1;

%% Combine windows

%initialize output
protrusionCombinedWindows = repmat(struct('eventInfo',[]),numTypeProt,1);

%go over all protrusion types
for iType = 1 : numTypeProt
    
    %get array storing protrusion events
    protEvents = sliceActivityGroup(iType).edgeClassInfo;
    numEvents = size(protEvents,1);
    
    %initialize structure array for output
    eventInfo = repmat(struct('onset',[],'befStatic',[],'befDynamic',[],...
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
        
        %put together windows upto ~1 um from cell edge at protrusion onset
        nBands = length(winPositions{onsetFrame,sliceIndx});
        nBands = min(nBands,numBands1um);
        indxWindowsOnset = {[(1:nBands)' repmat([sliceIndx onsetFrame onsetFrame],nBands,1)]};
        
        %put together windows before protrusion onset
        indxWindowsBefDynamic = cell(maxNegInc,1);
        indxWindowsBefStatic = cell(maxNegInc,1);
        for iInc = 1 : min(prevEventDur,maxNegInc)
            
            %specify current frame
            frameCurrent = onsetFrame - iInc;

            %dynamic windows that move with the cell edge
            nBands = length(winPositions{frameCurrent,sliceIndx});
            nBands = min(nBands,numBands1um);
            indxWindowsBefDynamic{iInc} = [(1:nBands)' repmat([sliceIndx frameCurrent frameCurrent],nBands,1)];
            
            %static windows fixed at location of protrusion onset
            tmp = indxWindowsOnset{1};
            tmp(:,4) = frameCurrent;
            indxWindowsBefStatic{iInc} = tmp;
            
        end
        
        %put together windows after protrusion onset
        indxWindowsAftStatic = cell(maxPosInc,1);
        indxWindowsAftDynamic = cell(maxPosInc);
        for iInc = 1 : eventDur - 1
            
            %specify current frame
            frameCurrent = onsetFrame + iInc;

            %fixed location of protrusion onset
            tmp = indxWindowsOnset{1};
            tmp(:,4) = frameCurrent;
            indxWindowsAftStatic{iInc} = tmp;
            
            %previous "new" cell-ECM contact areas
            for iNew = 1 : iInc - 1
                tmp = indxWindowsAftDynamic{iNew,iNew};
                if ~isempty(tmp)
                    tmp(:,4) = frameCurrent;
                    indxWindowsAftDynamic{iInc,iNew} = tmp;
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
            end
            
        end %(for iInc = 1 : eventDur - 1)
        
        %combine all new cell area events into one series with the events
        %aligned by their start frames
        indxWindowsAftDynamicComb = cell(maxPosInc,1);
        for iInc = 1 : eventDur - 1
            colIndx = (1 : eventDur-iInc)';
            rowIndx = colIndx + iInc - 1;
            linIndx = sub2ind([maxPosInc maxPosInc],rowIndx,colIndx);
            indxWindowsAftDynamicComb{iInc} = vertcat(indxWindowsAftDynamic{linIndx});
        end
        
        %store information for output
        eventInfo(iEvent).onset = indxWindowsOnset;
        eventInfo(iEvent).befStatic = indxWindowsBefStatic;
        eventInfo(iEvent).befDynamic = indxWindowsBefDynamic;
        eventInfo(iEvent).aftStatic = indxWindowsAftStatic;
        eventInfo(iEvent).aftDynamic = indxWindowsAftDynamic;
        eventInfo(iEvent).aftDynamicComb = indxWindowsAftDynamicComb;
        
    end %(for iEvent = 1 : numEvents)
    
    %final output
    protrusionCombinedWindows(iType).eventInfo = eventInfo;
    
end %(for iType = 1 : numTypeProt)

%% ~~~ the end ~~~

