function forceField = filterForceShortPeaks(forceField,bgPeakLifeTime)
% filterForceShortPeaks(MD) filter out short peak from the force field
% shown in a few frames.
% input: forceField : forceField 
%           bgPeakLifeTime: life time critera for noise peak (default: 20)
% output: filtered forceField
if nargin <2
    bgPeakLifeTime=20;
end
% filter forcefield temporally
% for each node
distinctiveness = 1.2; % if the temporal peak is twice more than the base, it'll be filtered out.
window = 3;
numNodes = length(forceField(1).pos);
numLMx = zeros(numNodes,1);
numLMy = zeros(numNodes,1);
numFrames = length(forceField);
forceMag = (forceField(1).vec(:,1).^2+forceField(1).vec(:,2).^2).^0.5;
forceMagSorted = sort(forceMag);
% Let's use the double criteria: low bgLevel and short lifetime, and a
% bit high 
bgForceLevel1 = forceMagSorted(round(0.3*length(forceMag)));
bgForceLevel2 = forceMagSorted(round(0.8*length(forceMag)));
% bgPeakLifeTime = 20; % frames
for k=1:numNodes
    % get the profile
    curVecX = arrayfun(@(x) x.vec(k,1),forceField);
    curVecY = arrayfun(@(x) x.vec(k,2),forceField);
    curMag = (curVecX.^2+curVecY.^2).^0.5;
    % find segments that exceed the bgForceLevel
    curMagExc = curMag > bgForceLevel1;
    [curSegLabel,numSegs] = bwlabel(curMagExc);
    % find the neighboring mean force
    for ii=1:numSegs
        if sum(curSegLabel==ii)<=bgPeakLifeTime/2
            neighIdxBefore = find(curSegLabel==ii,1,'first')-1;
            neighIdxAfter = find(curSegLabel==ii,1,'last')+1;
            if neighIdxBefore<1
                neighForceX = curVecX(neighIdxAfter);
                neighForceY = curVecY(neighIdxAfter);
            elseif neighIdxAfter>numFrames
                neighForceX = curVecX(neighIdxBefore);
                neighForceY = curVecY(neighIdxBefore);
            else
                neighForceX = mean(curVecX([neighIdxBefore neighIdxAfter]));
                neighForceY = mean(curVecY([neighIdxBefore neighIdxAfter]));
            end
            curVecX(curSegLabel==ii) = neighForceX;
            curVecX(curSegLabel==ii) = neighForceY;
        end
    end
    curMagExc = curMag > bgForceLevel2;
    [curSegLabel,numSegs] = bwlabel(curMagExc);
    % find the neighboring mean force
    for ii=1:numSegs
        if sum(curSegLabel==ii)<=bgPeakLifeTime
            neighIdxBefore = find(curSegLabel==ii,1,'first')-1;
            neighIdxAfter = find(curSegLabel==ii,1,'last')+1;
            if neighIdxBefore<1
                neighForceX = curVecX(neighIdxAfter);
                neighForceY = curVecY(neighIdxAfter);
            elseif neighIdxAfter>numFrames
                neighForceX = curVecX(neighIdxBefore);
                neighForceY = curVecY(neighIdxBefore);
            else
                neighForceX = mean(curVecX([neighIdxBefore neighIdxAfter]));
                neighForceY = mean(curVecY([neighIdxBefore neighIdxAfter]));
            end
            curVecX(curSegLabel==ii) = neighForceX;
            curVecX(curSegLabel==ii) = neighForceY;
        end
    end
    
    %Local maxima methods - it finds only peaks one by one, not in group
    % find the local maxima
    indLMx = locmax1d(curVecX,window);
    indLMy = locmax1d(curVecY,window);
    numLMx(k) = sum(indLMx);
    numLMy(k) = sum(indLMy);
    % see if the lm is above the distinctive range, and it only peaks for
    % one frame, i.e., the two neighboring points before and after the lm
    % should be in the baseline
    for iix = indLMx'
        if iix==1
            neighForceX = mean(curVecX(iix+1:iix+(window-1)/2));
        elseif iix == 2
            neighForceX = mean(curVecX([iix-1, iix+1, iix+(window-1)/2]));
        elseif iix == numFrames
            neighForceX = mean(curVecX(iix-(window-1)/2:iix-1));
        elseif iix == numFrames-1
            neighForceX = mean(curVecX([iix-(window-1)/2, iix-1, iix+1]));
        else
            neighForceX = mean(curVecX([iix-(window-1)/2, iix-1, iix+1, iix+(window-1)/2]));
        end
        if curVecX(iix) > distinctiveness * neighForceX
            curVecX(iix) = neighForceX;
        end
    end

    for iiy = indLMy'
        if iiy==1
            neighForceY = mean(curVecY(iiy+1:iiy+(window-1)/2));
        elseif iiy == 2
            neighForceY = mean(curVecY([iiy-1, iiy+1, iiy+(window-1)/2]));
        elseif iiy == numFrames
            neighForceY = mean(curVecY(iiy-(window-1)/2:iiy-1));
        elseif iiy == numFrames-1
            neighForceY = mean(curVecY([iiy-(window-1)/2, iiy-1, iiy+1]));
        else
            neighForceY = mean(curVecY([iiy-(window-1)/2, iiy-1, iiy+1, iiy+(window-1)/2]));
        end
        if curVecY(iiy) > distinctiveness * neighForceY
            curVecY(iiy) = neighForceY;
        end
    end

    for ii=1:numFrames
        forceField(ii).vec(k,:) = [curVecX(ii), curVecY(ii)];
    end
end
