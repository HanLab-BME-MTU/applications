function newForceField = filterForceShortPeaks(forceField)
% filterForceShortPeaks(MD) filter out short peak from the force field
% shown in a few frames. Now it uses Gaussian low-pass filtering.
% input: forceField : forceField 
%           bgPeakLifeTime: life time critera for noise peak (default: 20)
% output: filtered forceField
nWindow = 3;
numNodes = length(forceField(1).pos);
numFrames = length(forceField);
progressText(0,'Median filtering force field')
newForceField = forceField;
splineParam = 0.01;
filterWindow=3;
tRange = 1:numFrames;

for k=1:numNodes
    % get the profile
    curVecX = arrayfun(@(x) x.vec(k,1),forceField);
    curVecY = arrayfun(@(x) x.vec(k,2),forceField);
    % Smooth these out first

    curVecX_med = medfilt1(curVecX,filterWindow);
    curVecY_med = medfilt1(curVecY,filterWindow);
    curVecX_spline= csaps(tRange,curVecX_med,splineParam);
    curVecX_spline_discretized=ppval(curVecX_spline,tRange);
    curVecY_spline= csaps(tRange,curVecY_med,splineParam);
    curVecY_spline_discretized=ppval(curVecY_spline,tRange);

    for ii=1:numFrames
        newForceField(ii).vec(k,:) = [curVecX_spline_discretized(ii), curVecY_spline_discretized(ii)];
    end
    progressText(k/numNodes,'Median filtering force field')
end
% if nargin <2
%     bgPeakLifeTime=20;
% end
% filter forcefield temporally
% for each node
% distinctiveness = 1.2; % if the temporal peak is twice more than the base, it'll be filtered out.
% numLMx = zeros(numNodes,1);
% numLMy = zeros(numNodes,1);
% forceMag = (forceField(1).vec(:,1).^2+forceField(1).vec(:,2).^2).^0.5;
% forceMagSorted = sort(forceMag);
% Let's use the double criteria: low bgLevel and short lifetime, and a
% bit high 
% bgForceLevel1 = forceMagSorted(round(0.3*length(forceMag)));
% bgForceLevel2 = forceMagSorted(round(0.8*length(forceMag)));
% bgPeakLifeTime = 20; % frames
%     curMag = (curVecX.^2+curVecY.^2).^0.5;
    %% Remove too high force peak spatially
%     % find segments that exceed the bgForceLevel
%     curMagExc = curMag > bgForceLevel1;
%     [curSegLabel,numSegs] = bwlabel(curMagExc);
%     % find the neighboring mean force
%     for ii=1:numSegs
%         if sum(curSegLabel==ii)<=bgPeakLifeTime/2
%             neighIdxBefore = find(curSegLabel==ii,1,'first')-1;
%             neighIdxAfter = find(curSegLabel==ii,1,'last')+1;
%             if neighIdxBefore<1
%                 neighForceX = curVecX(neighIdxAfter);
%                 neighForceY = curVecY(neighIdxAfter);
%             elseif neighIdxAfter>numFrames
%                 neighForceX = curVecX(neighIdxBefore);
%                 neighForceY = curVecY(neighIdxBefore);
%             else
%                 neighForceX = mean(curVecX([neighIdxBefore neighIdxAfter]));
%                 neighForceY = mean(curVecY([neighIdxBefore neighIdxAfter]));
%             end
%             curVecX(curSegLabel==ii) = neighForceX;
%             curVecX(curSegLabel==ii) = neighForceY;
%         end
%     end
    %% do lowpass filtering
%     % FFT the signal
% %     subplot(3,1,1), plot(1:numFrames,curVecX)
%     NFFT = 2^nextpow2(numFrames);
%     curFreqX = fft(curVecX,NFFT)/numFrames;
%     curFreqY = fft(curVecY,NFFT)/numFrames;
% %     Fs = 100;
%     highFreq = 21;
% %     freq = Fs/2*linspace(0,1,NFFT/2+1);
% %     subplot(3,1,3), plot(freq,2*abs(curFreqX(1:NFFT/2+1)))
%     % Low-pass filtering
%     curFreqX(highFreq:end-highFreq+1) = 0;
%     curFreqY(highFreq:end-highFreq+1) = 0;
% %     subplot(3,1,2), plot(freq,2*abs(curFreqX(1:NFFT/2+1)))
%     % inverse FT
%     curFreqXfiltered = ifft(curFreqX,NFFT)*numFrames;
%     curFreqYfiltered = ifft(curFreqY,NFFT)*numFrames;
% %     subplot(3,1,3), plot(1:numFrames, curFreqXfiltered(1:numFrames))
%     curVecX = real(curFreqXfiltered(1:numFrames));
%     curVecY = real(curFreqYfiltered(1:numFrames));
%     curMagExc = curMag > bgForceLevel2;
%     [curSegLabel,numSegs] = bwlabel(curMagExc);
%     % find the neighboring mean force
%     for ii=1:numSegs
%         if sum(curSegLabel==ii)<=bgPeakLifeTime
%             neighIdxBefore = find(curSegLabel==ii,1,'first')-1;
%             neighIdxAfter = find(curSegLabel==ii,1,'last')+1;
%             if neighIdxBefore<1
%                 neighForceX = curVecX(neighIdxAfter);
%                 neighForceY = curVecY(neighIdxAfter);
%             elseif neighIdxAfter>numFrames
%                 neighForceX = curVecX(neighIdxBefore);
%                 neighForceY = curVecY(neighIdxBefore);
%             else
%                 neighForceX = mean(curVecX([neighIdxBefore neighIdxAfter]));
%                 neighForceY = mean(curVecY([neighIdxBefore neighIdxAfter]));
%             end
%             curVecX(curSegLabel==ii) = neighForceX;
%             curVecX(curSegLabel==ii) = neighForceY;
%         end
%     end
%     
    %% Local maxima methods - it finds only peaks one by one, not in group
%     % find the local maxima
%     indLMx = locmax1d(curVecX,window);
%     indLMy = locmax1d(curVecY,window);
%     numLMx(k) = sum(indLMx);
%     numLMy(k) = sum(indLMy);
%     % see if the lm is above the distinctive range, and it only peaks for
%     % one frame, i.e., the two neighboring points before and after the lm
%     % should be in the baseline
%     for iix = indLMx'
%         if iix==1
%             neighForceX = mean(curVecX(iix+1:iix+(window-1)/2));
%         elseif iix == 2
%             neighForceX = mean(curVecX([iix-1, iix+1, iix+(window-1)/2]));
%         elseif iix == numFrames
%             neighForceX = mean(curVecX(iix-(window-1)/2:iix-1));
%         elseif iix == numFrames-1
%             neighForceX = mean(curVecX([iix-(window-1)/2, iix-1, iix+1]));
%         else
%             neighForceX = mean(curVecX([iix-(window-1)/2, iix-1, iix+1, iix+(window-1)/2]));
%         end
%         if curVecX(iix) > distinctiveness * neighForceX
%             curVecX(iix) = neighForceX;
%         end
%     end
% 
%     for iiy = indLMy'
%         if iiy==1
%             neighForceY = mean(curVecY(iiy+1:iiy+(window-1)/2));
%         elseif iiy == 2
%             neighForceY = mean(curVecY([iiy-1, iiy+1, iiy+(window-1)/2]));
%         elseif iiy == numFrames
%             neighForceY = mean(curVecY(iiy-(window-1)/2:iiy-1));
%         elseif iiy == numFrames-1
%             neighForceY = mean(curVecY([iiy-(window-1)/2, iiy-1, iiy+1]));
%         else
%             neighForceY = mean(curVecY([iiy-(window-1)/2, iiy-1, iiy+1, iiy+(window-1)/2]));
%         end
%         if curVecY(iiy) > distinctiveness * neighForceY
%             curVecY(iiy) = neighForceY;
%         end
%     end

