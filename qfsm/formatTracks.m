function newTracks = formatTracks(tracks,detectedPts,nFrames,T)
% Format tracks structure into tracks with every frame
if nargin<4
    T=zeros(nFrames,2); % T is a translation matrix
end
% Get limits of transformation array
maxX = ceil(max(abs(T(:, 2))));
maxY = ceil(max(abs(T(:, 1))));


newTracks(numel(tracks),1) = struct('xCoord', [], 'yCoord', [],'iFrame',[],'presence',[],'amp',[],'bkgAmp',[]);
% BA: before adhesion, NA: nascent adh, FC: focal complex, FA: focal adh,
% ANA: after NA (failed to be matured.
for i = 1:numel(tracks)
    % Get the x and y coordinate of all compound tracks
    startNA = true;
    endNA = true;
    for  j = 1 : nFrames
        newTracks(i).iFrame(j) = j;
        if j<tracks(i).seqOfEvents(1,1)
            newTracks(i).xCoord(j) = NaN;
            newTracks(i).yCoord(j) = NaN;
            newTracks(i).presence(j) = false;
            newTracks(i).amp(j) = NaN;
        elseif j>tracks(i).seqOfEvents(2,1)
            newTracks(i).xCoord(j) = NaN;
            newTracks(i).yCoord(j) = NaN;
            newTracks(i).amp(j) = NaN;
            newTracks(i).presence(j) = false;
            if endNA
                newTracks(i).endingFrame = j-1;
                endNA = false;
            end
        elseif j==tracks(i).seqOfEvents(2,1)
            newTracks(i).xCoord(j) = tracks(i).tracksCoordAmpCG(1,1+8*(j-tracks(i).seqOfEvents(1,1)))+T(j,2)+maxX;
            newTracks(i).yCoord(j) = tracks(i).tracksCoordAmpCG(1,2+8*(j-tracks(i).seqOfEvents(1,1)))+T(j,1)+maxY;
            newTracks(i).amp(j) = tracks(i).tracksCoordAmpCG(1,4+8*(j-tracks(i).seqOfEvents(1,1)));
            if tracks(i).tracksFeatIndxCG(j-tracks(i).seqOfEvents(1,1)+1)==0
                newTracks(i).bkgAmp(j) = NaN;
            else
                newTracks(i).bkgAmp(j) = detectedPts(j-tracks(i).seqOfEvents(1,1)+1).bkg(tracks(i).tracksFeatIndxCG(j-tracks(i).seqOfEvents(1,1)+1));
                newTracks(i).sigma(j) = detectedPts(j-tracks(i).seqOfEvents(1,1)+1).sigmaX(tracks(i).tracksFeatIndxCG(j-tracks(i).seqOfEvents(1,1)+1));
            end
            newTracks(i).presence(j) = true;
            if endNA
                newTracks(i).endingFrame = j;
                endNA = false;
            end
        else
            newTracks(i).xCoord(j) = tracks(i).tracksCoordAmpCG(1,1+8*(j-tracks(i).seqOfEvents(1,1)))+T(j,2)+maxX;
            newTracks(i).yCoord(j) = tracks(i).tracksCoordAmpCG(1,2+8*(j-tracks(i).seqOfEvents(1,1)))+T(j,1)+maxY;
            newTracks(i).amp(j) = tracks(i).tracksCoordAmpCG(1,4+8*(j-tracks(i).seqOfEvents(1,1)));
            if tracks(i).tracksFeatIndxCG(j-tracks(i).seqOfEvents(1,1)+1)==0
                newTracks(i).bkgAmp(j) = NaN;
            else
                newTracks(i).bkgAmp(j) = detectedPts(j-tracks(i).seqOfEvents(1,1)+1).bkg(tracks(i).tracksFeatIndxCG(j-tracks(i).seqOfEvents(1,1)+1));
                newTracks(i).sigma(j) = detectedPts(j-tracks(i).seqOfEvents(1,1)+1).sigmaX(tracks(i).tracksFeatIndxCG(j-tracks(i).seqOfEvents(1,1)+1));
            end
            newTracks(i).presence(j) = true;
            if startNA
                newTracks(i).startingFrame = j;
                startNA = false;
            end
        end
            
        if isfield(tracks, 'label'),
            newTracks(iTrack).label = tracks(i).label;
        end
    end
    % go through frames again and fill NaNs with numbers at the gap
    % position
    for j=1:nFrames-1
        if j<nFrames-9 && sum(newTracks(i).presence(j:j+9))==10 ...
                && sum(isnan(newTracks(i).xCoord(j:j+9)))==10 
            gap = 10;
            for kk=1:gap
                newTracks(i).xCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).xCoord(j-1)+kk*newTracks(i).xCoord(j+gap))/(gap+1);
                newTracks(i).yCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).yCoord(j-1)+kk*newTracks(i).yCoord(j+gap))/(gap+1);
                newTracks(i).amp(j+kk-1) = ((gap+1-kk)*newTracks(i).amp(j-1)+kk*newTracks(i).amp(j+gap))/(gap+1);
                newTracks(i).bkgAmp(j+kk-1) = ((gap+1-kk)*newTracks(i).bkgAmp(j-1)+kk*newTracks(i).bkgAmp(j+gap))/(gap+1);
            end
        elseif j<nFrames-8 && sum(newTracks(i).presence(j:j+8))==9 ...
                && sum(isnan(newTracks(i).xCoord(j:j+8)))==9 
            gap = 9;
            for kk=1:gap
                newTracks(i).xCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).xCoord(j-1)+kk*newTracks(i).xCoord(j+gap))/(gap+1);
                newTracks(i).yCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).yCoord(j-1)+kk*newTracks(i).yCoord(j+gap))/(gap+1);
                newTracks(i).amp(j+kk-1) = ((gap+1-kk)*newTracks(i).amp(j-1)+kk*newTracks(i).amp(j+gap))/(gap+1);
                newTracks(i).bkgAmp(j+kk-1) = ((gap+1-kk)*newTracks(i).bkgAmp(j-1)+kk*newTracks(i).bkgAmp(j+gap))/(gap+1);
            end
        elseif j<nFrames-7 && sum(newTracks(i).presence(j:j+7))==8 ...
                && sum(isnan(newTracks(i).xCoord(j:j+7)))==8 
            gap = 8;
            for kk=1:gap
                newTracks(i).xCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).xCoord(j-1)+kk*newTracks(i).xCoord(j+gap))/(gap+1);
                newTracks(i).yCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).yCoord(j-1)+kk*newTracks(i).yCoord(j+gap))/(gap+1);
                newTracks(i).amp(j+kk-1) = ((gap+1-kk)*newTracks(i).amp(j-1)+kk*newTracks(i).amp(j+gap))/(gap+1);
                newTracks(i).bkgAmp(j+kk-1) = ((gap+1-kk)*newTracks(i).bkgAmp(j-1)+kk*newTracks(i).bkgAmp(j+gap))/(gap+1);
            end
        elseif j<nFrames-6 && sum(newTracks(i).presence(j:j+6))==7 ...
                && sum(isnan(newTracks(i).xCoord(j:j+6)))==7 
            gap = 7;
            for kk=1:gap
                newTracks(i).xCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).xCoord(j-1)+kk*newTracks(i).xCoord(j+gap))/(gap+1);
                newTracks(i).yCoord(j+kk-1) = ((gap+1-kk)*newTracks(i).yCoord(j-1)+kk*newTracks(i).yCoord(j+gap))/(gap+1);
                newTracks(i).amp(j+kk-1) = ((gap+1-kk)*newTracks(i).amp(j-1)+kk*newTracks(i).amp(j+gap))/(gap+1);
                newTracks(i).bkgAmp(j+kk-1) = ((gap+1-kk)*newTracks(i).bkgAmp(j-1)+kk*newTracks(i).bkgAmp(j+gap))/(gap+1);
            end
        elseif j<nFrames-5 && newTracks(i).presence(j) && newTracks(i).presence(j+1) && newTracks(i).presence(j+2) && newTracks(i).presence(j+3) ...
               && newTracks(i).presence(j+4) && newTracks(i).presence(j+5) && isnan(newTracks(i).xCoord(j)) ...
               && isnan(newTracks(i).xCoord(j+1)) && isnan(newTracks(i).xCoord(j+2)) && isnan(newTracks(i).xCoord(j+3))...
               && isnan(newTracks(i).xCoord(j+4)) && isnan(newTracks(i).xCoord(j+5))
            newTracks(i).xCoord(j) = (6*newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+6))/7;
            newTracks(i).yCoord(j) = (6*newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+6))/7;
            newTracks(i).amp(j) = (6*newTracks(i).amp(j-1)+newTracks(i).amp(j+6))/7;
            newTracks(i).bkgAmp(j) = (6*newTracks(i).bkgAmp(j-1)+newTracks(i).bkgAmp(j+6))/7;
            newTracks(i).xCoord(j+1) = (5*newTracks(i).xCoord(j-1)+2*newTracks(i).xCoord(j+6))/7;
            newTracks(i).yCoord(j+1) = (5*newTracks(i).yCoord(j-1)+2*newTracks(i).yCoord(j+6))/7;
            newTracks(i).amp(j+1) = (5*newTracks(i).amp(j-1)+2*newTracks(i).amp(j+6))/7;
            newTracks(i).bkgAmp(j+1) = (5*newTracks(i).bkgAmp(j-1)+2*newTracks(i).bkgAmp(j+6))/7;
            newTracks(i).xCoord(j+2) = (4*newTracks(i).xCoord(j-1)+3*newTracks(i).xCoord(j+6))/7;
            newTracks(i).yCoord(j+2) = (4*newTracks(i).yCoord(j-1)+3*newTracks(i).yCoord(j+6))/7;
            newTracks(i).amp(j+2) = (4*newTracks(i).amp(j-1)+3*newTracks(i).amp(j+6))/7;
            newTracks(i).bkgAmp(j+2) = (4*newTracks(i).bkgAmp(j-1)+3*newTracks(i).bkgAmp(j+6))/7;
            newTracks(i).xCoord(j+3) = (3*newTracks(i).xCoord(j-1)+4*newTracks(i).xCoord(j+6))/7;
            newTracks(i).yCoord(j+3) = (3*newTracks(i).yCoord(j-1)+4*newTracks(i).yCoord(j+6))/7;
            newTracks(i).amp(j+3) = (3*newTracks(i).amp(j-1)+4*newTracks(i).amp(j+6))/7;
            newTracks(i).bkgAmp(j+3) = (3*newTracks(i).bkgAmp(j-1)+4*newTracks(i).bkgAmp(j+6))/7;
            newTracks(i).xCoord(j+4) = (2*newTracks(i).xCoord(j-1)+5*newTracks(i).xCoord(j+6))/7;
            newTracks(i).yCoord(j+4) = (2*newTracks(i).yCoord(j-1)+5*newTracks(i).yCoord(j+6))/7;
            newTracks(i).amp(j+4) = (2*newTracks(i).amp(j-1)+5*newTracks(i).amp(j+6))/7;
            newTracks(i).bkgAmp(j+4) = (2*newTracks(i).bkgAmp(j-1)+5*newTracks(i).bkgAmp(j+6))/7;
            newTracks(i).xCoord(j+5) = (newTracks(i).xCoord(j-1)+6*newTracks(i).xCoord(j+6))/7;
            newTracks(i).yCoord(j+5) = (newTracks(i).yCoord(j-1)+6*newTracks(i).yCoord(j+6))/7;
            newTracks(i).amp(j+5) = (newTracks(i).amp(j-1)+6*newTracks(i).amp(j+6))/7;
            newTracks(i).bkgAmp(j+5) = (newTracks(i).bkgAmp(j-1)+6*newTracks(i).bkgAmp(j+6))/7;
        elseif j<nFrames-4 && newTracks(i).presence(j) && newTracks(i).presence(j+1) && newTracks(i).presence(j+2) && newTracks(i).presence(j+3) ...
                && newTracks(i).presence(j+4) && isnan(newTracks(i).xCoord(j)) && isnan(newTracks(i).xCoord(j+1)) ...
                && isnan(newTracks(i).xCoord(j+2)) && isnan(newTracks(i).xCoord(j+3)) && isnan(newTracks(i).xCoord(j+4))
            newTracks(i).xCoord(j) = (5*newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+5))/6;
            newTracks(i).yCoord(j) = (5*newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+5))/6;
            newTracks(i).amp(j) = (5*newTracks(i).amp(j-1)+newTracks(i).amp(j+5))/6;
            newTracks(i).bkgAmp(j) = (5*newTracks(i).bkgAmp(j-1)+newTracks(i).bkgAmp(j+5))/6;
            newTracks(i).xCoord(j+1) = (4*newTracks(i).xCoord(j-1)+2*newTracks(i).xCoord(j+5))/6;
            newTracks(i).yCoord(j+1) = (4*newTracks(i).yCoord(j-1)+2*newTracks(i).yCoord(j+5))/6;
            newTracks(i).amp(j+1) = (4*newTracks(i).amp(j-1)+2*newTracks(i).amp(j+5))/6;
            newTracks(i).bkgAmp(j+1) = (4*newTracks(i).bkgAmp(j-1)+2*newTracks(i).bkgAmp(j+5))/6;
            newTracks(i).xCoord(j+2) = (3*newTracks(i).xCoord(j-1)+3*newTracks(i).xCoord(j+5))/6;
            newTracks(i).yCoord(j+2) = (3*newTracks(i).yCoord(j-1)+3*newTracks(i).yCoord(j+5))/6;
            newTracks(i).amp(j+2) = (3*newTracks(i).amp(j-1)+3*newTracks(i).amp(j+5))/6;
            newTracks(i).bkgAmp(j+2) = (3*newTracks(i).bkgAmp(j-1)+3*newTracks(i).bkgAmp(j+5))/6;
            newTracks(i).xCoord(j+3) = (2*newTracks(i).xCoord(j-1)+4*newTracks(i).xCoord(j+5))/6;
            newTracks(i).yCoord(j+3) = (2*newTracks(i).yCoord(j-1)+4*newTracks(i).yCoord(j+5))/6;
            newTracks(i).amp(j+3) = (2*newTracks(i).amp(j-1)+4*newTracks(i).amp(j+5))/6;
            newTracks(i).bkgAmp(j+3) = (2*newTracks(i).bkgAmp(j-1)+4*newTracks(i).bkgAmp(j+5))/6;
            newTracks(i).xCoord(j+4) = (newTracks(i).xCoord(j-1)+5*newTracks(i).xCoord(j+5))/6;
            newTracks(i).yCoord(j+4) = (newTracks(i).yCoord(j-1)+5*newTracks(i).yCoord(j+5))/6;
            newTracks(i).amp(j+4) = (newTracks(i).amp(j-1)+5*newTracks(i).amp(j+5))/6;
            newTracks(i).bkgAmp(j+4) = (newTracks(i).bkgAmp(j-1)+5*newTracks(i).bkgAmp(j+5))/6;
        elseif j<nFrames-3 && newTracks(i).presence(j) && newTracks(i).presence(j+1) && newTracks(i).presence(j+2) && newTracks(i).presence(j+3) ...
                && isnan(newTracks(i).xCoord(j)) && isnan(newTracks(i).xCoord(j+1)) && isnan(newTracks(i).xCoord(j+2)) && isnan(newTracks(i).xCoord(j+3))
            newTracks(i).xCoord(j) = (4*newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+4))/5;
            newTracks(i).yCoord(j) = (4*newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+4))/5;
            newTracks(i).amp(j) = (4*newTracks(i).amp(j-1)+newTracks(i).amp(j+4))/5;
            newTracks(i).bkgAmp(j) = (4*newTracks(i).bkgAmp(j-1)+newTracks(i).bkgAmp(j+4))/5;
            newTracks(i).xCoord(j+1) = (3*newTracks(i).xCoord(j-1)+2*newTracks(i).xCoord(j+4))/5;
            newTracks(i).yCoord(j+1) = (3*newTracks(i).yCoord(j-1)+2*newTracks(i).yCoord(j+4))/5;
            newTracks(i).amp(j+1) = (3*newTracks(i).amp(j-1)+2*newTracks(i).amp(j+4))/5;
            newTracks(i).bkgAmp(j+1) = (3*newTracks(i).bkgAmp(j-1)+2*newTracks(i).bkgAmp(j+4))/5;
            newTracks(i).xCoord(j+2) = (2*newTracks(i).xCoord(j-1)+3*newTracks(i).xCoord(j+4))/5;
            newTracks(i).yCoord(j+2) = (2*newTracks(i).yCoord(j-1)+3*newTracks(i).yCoord(j+4))/5;
            newTracks(i).amp(j+2) = (2*newTracks(i).amp(j-1)+3*newTracks(i).amp(j+4))/5;
            newTracks(i).bkgAmp(j+2) = (2*newTracks(i).bkgAmp(j-1)+3*newTracks(i).bkgAmp(j+4))/5;
            newTracks(i).xCoord(j+3) = (newTracks(i).xCoord(j-1)+4*newTracks(i).xCoord(j+4))/5;
            newTracks(i).yCoord(j+3) = (newTracks(i).yCoord(j-1)+4*newTracks(i).yCoord(j+4))/5;
            newTracks(i).amp(j+3) = (newTracks(i).amp(j-1)+4*newTracks(i).amp(j+4))/5;
            newTracks(i).bkgAmp(j+3) = (newTracks(i).bkgAmp(j-1)+4*newTracks(i).bkgAmp(j+4))/5;
        elseif j<nFrames-2 &&newTracks(i).presence(j) && newTracks(i).presence(j+1) && newTracks(i).presence(j+2) ...
                && isnan(newTracks(i).xCoord(j)) && isnan(newTracks(i).xCoord(j+1)) && isnan(newTracks(i).xCoord(j+2))
            newTracks(i).xCoord(j) = (3*newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+3))/4;
            newTracks(i).yCoord(j) = (3*newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+3))/4;
            newTracks(i).amp(j) = (3*newTracks(i).amp(j-1)+newTracks(i).amp(j+3))/4;
            newTracks(i).bkgAmp(j) = (3*newTracks(i).bkgAmp(j-1)+newTracks(i).bkgAmp(j+3))/4;
            newTracks(i).xCoord(j+1) = (2*newTracks(i).xCoord(j-1)+2*newTracks(i).xCoord(j+3))/4;
            newTracks(i).yCoord(j+1) = (2*newTracks(i).yCoord(j-1)+2*newTracks(i).yCoord(j+3))/4;
            newTracks(i).amp(j+1) = (2*newTracks(i).amp(j-1)+2*newTracks(i).amp(j+3))/4;
            newTracks(i).bkgAmp(j+1) = (2*newTracks(i).bkgAmp(j-1)+2*newTracks(i).bkgAmp(j+3))/4;
            newTracks(i).xCoord(j+2) = (newTracks(i).xCoord(j-1)+3*newTracks(i).xCoord(j+3))/4;
            newTracks(i).yCoord(j+2) = (newTracks(i).yCoord(j-1)+3*newTracks(i).yCoord(j+3))/4;
            newTracks(i).amp(j+2) = (newTracks(i).amp(j-1)+3*newTracks(i).amp(j+3))/4;
            newTracks(i).bkgAmp(j+2) = (newTracks(i).bkgAmp(j-1)+3*newTracks(i).bkgAmp(j+3))/4;
        elseif j<nFrames-1 &&newTracks(i).presence(j) && newTracks(i).presence(j+1) && isnan(newTracks(i).xCoord(j)) && isnan(newTracks(i).xCoord(j+1))
            newTracks(i).xCoord(j) = (2*newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+2))/3;
            newTracks(i).yCoord(j) = (2*newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+2))/3;
            newTracks(i).amp(j) = (2*newTracks(i).amp(j-1)+newTracks(i).amp(j+2))/3;
            newTracks(i).bkgAmp(j) = (2*newTracks(i).bkgAmp(j-1)+newTracks(i).bkgAmp(j+2))/3;
            newTracks(i).xCoord(j+1) = (newTracks(i).xCoord(j-1)+2*newTracks(i).xCoord(j+2))/3;
            newTracks(i).yCoord(j+1) = (newTracks(i).yCoord(j-1)+2*newTracks(i).yCoord(j+2))/3;
            newTracks(i).amp(j+1) = (newTracks(i).amp(j-1)+2*newTracks(i).amp(j+2))/3;
            newTracks(i).bkgAmp(j+1) = (newTracks(i).bkgAmp(j-1)+2*newTracks(i).bkgAmp(j+2))/3;
        elseif newTracks(i).presence(j) && isnan(newTracks(i).xCoord(j))
            newTracks(i).xCoord(j) = (newTracks(i).xCoord(j-1)+newTracks(i).xCoord(j+1))/2;
            newTracks(i).yCoord(j) = (newTracks(i).yCoord(j-1)+newTracks(i).yCoord(j+1))/2;
            newTracks(i).amp(j) = (newTracks(i).amp(j-1)+newTracks(i).amp(j+1))/2;
            newTracks(i).bkgAmp(j) = (newTracks(i).bkgAmp(j-1)+newTracks(i).bkgAmp(j+1))/2;
        end
    end
end
end
