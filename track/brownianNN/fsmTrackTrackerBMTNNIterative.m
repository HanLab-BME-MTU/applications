function M=fsmTrackTrackerBMTNNIterative(initM,I,J,threshold,influence)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If an initializer for the tracker exists, use it to propagate speckle positions at frame I
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(initM)
    
    % Sort I and make sure there are no empty lines
    I=sortrows(I(find(I(:,1)~=0),:));
    
    % Propagate speckle positions
    spPos=I(:,1:2); % Extract only positions
    pSpPos=fsmTrackPropSpecklePos(spPos,initM,'FORWARD');
    
    % Now I is the propagated version of the original I
    I=cat(2,pSpPos,I(:,3)); % Add original intensities to the propagated positions
    clear Itmp;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Track with the brownian motion tracker with neural network
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mone=fsmTrackTrackerBMTNN(I,J,threshold,influence);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If needed, correct back the propagated coordinates of the first frame
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(initM)

    % Create a copy of Mone
    copyMone=Mone;
    
    % Correct back the coordinates of the first frame (= currentM(:,1:2))
    for i=1:size(spPos,1)
        indx=find(Mone(:,1)==pSpPos(i,1) & Mone(:,2)==pSpPos(i,2));
        if length(indx)~=1
            error('Propagated position not univocally found in Mone');
        end
        copyMone(indx,1:2)=spPos(i,:);
    end
    M=copyMone;

else
    
    % Return the result of the first tracking
    M=Mone;
    
end
