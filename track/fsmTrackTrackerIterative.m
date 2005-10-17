function M=fsmTrackTrackerIterative(initM,I,J,threshold,influence, TRACKER, initCorLen)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If an initializer for the tracker exists, use it to propagate speckle positions at frame I
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7
    initCorLen = Inf;
end

if ~isempty(initM)
    
    % Sort I and make sure there are no empty lines
    I=sortrows(I(find(I(:,1)~=0),:));
    
    % Propagate speckle positions
    spPos=I(:,1:2); % Extract only positions
    pSpPos=fsmTrackPropSpecklePos(spPos,initM,'FORWARD',initCorLen);
    
    % Now I is the propagated version of the original I
    I=cat(2,pSpPos,I(:,3)); % Add original intensities to the propagated positions
    clear Itmp;
    
end

if TRACKER == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Track with the brownian motion tracker with neural network
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Mone=fsmTrackTrackerBMTNN(I,J,threshold,influence);
elseif TRACKER == 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Track with linear assignment code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NONLINK_MARKER = -1;    % value to indicate that two points cannot be linked.
    extendedTesting = 0;
    augmentCC = 1;          % automatically generate ghost marker to account for birth and deaths
    

    % Extract speckle coordinates
    posI=I(:,1:2);
    posJ=J(:,1:2);
    % create cost matrix for the linear assignment code
    % in the most general case that is equal to the distance matrix 
    cc = createSparseDistanceMatrix(posI, posJ, threshold);
    [all_links_1, all_links_2] = lap(cc, NONLINK_MARKER, extendedTesting, augmentCC);
    
    I_N = length(posI);
    J_N = length(posJ);
    
    % first get points that are in frame 1 and have a link in frame 2
    all_links_1_2 = all_links_1(1:I_N);
    frame_1_indices  = find(all_links_1_2 <= J_N);
    frame_2_indices  = all_links_1(frame_1_indices);
    Mone = [posI(frame_1_indices,1:2),  posJ(frame_2_indices,1:2)];
    
    % Add non-paired speckles from frame 1
    frame_1_indices_noLink  = find(all_links_1_2 > J_N);
    lenNPI=size(frame_1_indices_noLink,1);
    MNPI=zeros(lenNPI,4);
    MNPI(1:lenNPI,1:2)=posI(frame_1_indices_noLink,1:2);
    Mone=cat(1,Mone,MNPI);    

    % Add non-paired speckles from frame 2
    all_links_2_1 = all_links_2(1:J_N);
    frame_2_indices_noLink  = find(all_links_2_1 > I_N);
    lenNPJ=size(frame_2_indices_noLink,1);
    MNPJ=zeros(lenNPJ,4);
    MNPJ(1:lenNPJ,3:4)=posJ(frame_2_indices_noLink,1:2);
    Mone=cat(1,Mone,MNPJ);    
end

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
