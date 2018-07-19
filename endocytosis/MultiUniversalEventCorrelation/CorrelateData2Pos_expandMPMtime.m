function [results] = CorrelateData2Pos_expandMPMtime(positions, data, rdist);
%% under construction - don't use at this point

%function [M_ext]=extendMat(M,timewindow)
% this function expands matrices like Mx or My to extend the first and 
% last value in a trajectory to a specified number of frames before and 
% after the trajectory appearance/disappearance


M = full(M);
[sx,sy]=size(M);

h= waitbar(0,'extending matrix');
for i=1:sx
    currentTraj = M(i,:);
    extTraj = currentTraj;
    TrajFirstPoint = min(find(currentTraj>0));
    % if there's any empty frames before appearance
    if ( (~isempty(TrajFirstPoint)) & (TrajFirstPoint>1) )
        TrajFirstPointVal = currentTraj(TrajFirstPoint);
        Eleft = max(TrajFirstPoint-timewindow,1);
        extTraj(Eleft:TrajFirstPoint-1) = TrajFirstPointVal;
    end
    
    TrajLastPoint = max(find(currentTraj>0));
    % if there's any empty frames after disappearance
    if ( (~isempty(TrajLastPoint)) & (TrajLastPoint<sy) )
        TrajLastPointVal = currentTraj(TrajLastPoint);
        Eright = min(TrajLastPoint+timewindow,sy);
        extTraj(TrajLastPoint+1:Eright) =  TrajLastPointVal;
    end
   
   M_ext(i,:) = extTraj;
   
   if (mod(i,round(sx/50))==0)
        waitbar(i/sx); 
   end
    
end % of for loop
close(h);
   
end % of subfunction