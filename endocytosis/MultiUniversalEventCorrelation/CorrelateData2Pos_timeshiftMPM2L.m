function [MPMshift] = CorrelateData2Pos_timeshiftMPM2L(MPMbasis, tvec);
% CorrelateData2Pos_timeshiftMPM creates a three-dimensional MPM file,
% where time-shifted version of the basis MPM are contained in the
% subsequent layers
%
% SYNOPSIS [results] = CorrelateData2Pos(positions, data, rdist, type)
%
% INPUT     MPMbasis:   mpm file containing positions of tracked features 
%                       in x,y-columns
%           tvec:       time vector (frames)  
%
% OUTPUT:   MPMshift:        
%
% last modified: Dinah Loerke, Aug 28, 2008
% 
% NOTE: This function only creates a simple time shift, i.e. the identical
% positions are shifted to earlier and later time points. If tracking of
% points in a trajectory is desired - e.g. you want to examine an earlier
% point in a moving trajectory, as opposed to the fixed final location at
% an earlier time point, then use the function
% CorrelateData2Pos_expandMPMtime

tlen = length(tvec);
[sx,sy] = size(MPMbasis);
MPMshift = MPMbasis;
MPMtvec = nan*MPMshift;

for i=1:tlen
        
    % current time shift
    tau = tvec(i);
    
    cmpm = MPMbasis;
    lmpm = MPMbasis;
    
    if tau<0
        lmpm(:) = nan;
        cmpm(:,1:2*abs(tau)) = [];
        [csx,csy] = size(cmpm);
        lmpm(:,1:csy) = cmpm;
              
    elseif tau>0
        lmpm(:) = nan;
        cmpm(:,sy-2*tau+1:sy) = [];
        [csx,csy] = size(cmpm);
        lmpm(:,sy-csy+1:sy) = cmpm;
        
    end

    upos = find(isfinite(lmpm));
    MPMshift(upos) = lmpm(upos);
    MPMtvec(upos) = i;
    
end

MPMshift(:,:,2) = MPMtvec;

end % of function
