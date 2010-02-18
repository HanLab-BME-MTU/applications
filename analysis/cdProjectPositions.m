function [n_spindleVector, cenPosNorm, goodTimes, s1s2c1c2int] = cdProjectPositions(idlist,tagCorrection,randomAssignment)
%CDPROJECTPOSITIOS projects the cen positions of an idlist onto the spb axis
%
% n_spindleVector : spindle length in microns
% cenPosNorm : relative cen positions projected onto axis
% goodTimes : frames where idlist.linklist is not empty
% s1s2c1c2int : spot intensities

if nargin == 1 || isempty(tagCorrection)
    tagCorrection = 0;
end
if nargin < 3 || isempty(randomAssignment)
    randomAssignment = false;
end

% find indices of spb, cen
spb1idx = find(strncmpi('spb1',idlist(1).stats.labelcolor,4)); % allow sbp1*
spb2idx = find(strcmpi(idlist(1).stats.labelcolor,'spb2'));
cen1idx = find(strcmpi(idlist(1).stats.labelcolor,'cen1'));
cen2idx = find(strcmpi(idlist(1).stats.labelcolor,'cen2'));
if isempty(cen2idx)
    if isempty(cen1idx)
        cen1idx = find(strcmpi(idlist(1).stats.labelcolor,'cen1*'));
    else
        cen2idx = cen1idx;
    end
end



% loop through idlist, find spb-axis, spb-cen vectors
tmax = length(idlist);
spindleVector = zeros(tmax,3);
s1c1Vector = zeros(tmax,3);
s2c2Vector = NaN(tmax,3);
s1s2c1c2int = NaN(tmax,4);
goodTimeL = false(tmax,1);

for t = 1:tmax
    if ~isempty(idlist(t).linklist)

        % find vector in direction of spindle
        spindleVector(t,:) =...
            diff(idlist(t).linklist([spb1idx, spb2idx],9:11));

        % spb - cen vectors
        s1c1Vector(t,:) = ...
            diff(idlist(t).linklist([spb1idx, cen1idx],9:11));
        if ~isempty(cen2idx)
        s2c2Vector(t,:) = ...
            diff(idlist(t).linklist([spb2idx, cen2idx],9:11));
        end

        % store intensities
        if ~isempty(cen2idx)
        s1s2c1c2int(t,:) = ...
            idlist(t).linklist([spb1idx, spb2idx, cen1idx, cen2idx], 8)';
        else
            s1s2c1c2int(t,1:3) = ...
            idlist(t).linklist([spb1idx, spb2idx, cen1idx], 8)';
        end

        % remember goodTime
        goodTimeL(t) = true;

    else
        % do nothing
    end
end % for t = 1:tmax

% shorten Vectors
spindleVector(~goodTimeL,:) = [];
s1c1Vector(~goodTimeL,:) = [];
s2c2Vector(~goodTimeL,:) = [];
s1s2c1c2int(~goodTimeL,:) = [];

% normalize spindleVector
[n_spindleVector, e_spindleVector] = normList(spindleVector);

% In case we want to correct for the centromere, we need to have the
% sxcx vectors normed, too
[n_s1c1Vector, e_s1c1Vector] = normList(s1c1Vector);
[n_s2c2Vector, e_s2c2Vector] = normList(s2c2Vector);

% correct - but don't make the length negative!
s1c1Vector = s1c1Vector - repmat(min(tagCorrection,n_s1c1Vector-0.08),1,3) .* e_s1c1Vector;
s2c2Vector = s2c2Vector - repmat(min(tagCorrection,n_s2c2Vector-0.08),1,3) .* e_s2c2Vector;


% project spb-cen vectors. Distance from spb1
cen1Dist = dot(s1c1Vector, e_spindleVector, 2);
cen2Dist = n_spindleVector + dot(s2c2Vector, e_spindleVector, 2);

% get normalized cen positions.
cenPosNorm = [cen1Dist./n_spindleVector, cen2Dist./n_spindleVector];

goodTimes = find(goodTimeL);


% potentially: make lower spb into spb1
if randomAssignment
    lowerSpb2 = idlist(goodTimes(1)).linklist(spb1idx,11) > ...
        idlist(goodTimes(1)).linklist(spb1idx,11);
    if lowerSpb2
        cenPosNorm = cenPosNorm(:,[2,1]);
    end
end


