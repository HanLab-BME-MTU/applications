function [ tz, edge_tz, mean_tz ] = findFocalPlane( MD , framePenalty)
%findFocalPlane Finds the focal plane(s) of a movie

if(nargin < 2)
    framePenalty = 0.95;
end

R = CellReader(CachedReader(MD.getReader()));
eC = mat2gray(edgeCritereon(R));
mC = mat2gray(meanCritereon(R));
% fC = 1:-1/(R.getSizeT*R.getSizeZ-1):0;
fC = framePenalty.^(0:(R.getSizeT*R.getSizeZ-1));
fC = repmat(fC,R.getSizeC,1);

criteria = eC.*mC.*fC;
[~,tz] = max(criteria,[],2);

[~,edge_tz] = max(eC,[],2);
[~,mean_tz] = max(mC,[],2);

if(any(tz == 1))
    if(tz == 1)
        tz(tz == 1) = mean_tz(tz == 1);
    end
end

if(any(tz == length(criteria)))
    [lower.tz, lower.edge_tz, lower.mean_tz] = lamins.functions.findFocalPlane(MD,framePenalty - 0.1);
    tz = lower.tz;
    edge_tz = lower.edge_tz;
    mean_tz = lower.mean_tz;
end



% cleanup
R.delete();

    function out = edgeCritereon(R)
        h = -fspecial('log',[25],5);
        out = cellfun(@(x) sum(joinColumns(imfilter(double(x),h))),R(:,:));
    end
    function out = meanCritereon(R)
        out = cellfun(@(x) mean(x(:)),R(:,:));
    end

end

