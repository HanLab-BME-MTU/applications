function subsamplePoints(obj,fraction)

nPointsBefore = obj.data.nPoints;

[~,frameIdx] = sort(obj.data.frame);

if fraction >= 0 && fraction < 1
    obj.data.points = obj.data.points(frameIdx(1:round(fraction*obj.data.nPoints)),:);
elseif fraction > 1
    obj.data.points = obj.data.points(frameIdx(1:fraction),:);
elseif fraction == 1
    disp('Process: fraction == 1, data will not be subsampled!')
else
    disp('Process: This subsample fraction is not valid!')
end

nPointsAfter = obj.data.nPoints;

fprintf('Process: Data points subsampled: %d -> %d\n',nPointsBefore,nPointsAfter);

end


