function movieInfo = wyssFileConvertforTracking(fulldata)
% takes a PointList.fulldata from wyssFileImport and converts it into a movieInfo
% type cell array for using rackCloseGapsKalmanSparse
%
%


n = max(fulldata(:,1));
frames = fulldata(:,1);
movieInfo = repmat(struct('xCoord',[],'yCoord',[],'amp',[]),[n,1]);


for i = 1:n;
    ind = find(frames==i);
    movieInfo(i).xCoord = [fulldata(ind,[2,5])];
    movieInfo(i).yCoord = [fulldata(ind,[3,6])];
    movieInfo(i).amp = [fulldata(ind,[10,9])];
end


end
