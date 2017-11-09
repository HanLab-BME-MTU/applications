function pointSourceSpotsPlot(MD, movieinfo)
%EB3SpotsPlot plots max projection image and detected spots from EB3
%   Detailed explanation goes here

vol = double(MD.getChannel(1).loadStack(1));
figure, imshow(max(vol,[],3),[])
hold on
for i = 1:(numel(movieinfo.xCoord)/2)
    plot(movieinfo.xCoord(i), movieinfo.yCoord(i), 'r+');
end

end

