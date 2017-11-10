function spotsPlot1(movieName,channelNum,cordVar)

%SPOTSPLOT plots detected spots with x,y coordinates on averaged z-stack layer
%
% INPUT cordVar :   variable contains spots coordinates
%
% OUTPUT Figure with spots plot

% 07/14 Ning


spotsNum = size(cordVar.sp,2);
spotsCord = zeros(spotsNum,3);

for i = 1:spotsNum
    for j = 1:3
        spotsCord(i,j) = cordVar.sp(i).cord(j);
    end
end
avgLayer = round(mean(spotsCord(:,3)));

% Take movie data and process to get 3D matrix
Path = strcat('', movieName);
MD = MovieData.load(Path);
frameN = 1;
layerMat = MD.getChannel(channelNum).loadImage(frameN,avgLayer);

% xSize = MD.imSize_(1);
% ySize = MD.imSize_(2);
% nDepth = MD.zSize_;
% zStack3D  = zeros(xSize, ySize, nDepth);
% for ii = 1:nDepth
%     zStack3D(:,:,ii) = MD.getChannel(channelNum).loadImage(frameN,ii);
% end

figName = strcat(inputname(3), '; layer:', num2str(avgLayer), '; spotsFound:', num2str(spotsNum));
figure('name', figName);
imshow(layerMat,[]);
hold on;

for singleSpot = 1:spotsNum
    plot(spotsCord(singleSpot,1),spotsCord(singleSpot,2),'mo');
end

Y=pdist(spotsCord);
Z=linkage(Y);
T=cluster(Z,'cutoff',33,'criterion','distance');

aa = max(T);
realcord = zeros(aa,3);
for i = 1:aa
    temp = zeros(size(T,1),3);
    for j = 1:size(T,1)
        if T(j) == i
            temp(j,:) = spotsCord(j,:);
        end
    end
    temp(~any(temp,2),:) = [];
    row = size(temp,1);
    if row == 1
        realcord(i,:) = temp;
    else
        realcord(i,:) = mean(temp);
    end
    hold on
    plot(realcord(i,1),realcord(i,2),'c+')
end


hold off;
figure;
dendrogram(Z);