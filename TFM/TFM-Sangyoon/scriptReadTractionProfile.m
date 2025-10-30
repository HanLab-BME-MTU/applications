% make the channel
tChan = Channel('./Cell');
numFrames = numel(tChan.getImageFileNames);
tStack = tChan.loadStack(1:numFrames);
%% Pick the points
figure, imshow(tStack(:,:,1), [0 700]); colormap jet
[X,Y] = ginput; close
figure,  imshow(tStack(:,:,1), [0 700]); hold on
myColors = distinguishable_colors(numel(X),'b');
%% Read from the stack using X and Y
for ii=1:numel(X)
    plot(X(ii),Y(ii),'*','Color',myColors(ii,:))
    tProfiles(:,ii) = tStack(round(Y(ii)),round(X(ii)),:);
end
%% Plot tProfiles
figure, hold on
tInt = 5/60; %in min
myX = (1:numFrames)*tInt; % min
for ii=1:numel(X)
    plot(myX, tProfiles(:,ii),'.-','Color',myColors(ii,:))
end
xlabel('Time (min)')
ylabel('Traction mag (Pa)')
