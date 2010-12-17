function checkFlow(flow,nbins)
subplot(1,3,1)
quiver(flow(:,2),flow(:,1),flow(:,4)-flow(:,2),flow(:,3)-flow(:,1))
set(gca,'Ydir','reverse')
title('Flow')

subplot(1,3,2)
hist(flow(:,4)-flow(:,2),nbins)
title('Histogram of the x-component of the Flow')

subplot(1,3,3)
hist(flow(:,3)-flow(:,1),nbins)
title('Histogram of the y-component of the Flow')