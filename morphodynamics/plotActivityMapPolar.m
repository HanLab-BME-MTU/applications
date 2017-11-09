function handle = plotActivityMapPolar(mapIn)

%UNDER CONSTRUCTION BLA BLA BLA


mapIn(isnan(mapIn)) = min(mapIn(:));

polarplot3d(mapIn(:,end:-1:1)','TickSpacing',0,'PolarGrid',{1 1});
view(2)



% [M,N] = size(mapIn);
% 
% r = 
% 
% [X,Y] = meshgrid(-r-1:r+1,-r-1:r+1);
% 
% m = sqrt(X.^2 + Y.^2) < r;
% 
% w = getMaskWindows
% 
% 
% 
