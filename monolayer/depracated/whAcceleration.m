function [] = whAcceleration(params,dirs)

time = 2 : params.nTime - params.frameJump - 1; % -1 because segmentation is based on motion

for t = time   
    accelerationFname = [dirs.acceleration pad(t,3) '_acceleration.mat'];
    
    if exist(accelerationFname,'file') && ~params.always
        continue;
    end
    
    mfFname = [dirs.mfData pad(t,3) '_mf.mat']; % dxs, dys
        
    load(mfFname);
    
    [sizeY,sizeX] = size(dxs);
    [dxs] = imresize(dxs,1./params.patchSize);
    [dys] = imresize(dys,1./params.patchSize);      
    
    [directionY,directionX] = quantizeAngles(dys,dxs);
    
    clear dxs; clear dys;
    load([dirs.mfData pad(t-1,3) '_mf.mat']);
    [dxs] = imresize(dxs,1./params.patchSize);
    [dys] = imresize(dys,1./params.patchSize);
    speedBefore = sqrt(dxs.^2 + dys.^2);
    
    clear dxs; clear dys;
    load([dirs.mfData pad(t+1,3) '_mf.mat']);
    [dxs] = imresize(dxs,1./params.patchSize);
    [dys] = imresize(dys,1./params.patchSize);
    speedAfter = sqrt(dxs.^2 + dys.^2);
    clear dxs; clear dys;
    
    acc = getAcceleration(directionX,directionY,speedBefore,speedAfter);
    acc = imresize(acc,[sizeY,sizeX],'nearest');        
    
    save(accelerationFname,'acc');
end

end

%% get map to calculate strain rate
% xs = [1,  1,  0, -1,  -1,  -1,   0,   1];
% ys = [0,  1,  1,  1,   0,  -1,  -1,  -1];
%       0   45  90 -45   0    45  -90  -45
% >> atand(ys./xs)
% 
% ans =
% 
%      0    45    90   -45     0    45   -90   -45
function [directionsYs,directionsXs] = quantizeAngles(dys,dxs)
angles = atand(abs(dys./dxs)); 
quantAngles = 0:45:90;
nQuantAngles = length(quantAngles);

% Find closest angles
diff = nan(size(angles,1),size(angles,2),nQuantAngles);
for i = 1 : nQuantAngles
    diff(:,:,i) = abs(quantAngles(i) - angles);
end

[match,positions] = min(diff,[],3);

% % this is the closest quantified angle
% anglesQnatized = nan(size(angles));
% for i = 1 : nQuantAngles
%     anglesQnatized(positions == i) = quantAngles(i);    
% end

directionsXs = nan(size(angles));
directionsYs = nan(size(angles));

MASK = (positions == 1); % 0
directionsYs(MASK) = 0;
directionsXs(MASK & dxs >= 0) = 1;
directionsXs(MASK & dxs < 0) = -1;

MASK = (positions == 2); % 90
directionsXs(MASK & dxs >= 0) = 1;
directionsXs(MASK & dxs < 0) = -1;
directionsYs(MASK & dys >= 0) = 1;
directionsYs(MASK & dys < 0) = -1;

MASK = (positions == 3); % 90
directionsXs(MASK) = 0;
directionsYs(MASK & dys >= 0) = 1;
directionsYs(MASK & dys < 0) = -1;

% MASK = (positions == 1); % -90
% directionsXs(MASK) = 0;
% directionsYs(MASK) = -1;
% 
% MASK = (positions == 2); % -45
% directionsXs(MASK & dxs >= 0) = 1;
% directionsYs(MASK & dxs >= 0) = -1;
% directionsXs(MASK & dxs < 0) = -1;
% directionsYs(MASK & dxs < 0) = 1;
% 
% MASK = (positions == 3); % 0
% directionsYs(MASK) = 0;
% directionsXs(MASK & dxs >= 0) = 1;
% directionsXs(MASK & dxs < 0) = -1;
% 
% MASK = (positions == 4); % 45
% directionsXs(MASK & dxs >= 0) = 1;
% directionsYs(MASK & dxs >= 0) = 1;
% directionsXs(MASK & dxs < 0) = -1;
% directionsYs(MASK & dxs < 0) = -1;
% 
% MASK = (positions == 5); % 90
% directionsXs(MASK) = 0;
% directionsYs(MASK) = 1;
end

function [acceleration] = getAcceleration(directionX,directionY,speedBefore,speedAfter)
[sizeY, sizeX] = size(speedBefore);
[XGrid,YGrid] = meshgrid(1:size(speedBefore,2),1:size(speedBefore,1));

XGridBefore = XGrid -  directionX;
XGridAfter = XGrid +  directionX;
YGridBefore = YGrid -  directionY;
YGridAfter = YGrid +  directionY;
MASK = ...
    XGridBefore > 0 & XGridBefore < sizeX &...
    XGridAfter > 0 & XGridAfter < sizeX &...
    YGridBefore > 0 & YGridBefore < sizeY &...
    YGridAfter > 0 & YGridAfter < sizeY;

acceleration = nan(size(speedBefore));
for y = 1 : sizeY
    for x = 1 : sizeX
        if MASK(y,x)
            acceleration(y,x) = speedAfter(YGridAfter(y,x),XGridAfter(y,x)) - speedBefore(YGridBefore(y,x),XGridBefore(y,x));
        end
    end
end

end





