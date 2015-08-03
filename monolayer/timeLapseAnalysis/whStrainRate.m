function [] = whStrainRate(params,dirs)

time = 1 : params.nTime - params.frameJump;% - 1; % -1 because segmentation is based on motion

for t = time   
    strainRateFname = [dirs.strainRate pad(t,3) '_strainRate.mat'];
    
    if exist(strainRateFname,'file') && ~params.always
        continue;
    end
    
    mfFname = [dirs.mfData pad(t,3) '_mf.mat']; % dxs, dys
        
    load(mfFname);
    
    [sizeY,sizeX] = size(dxs);
    [dxs] = imresize(dxs,1./params.patchSize);
    [dys] = imresize(dys,1./params.patchSize);
    
    speed = sqrt(dxs.^2 + dys.^2);             
    
    [directionY,directionX] = quantizeAngles(dys,dxs);
    strainRate = getStrainRate(directionX,directionY,speed);
    sr = imresize(strainRate,[sizeY,sizeX],'nearest');
    
    save(strainRateFname,'sr');
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
function [strainDirectionsYs,strainDirectionsXs] = quantizeAngles(dys,dxs)
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

strainDirectionsXs = nan(size(angles));
strainDirectionsYs = nan(size(angles));

MASK = (positions == 1); % 0
strainDirectionsYs(MASK) = 0;
strainDirectionsXs(MASK & dxs >= 0) = 1;
strainDirectionsXs(MASK & dxs < 0) = -1;

MASK = (positions == 2); % 90
strainDirectionsXs(MASK & dxs >= 0) = 1;
strainDirectionsXs(MASK & dxs < 0) = -1;
strainDirectionsYs(MASK & dys >= 0) = 1;
strainDirectionsYs(MASK & dys < 0) = -1;

MASK = (positions == 3); % 90
strainDirectionsXs(MASK) = 0;
strainDirectionsYs(MASK & dys >= 0) = 1;
strainDirectionsYs(MASK & dys < 0) = -1;

% MASK = (positions == 1); % -90
% strainDirectionsXs(MASK) = 0;
% strainDirectionsYs(MASK) = -1;
% 
% MASK = (positions == 2); % -45
% strainDirectionsXs(MASK & dxs >= 0) = 1;
% strainDirectionsYs(MASK & dxs >= 0) = -1;
% strainDirectionsXs(MASK & dxs < 0) = -1;
% strainDirectionsYs(MASK & dxs < 0) = 1;
% 
% MASK = (positions == 3); % 0
% strainDirectionsYs(MASK) = 0;
% strainDirectionsXs(MASK & dxs >= 0) = 1;
% strainDirectionsXs(MASK & dxs < 0) = -1;
% 
% MASK = (positions == 4); % 45
% strainDirectionsXs(MASK & dxs >= 0) = 1;
% strainDirectionsYs(MASK & dxs >= 0) = 1;
% strainDirectionsXs(MASK & dxs < 0) = -1;
% strainDirectionsYs(MASK & dxs < 0) = -1;
% 
% MASK = (positions == 5); % 90
% strainDirectionsXs(MASK) = 0;
% strainDirectionsYs(MASK) = 1;
end

function [strainRate] = getStrainRate(directionX,directionY,speed)
[sizeY, sizeX] = size(speed);
[XGrid,YGrid] = meshgrid(1:size(speed,2),1:size(speed,1));

XGridBefore = XGrid -  directionX;
XGridAfter = XGrid +  directionX;
YGridBefore = YGrid -  directionY;
YGridAfter = YGrid +  directionY;
MASK = ...
    XGridBefore > 0 & XGridBefore < sizeX &...
    XGridAfter > 0 & XGridAfter < sizeX &...
    YGridBefore > 0 & YGridBefore < sizeY &...
    YGridAfter > 0 & YGridAfter < sizeY;

strainRate = nan(size(speed));
for y = 1 : sizeY
    for x = 1 : sizeX
        if MASK(y,x)
            strainRate(y,x) = speed(YGridAfter(y,x),XGridAfter(y,x)) - speed(YGridBefore(y,x),XGridBefore(y,x));
        end
    end
end

end



