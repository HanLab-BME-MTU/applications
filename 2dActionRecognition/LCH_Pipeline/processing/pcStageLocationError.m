
function [] = pcStageLocationError(MD,params,dirs)

% Assaf Zaritsky, June 2017
searchR = 4;

if ~isfield(params,'sTime')
    params.sTime = 1;
end

stageLocationErrorFname = [dirs.mfData 'stageLocationError.mat'];

if exist(stageLocationErrorFname,'file') && ~params.always   
    return;
end


time = params.sTime : params.nTime - params.frameJump;

stageShift = nan(1,length(time));

for t = time    
    I0 = MD.getChannel(1).loadImage(t);
    I1 = MD.getChannel(1).loadImage(t + params.frameJump);    
    
    stageShift(t-time(1)+1) = calcShift(I0, I1, searchR);     
end
save(stageLocationErrorFname,'stageShift');
end

%%

function [shift] = calcShift(im1, im2, searchR)

maskPixels = 200;
nrep = 20;
dxs = nan(1,nrep);
dys = nan(1,nrep);
for i = 1 : nrep
    [curDx,curDy] = calcShiftRand(im1,im2,maskPixels,searchR);
    dxs(i) = curDx;
    dys(i) = curDy;
end
shift = sqrt(mean(dxs).^2 + mean(dys).^2);
end

%%
function [dx,dy] = calcShiftRand(im1,im2,maskPixels,searchR)

im1 = double(im1);
im2 = double(im2);

% initialization
[r1, c1] = size(im1);

a = searchR+1;
bx = c1 - maskPixels - searchR - 1;
by = r1 - maskPixels - searchR - 1; 
randX = floor((bx-a)*rand()+a);
randY = floor((by-a)*rand()+a);

bb1 = im1(randY:randY+maskPixels,randX:randX+maskPixels);
bb1 = bb1 ./ sum(bb1(:));

dx = nan;
dy = nan;
curScore = -inf;

for idy = -searchR : searchR
    for idx =  -searchR : searchR
        bb2 = im2((randY+idy):(randY+maskPixels+idy),(randX+idx):(randX+maskPixels+idx));
        bb2 = bb2 ./ sum(bb2(:));
        corr12 = sqrt(bb1 .* bb2); % Bhattacharyya coefficient
        corr12 = sum(corr12(:));
        if corr12 > curScore
            curScore = corr12;
            dx = idx;
            dy = idy;
        end
    end
end
end

