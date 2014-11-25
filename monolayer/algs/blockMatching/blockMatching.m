% im1, im2 - images
% width - block's width
% searchR - how far to look for a match
function [dydx, dys, dxs, scores] = blockMatching(im1, im2, width, searchR, ROI)

if nargin < 5
    ROI = true(size(im1));
end

im1 = double(im1);
im2 = double(im2);

% initialization
[r1, c1] = size(im1);
[r2, c2] = size(im2);
if r1 ~= r2 || c1 ~= c2
    error('The images must be in a same size')
end

yBins = int16(floor((r1-2*searchR) ./ width));
xBins = int16(floor((c1-2*searchR) ./ width));

dxs = nan(r1,c1);
dys = nan(r1,c1);
scores = nan(r1,c1);
countLocalMatches = 0;
for y = 1 : yBins
    for x = 1 : xBins
        sy = (y-1) * width + searchR + 1;
        fy = y*width + searchR;
        sx = (x-1) * width + searchR + 1;
        fx = x*width + searchR;
        if sum(sum(ROI(sy:fy,sx:fx))) < 0.6*width*width
            continue;
        end
        bb = im1(sy:fy,sx:fx);
        bb = bb ./ sum(bb(:));
        %         fprintf('search %d,%d->%d,%d\n\n',sy,sx,fy,fx);
        [dy, dx, score] = localMatch(bb,im2,width,searchR,floor((sx+fx)/2),floor((sy+fy)/2));
        dys(sy:fy,sx:fx) = dy;
        dxs(sy:fy,sx:fx) = dx;
        scores(sy:fy,sx:fx) = score;
        countLocalMatches = countLocalMatches + 1;
    end
end

precentLocalMatches = countLocalMatches /  double(double(yBins) * double(xBins));

dydx = [mean(dys(:)),mean(dxs(:))];

end

function [dy, dx, score] = localMatch(bb,im2,width,searchR,centerX,centerY)
bsize = 2*searchR + 1;
matchScores = zeros(bsize,bsize);
h = floor(width/2);
initY = centerY - h - searchR;
initX = centerX - h - searchR;

for y = 1 : bsize
    for x = 1 : bsize
        sy = initY + (y-1);
        fy = sy + width - 1;
        sx = initX + (x-1);
        fx = sx + width - 1;
        %         fprintf('%d,%d->%d,%d\n',sy,sx,fy,fx);
        bb2 = im2(sy:fy,sx:fx);
        bb2 = bb2 ./ sum(bb2(:));
        corr = sqrt(bb .* bb2); % Bhattacharyya coefficient
        matchScores(y,x) = sum(corr(:));
    end
end

% scores = matchScores(:);
% if (max(scores) > mean(scores) + std(scores))
score = 0;
dx = -1;
dy = -1;
for i = 1 : bsize
    for j = 1 : bsize
        if (score < matchScores(i,j))
            score = matchScores(i,j);
            dy = i - searchR - 1;
            dx = j - searchR - 1;
        end
    end
end
% else
%     dx = 0;
%     dy = 0;
% end
% fprintf('------------------\n');

end