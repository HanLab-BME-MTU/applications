function [pixelPath,piece]=createPixelPath(pts)
pixelPath=[];
for i = 1:size(pts,1)-1
    p1 = pts(i  ,:);
    p2 = pts(i+1,:);
    piece(i).pix = bresenham(p1,p2,4); % pList is a Mx2 matrix
    if i==1
        pixelPath=vertcat(pixelPath,piece(i).pix(1:end,:)); 
    else
        % To avoid double points:
        pixelPath=vertcat(pixelPath,piece(i).pix(2:end,:));
    end
end
pixelPath=removeDoublePoints(pixelPath);