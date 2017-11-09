function [ h, x, y, u, v] = drawOrientation( theta , arrow_scale)
%drawOrientation Draw quivers from normal theta encoded in matrix
    
    if(nargin < 2)
        arrow_scale = 5;
    end

    idx =  find(~isnan(theta));
    % grab a random permutation of arrows, reduce the value of arrows
%     idx = idx(randperm(length(idx),floor(length(idx)/arrows)));
    % grab the subscripts for the coordinates
    [r,c] = ind2sub(size(theta),idx);
    % xy are reversed in ij space
    % theta is the normal vector
    x = c;
    y = r;
    u = sin(theta(idx));
    v = -cos(theta(idx));
    h = quiver(c,r,sin(theta(idx)).*arrow_scale,-cos(theta(idx)).*arrow_scale,0,'color','white');

end

