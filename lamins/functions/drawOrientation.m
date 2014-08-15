function [ h ] = drawOrientation( theta )
%drawOrientation Draw quivers from normal theta encoded in matrix

    idx =  find(~isnan(theta));
    % grab a random permutation of arrows, reduce the value of arrows
%     idx = idx(randperm(length(idx),floor(length(idx)/arrows)));
    % grab the subscripts for the coordinates
    [r,c] = ind2sub(size(theta),idx);
    % xy are reversed in ij space
    % theta is the normal vector
    h = quiver(c,r,sin(theta(idx)),-cos(theta(idx)),0,'color','white');

end

