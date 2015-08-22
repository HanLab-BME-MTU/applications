function [ h ] = imRadialLine( varargin )
%imRadialLine Create an imline object constrained to be radial from the
%center

dim.x = xlim;
dim.y = ylim;

center = [mean(dim.x) mean(dim.y)];

h = imline(varargin{:});
% setPosition(h,[center; dim.x(2)*3/4 center(2)]);
oldPosition = getPosition(h);

setPositionConstraintFcn(h,@radialConstraint);

    function constrainedPos = radialConstraint(pos)
        pos(pos(:,1) > dim.x(2),1) = dim.x(:,2);
        pos(pos(:,2) > dim.y(2),2) = dim.y(:,2);
        pos(pos(:,1) < dim.x(1),1) = dim.x(:,1);
        pos(pos(:,2) < dim.y(1),2) = dim.y(:,1);
        constrainedPos = pos;
%         constrainedPos(1,:) = center;
        [oldTheta(1),oldRho(1)] = cart2pol(oldPosition(1,1)-center(1),oldPosition(1,2)-center(2));
        [oldTheta(2),oldRho(2)] = cart2pol(oldPosition(2,1)-center(1),oldPosition(2,2)-center(2));
        [theta(1),rho(1)] = cart2pol(pos(1,1)-center(1),pos(1,2)-center(2));
        [theta(2),rho(2)] = cart2pol(pos(2,1)-center(1),pos(2,2)-center(2));
        
        constrainedPt = (rho(1) > rho(2)) + 1;
        referencePt = 3 - constrainedPt;

%         posDiff = pos - oldPosition;
%         ptDiff = sum(posDiff,2) > 0;
%         disp(ptDiff);
        
        newRho = rho(constrainedPt);
        
%         if(all(ptDiff))
%             newRho = rho(referencePt) - abs(oldRho(2)-oldRho(1));
%             newRho = 0;
%         end
        
        [x,y] = pol2cart(theta(referencePt),newRho);
        constrainedPos(constrainedPt,:) = [x y] + center;
        oldPosition = constrainedPos;
        
    end


end

